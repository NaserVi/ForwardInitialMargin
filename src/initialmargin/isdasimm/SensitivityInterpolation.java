package initialmargin.isdasimm;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

import initialmargin.isdasimm.SIMMPortfolio.PortfolioInstrument;
import initialmargin.isdasimm.SIMMPortfolio.SensitivityMode;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.BermudanSwaption;
import initialmargin.isdasimm.changedfinmath.products.SimpleSwap;
import initialmargin.isdasimm.changedfinmath.products.Swap;
import initialmargin.isdasimm.changedfinmath.products.Swaption;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

/** This class performs linear sensitivity melting on SIMM buckets (possibly with a reset of the sensitivities
 *  to a the true sensitivity values obtained by AAD. Moreover, linear interpolation of sensitivities on the SIMM
 *  buckets may be done with this class.
 * 
 * @author Mario Viehmann
 *
 */

public class SensitivityInterpolation {

	private HashMap<Double /*time*/,HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> sensiMap = new HashMap<Double, HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>>();
	private HashMap<Double /*time*/,HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> bermudanSwapSensiMap = new HashMap<Double, HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>>();

    private PortfolioInstrument product;
    private SIMMPortfolio       portfolio;
    private double              lastEvaluationTime;
    private SensitivityMode     mode;
    private double              finalInterpolationTime;
 
 
	
	public SensitivityInterpolation(PortfolioInstrument product) throws CalculationException{
		this.product   = product;
		this.portfolio = product.portfolio;
		this.mode      = portfolio.getSensitivityMode();
		if(mode == SensitivityMode.Interpolation) this.finalInterpolationTime = getFinalInterpolationTime();
	}
	
	/** Calculate the delta SIMM sensitivities for a given risk class and index curve at a given evaluation time with the specified Libor market model.
	 *  The sensitivities are calculated by interpolation or melting (a particular case of interpolation).
	 * 		
	 * @param riskClass The risk class of the product 
	 * @param curveIndexName The name of the index curve
	 * @param evaluationTime The time at which the sensitivities should be calculated
	 * @param model The Libor market model
	 * @return The sensitivities at evaluationTime on the SIMM buckets (or on Libor Buckets)
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getDeltaSensitivities(String riskClass, 
														   String curveIndexName,
														   double evaluationTime, 
														   LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{
		
		RandomVariableInterface[] maturityBucketSensis = null;
			
		if(evaluationTime > finalInterpolationTime) this.mode = SensitivityMode.LinearMelting; // Cannot interpolate anymore. Use linear Melting instead (only for Bermudan Swaption: Here we interpolate on the not-yet exercised paths)    		    
		
		switch(mode){
		
		   case LinearMelting:
			   
			    HashMap<Double /*time*/,HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> sensiCache = (product.getProduct() instanceof BermudanSwaption && evaluationTime > ((BermudanSwaption)product.getProduct()).getLastValuationExerciseTime().getMin()) ? this.bermudanSwapSensiMap : this.sensiMap;
		        
		        double initialMeltingTime = getInitialMeltingTime(evaluationTime);
		
		        if(isEmptyMap(initialMeltingTime, riskClass, curveIndexName, sensiCache)) setMap(riskClass, curveIndexName, initialMeltingTime, model, sensiCache);
				 		        
		        RandomVariableInterface[] deltaSensis = sensiCache.get(new Double(initialMeltingTime)).get(riskClass).stream().filter(n -> n.containsKey(curveIndexName)).findFirst().get().get(curveIndexName);			        			         
		
		        maturityBucketSensis = getMeltedSensitivities(initialMeltingTime,evaluationTime, deltaSensis, riskClass);
		        
		        if(product.getProduct() instanceof BermudanSwaption && evaluationTime > ((BermudanSwaption)product.getProduct()).getLastValuationExerciseTime().getMin()){
		        	// We melt down the sensitivities on the paths which we have exercised but we continue the interpolation 
		        	// of the sensis on the not-yet-exercised paths between two adjacent exercise times
		        	RandomVariableInterface[] sensisNotExercisedPaths = getSensitivityInterpolation(riskClass, curveIndexName, evaluationTime, model);
		        	for(int i=0; i < maturityBucketSensis.length; i++) maturityBucketSensis[i] = maturityBucketSensis[i].add(sensisNotExercisedPaths[i]);
		        }
		            	
		        break;
		
		   case Interpolation:		   		   
			   
			    maturityBucketSensis = getSensitivityInterpolation(riskClass, curveIndexName, evaluationTime, model);
			    
			    break;
			    
		default:
			    break;
			   			    			 	
		}
		
		
		this.lastEvaluationTime = evaluationTime;		
		
		return maturityBucketSensis;
	}
			
	/** Set the sensitivity map <code> sensiMap <code>: Save for a given time, index curve, risk class and model the sensitivities of the product.
	 *  The sensitivities may be later obtained from the map instead of being re-calculated over and over. 
	 * 
	 * @param riskClass The risk class of the product 
	 * @param curveIndexName The name of the index curve
	 * @param time The time at which the sensitivities should be set
	 * @param model The Libor Market model
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	private void setMap(String riskClass, String curveIndexName, double time, LIBORModelMonteCarloSimulationInterface model, 
			            HashMap<Double /*time*/,HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> cache) throws SolverException, CloneNotSupportedException, CalculationException{
		
		if(cache == this.bermudanSwapSensiMap){
			setBermudanExercisedSwapSensitivityMap(time, product.getProduct(), model);
		} else {
		// Calculate the sensitivities 
		RandomVariableInterface[] deltaSensis = portfolio.doCalculateDeltaSensitivitiesIR(curveIndexName, product, time, model);
		
		// Map the sensitivities onto the SIMM Buckets
		//if(mode == SensitivityMode.LinearMelting || mode == SensitivityMode.Interpolation) deltaSensis = portfolio.getSensitivitiesOnBuckets(deltaSensis, riskClass, null);
		
		// Create a new element of the curveIndex List for given risk class		         
        HashMap<String,RandomVariableInterface[]> curveIndexNameSensiMap = new HashMap<String,RandomVariableInterface[]>();
        curveIndexNameSensiMap.put(curveIndexName,deltaSensis);
        
        // Check if the list of curveIndexNames in the HashMap already exist
        if(sensiMap.containsKey(new Double(time)) && sensiMap.get(new Double(time)).containsKey(riskClass)) {
           sensiMap.get(new Double(time)).get(riskClass).add(curveIndexNameSensiMap);
        }
        else {
        	List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>> list = new ArrayList<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>();
        	list.add(curveIndexNameSensiMap);
        	HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>> riskClassMap = new HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>();
        	riskClassMap.put(riskClass, list);
        	sensiMap.put(new Double(time), riskClassMap);
        }
		}
		
	}
		
	
	 /**Provides the time at which the melting starts. This depends on the product.
	  * 
	  * @param evaluationTime
	  * @return
	  * @throws CalculationException
	  */
	 public double getInitialMeltingTime(double evaluationTime) throws CalculationException{
		   
		   // Determine the initialMeltingTime based on the reset step specified by the SIMMPortfolio
		   //TimeDiscretizationInterface initialMeltingTimes = new TimeDiscretization(0,50,portfolio.getSensiResetStep()); // We need this if we want to reset the melting. Instead, we do interpolation
		   double initialMeltingTime = 0; //initialMeltingTimes.getTime(initialMeltingTimes.getTimeIndexNearestLessOrEqual(evaluationTime));
		   
		   AbstractLIBORMonteCarloProduct product = this.product.getProduct();
		   
		   // Adjust initialMeltingTime and / or clear sensiMap depending on the product
		   if(product instanceof Swaption){
			      double exerciseDate = ((Swaption)product).getExerciseDate();
			      // Reset initial Melting Time before exercise only
			      if(evaluationTime==initialMeltingTime && lastEvaluationTime!=evaluationTime && evaluationTime < exerciseDate) sensiMap.clear();
			      
			      if(((Swaption) product).getDeliveryType()=="Physical" && evaluationTime>=exerciseDate){
			    	  if(lastEvaluationTime < exerciseDate) sensiMap.clear();
			    	  initialMeltingTime = exerciseDate;
			      }
			      			      		   
		   }
		   
		   else if(product instanceof SimpleSwap || product instanceof Swap){
			   double forwardStartTime = 0; // Swap always starting at time zero
			   if(product instanceof SimpleSwap) forwardStartTime = ((SimpleSwap)product).getStartTime();		   
			   if(evaluationTime==initialMeltingTime && lastEvaluationTime!=evaluationTime && evaluationTime < forwardStartTime) sensiMap.clear();
			      
			   if(evaluationTime>=forwardStartTime){			    
			    	  initialMeltingTime = forwardStartTime;
			    	  if(lastEvaluationTime < forwardStartTime) sensiMap.clear();
			   }
			   
		   }
		   
		   else if(product instanceof BermudanSwaption){
			   double firstExerciseTime = ((BermudanSwaption)product).getLastValuationExerciseTime().getMin();
			   if(evaluationTime==initialMeltingTime && lastEvaluationTime!=evaluationTime && evaluationTime <firstExerciseTime) sensiMap.clear(); // Melting reset
		       if(evaluationTime >= firstExerciseTime) {
		    	  double[] exerciseTimes = ((BermudanSwaption)product).getExerciseTimes();
		 		  int      lastExerciseIndex = new TimeDiscretization(exerciseTimes).getTimeIndexNearestLessOrEqual(evaluationTime);
		 		           initialMeltingTime = new TimeDiscretization(exerciseTimes).getTime(lastExerciseIndex);

		       }
		   }
		   
		   return initialMeltingTime;
	  }
	  
	  
	  /** Set the sensitivity map containing the sensitivities of the swaps on the paths on which we have exercised the Bermudan Swaption. 
	   * 
	   * @param evaluationTime The time of evaluation
	   * @param product The product (Bermudan Swaption)
	   * @param model The Libor market model
	   * @return The initial melting time
	   * @throws CalculationException
	   * @throws CloneNotSupportedException 
	   * @throws SolverException 
	   */
	  private void setBermudanExercisedSwapSensitivityMap(double evaluationTime, AbstractLIBORMonteCarloProduct product, LIBORModelMonteCarloSimulationInterface model) throws CalculationException, SolverException, CloneNotSupportedException{
		  double[] exerciseTimes = ((BermudanSwaption)product).getExerciseTimes();
		  double firstExerciseTime = ((BermudanSwaption)product).getLastValuationExerciseTime().getMin();
		  int lastExerciseIndex = new TimeDiscretization(exerciseTimes).getTimeIndexNearestLessOrEqual(evaluationTime);
		  double lastExerciseTime = new TimeDiscretization(exerciseTimes).getTime(lastExerciseIndex);
		  double previousExerciseTime = lastExerciseIndex==0 ? lastExerciseTime : new TimeDiscretization(exerciseTimes).getTime(lastExerciseIndex-1);

		  RandomVariableInterface[] meltedSensisOIS = null;
		  RandomVariableInterface[] meltedSensisLibor = null;
		  
		  // Reset Melting Maps if evaluationTime for the first time exceeds the exercise times of the Bermudan 
		  if((lastEvaluationTime == firstExerciseTime) || (evaluationTime>=lastExerciseTime && lastEvaluationTime < lastExerciseTime && lastExerciseTime!=firstExerciseTime))  {		
			 
			 double evalTimeAdjustment = lastEvaluationTime == firstExerciseTime ? (evaluationTime-lastEvaluationTime) : 0;
			 // Create a Swap 
			 double[] fixingDates  = ((BermudanSwaption)product).getFixingDates(evaluationTime);
			 double[] paymentDates = IntStream.range(0, fixingDates.length).mapToDouble(i->fixingDates[i]+0.5).toArray();
			 double[] swapRates    = ((BermudanSwaption)product).getSwapRates();
			 double[] periodNotional = ((BermudanSwaption)product).getPeriodNotionals();
			
			 SimpleSwap swap = new SimpleSwap(fixingDates,paymentDates,swapRates,periodNotional);
			 SIMMClassifiedProduct swapClassified = new SIMMClassifiedProduct(swap,"RatesFX",new String[] {"InterestRate"}, new String[] {"OIS","Libor6m"},"EUR",null,false, false);
		     PortfolioInstrument swapCache = portfolio.createInstrumentCache(swapClassified);
			 
		     // Calculate sensitivities of the swap dV/dS analytically. The sensitivities are already mapped onto SIMM buckets
			 RandomVariableInterface[] swapSensisLibor = portfolio.doCalculateDeltaSensitivitiesIR("Libor6m",swapCache,evaluationTime-evalTimeAdjustment, model);//(evaluationTime, fixingDates, null, periodLength[0], portfolio.getModel(),"Libor");
		     RandomVariableInterface[] swapSensisOIS   = portfolio.doCalculateDeltaSensitivitiesIR("OIS",swapCache,evaluationTime-evalTimeAdjustment, model);
		     
		     if(lastExerciseTime==firstExerciseTime) { // The melting map is created at the first exercise time. We need to clear it.
				 bermudanSwapSensiMap.clear();
				 meltedSensisLibor = new RandomVariableInterface[swapSensisLibor.length];
				 meltedSensisOIS   = new RandomVariableInterface[swapSensisLibor.length];
			 } else {
				 RandomVariableInterface[] lastExactSensisOIS   = bermudanSwapSensiMap.get(new Double(previousExerciseTime)).get("InterestRate").stream().filter(n -> n.containsKey("OIS")).findFirst().get().get("OIS");
				 RandomVariableInterface[] lastExactSensisLibor = bermudanSwapSensiMap.get(new Double(previousExerciseTime)).get("InterestRate").stream().filter(n -> n.containsKey("Libor6m")).findFirst().get().get("Libor6m");
			     // Melt sensis to evaluationDate
				 meltedSensisLibor = getMeltedSensitivities(exerciseTimes[lastExerciseIndex-1],evaluationTime, lastExactSensisLibor, "InterestRate");
				 meltedSensisOIS   = getMeltedSensitivities(exerciseTimes[lastExerciseIndex-1],evaluationTime, lastExactSensisOIS, "InterestRate");
			 }
		     
		     // We perform the melting of the Swaps only on the paths at which we have exercised
		     RandomVariableInterface pathExerciseTimes = ((BermudanSwaption)product).getLastValuationExerciseTime();
			 for(int i=0;i<swapSensisLibor.length;i++) { 
				 // Get sensis only on paths on which we have exercised
				 swapSensisLibor[i] = swapSensisLibor[i].barrier(new RandomVariable(pathExerciseTimes.sub(evaluationTime+0.0001)), new RandomVariable(0.0), swapSensisLibor[i]);
				 swapSensisOIS[i] = swapSensisOIS[i].barrier(new RandomVariable(pathExerciseTimes.sub(evaluationTime+0.0001)), new RandomVariable(0.0), swapSensisOIS[i]);
				 if(lastExerciseTime!=firstExerciseTime) {
					 swapSensisLibor[i] = swapSensisLibor[i].barrier(new RandomVariable(pathExerciseTimes.sub(previousExerciseTime+0.0001)), swapSensisLibor[i], new RandomVariable(0.0));
					 swapSensisOIS[i] = swapSensisOIS[i].barrier(new RandomVariable(pathExerciseTimes.sub(previousExerciseTime+0.0001)), swapSensisOIS[i], new RandomVariable(0.0));
				 }
				 // 
				 meltedSensisLibor[i] = lastExerciseTime==firstExerciseTime ? swapSensisLibor[i] : meltedSensisLibor[i].add(swapSensisLibor[i]);
				 meltedSensisOIS[i] = lastExerciseTime==firstExerciseTime ? swapSensisOIS[i] : meltedSensisOIS[i].add(swapSensisOIS[i]);
			 }
			 
			 // Put sensitivities to the bermudanSwapSensiMap
			 // Create a new element of the curveIndex List for given risk class		         
	         HashMap<String,RandomVariableInterface[]> curveIndexNameSensiMap = new HashMap<String,RandomVariableInterface[]>();
	         curveIndexNameSensiMap.put("Libor6m",meltedSensisLibor);
	         curveIndexNameSensiMap.put("OIS",meltedSensisOIS);
	         
	         List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>> list = new ArrayList<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>();
	         list.add(curveIndexNameSensiMap);
	         HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>> riskClassMap = new HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>();
	         riskClassMap.put("InterestRate", list);
	         bermudanSwapSensiMap.put(new Double(lastExerciseTime), riskClassMap); // lastExerciseTime == initialMeltingTime
	         
		  }
	  }
	  
	
	  /**Linear melting of the sensitivities given on the SIMM Buckets or - not mapped to SIMM Bucket i.e. on the LiborPeriodDiscretization
		 * 
		 * @param initialMeltingTime The time at which the melting should start, i.e. time zero
		 * @param evaluationTime The time at which the melted sensitivites are calculated
		 * @param sensitivities The sensitivities on SIMM buckets or LiborPeriodDiscretization to be melted
		 * @param riskClass The SIMM risk class of the product whose sensitivities we consider
		 * @return The melted sensitivities
		 */
	  private RandomVariableInterface[] getMeltedSensitivities(double initialMeltingTime, double evaluationTime, RandomVariableInterface[] sensitivities, String riskClass){
		  
		  RandomVariableInterface[] sensisOnBuckets = null;
		  
		  switch(mode){
		    case LinearMelting:
		    	
		       if(this.product.getProduct() instanceof Swaption && ((Swaption)product.getProduct()).getDeliveryType()!="Physical"){
		    	   sensisOnBuckets = portfolio.getSensitivitiesOnBuckets(new RandomVariableInterface[]{new RandomVariable(0.0)}, riskClass, new int[]{12/*random*/});
		    	   break;
		       }
		       
			   int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};	
			   
			   // If sensitivities are given on LiborPeriodDiscretization, map them on SIMM Buckets 
			   if(sensitivities.length!=riskFactorsSIMM.length) sensitivities = portfolio.getSensitivitiesOnBuckets(sensitivities, riskClass, null); 
			   
			   // Get new riskFactor times
			   int[] riskFactorDays = Arrays.stream(riskFactorsSIMM).filter(n -> n > (int)Math.round(365*(evaluationTime-initialMeltingTime))).map(n -> n-(int)Math.round(365*(evaluationTime-initialMeltingTime))).toArray();
		       
			   // Find first bucket later than evaluationTime
			   int firstIndex = IntStream.range(0, riskFactorsSIMM.length)
			                          .filter(i -> riskFactorsSIMM[i]>(int)Math.round(365*(evaluationTime-initialMeltingTime))).findFirst().getAsInt();
			   
			   //Calculate melted sensitivities
			   RandomVariableInterface[] meltedSensis = new RandomVariableInterface[sensitivities.length-firstIndex];
			   
			   for(int i=0;i<meltedSensis.length;i++){
				   meltedSensis[i]=sensitivities[i+firstIndex].mult(1.0-(double)Math.round(365*(evaluationTime-initialMeltingTime))/(double)riskFactorsSIMM[i+firstIndex]);
			   }
			   
			   sensisOnBuckets = portfolio.getSensitivitiesOnBuckets(meltedSensis, riskClass, riskFactorDays); 
			   break;
		       
		    case LinearMeltingOnLiborBuckets:
		
			   TimeDiscretizationInterface times = portfolio.getModel().getLiborPeriodDiscretization();
		       
			   // Find first bucket later than evaluationTime
			   firstIndex = evaluationTime-initialMeltingTime==0 ? 1 : times.getTimeIndexNearestGreaterOrEqual(evaluationTime-initialMeltingTime);

			   //Calculate melted sensitivities
			   meltedSensis = new RandomVariableInterface[sensitivities.length-firstIndex+1];
			   
			   for(int i=0;i<meltedSensis.length;i++){
				  double time = (i+firstIndex)*portfolio.getModel().getLiborPeriodDiscretization().getTimeStep(0);
				  meltedSensis[i]=sensitivities[i+firstIndex-1].mult(1.0-(evaluationTime-initialMeltingTime)/time);
			   }
			   
			   sensisOnBuckets = portfolio.getSensitivitiesOnBuckets(meltedSensis, riskClass, null);  
			   break;
			   
		    default:
				  break;
		    }
		  
			
			return sensisOnBuckets;
		}
	  
	  
	  /**Gives the maximal time until one can use interpolation for the product
	   * 
	   * @return The maximal time of the interpolation
	   * @throws CalculationException
	   */
	  public double getFinalInterpolationTime() throws CalculationException{
		   
		  double finalInterpolationTime = 87^3;
		   		 		   
		  AbstractLIBORMonteCarloProduct product = this.product.getProduct();
		  
		  if(product instanceof Swap) finalInterpolationTime = 0;
		   
		  if(product instanceof Swaption) finalInterpolationTime = ((Swaption)product).getExerciseDate();
			      
		  if(product instanceof BermudanSwaption) finalInterpolationTime = 100;//((BermudanSwaption)product).getLastValuationExerciseTime().getMin();
			
		  return finalInterpolationTime;
	  }
	   
	   
	  /** Interpolates sensitivities on SIMM buckets linearly between two exact sensitivities obtained by AAD. 
	   *  Information of future sensitivities (after evaluation time) is used. 
	   * 
	   * @param riskClass The risk class of the product 
	   * @param curveIndexName The name of the index curve 
	   * @param evaluationTime The time of evaluation
	   * @param model The Libor market model
	   * @return The interpolated sensitivities on SIMM buckets at evaluation time
	   * @throws SolverException
	   * @throws CloneNotSupportedException
	   * @throws CalculationException
	   */
	  public RandomVariableInterface[] getSensitivityInterpolation(String riskClass, 
			 													   String curveIndexName,
			 													   double evaluationTime, 
			 													   LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{
		 
		 double initialTime;
		 double finalTime;
		 AbstractLIBORMonteCarloProduct product = this.product.getProduct();
		 RandomVariableInterface exerciseIndicator = new RandomVariable(1.0);
		 
		 if(product instanceof BermudanSwaption && evaluationTime >= finalInterpolationTime){
			 double[] exerciseTimes = ((BermudanSwaption)product).getExerciseTimes();
			 double   finalMaturity = ((BermudanSwaption)product).getFinalMaturity();
			 int initialIndex =  new TimeDiscretization(exerciseTimes).getTimeIndexNearestLessOrEqual(evaluationTime);
			 
			 // return zero after final maturity
			 if(evaluationTime >=finalMaturity) return portfolio.getSensitivitiesOnBuckets(new RandomVariableInterface[]{new RandomVariable(0.0)},riskClass,new int[]{1});
			 
			 initialTime = new TimeDiscretization(exerciseTimes).getTime(initialIndex);
			 // final time cannot be later than the final maturity
			 finalTime =  initialIndex==exerciseTimes.length-1 ? finalMaturity : new TimeDiscretization(exerciseTimes).getTime(initialIndex+1); 
			 // no interpolation at first exercise time ( = "finalInterpolationTime"). We still use the sensis on all paths.
			 if(evaluationTime==finalInterpolationTime) finalTime=initialTime;
			 // After the first exercise time, we interpolate only on the paths on which we have not yet exercised. On the other paths we melt the sensis to zero.
			 if(evaluationTime > finalInterpolationTime) {
				 RandomVariableInterface pathExerciseTimes = ((BermudanSwaption)product).getLastValuationExerciseTime();
				 // exerciseIndicator: 1 if not yet exercised, 0 if already exercised at evaluationTime
				 exerciseIndicator = new RandomVariable(1.0).barrier(new RandomVariable(pathExerciseTimes.sub(evaluationTime+0.0001)), new RandomVariable(1.0), new RandomVariable(0.0));			 
			 }
						 
		 } else {
		     // time of initial and final sensitivities
			 TimeDiscretizationInterface exactSensiTimes = new TimeDiscretization(0,50,portfolio.getSensiResetStep());
			 int initialIndex = exactSensiTimes.getTimeIndexNearestLessOrEqual(evaluationTime);
			 initialTime = exactSensiTimes.getTime(initialIndex); // always smaller than finalInterpolationTime
			 finalTime   = Math.min(getFinalInterpolationTime(),exactSensiTimes.getTime(initialIndex+1));
		 }
		 
         // Set sensiMap for intial and final time if applicable
         if(!sensiMap.containsKey(new Double(initialTime)) || !sensiMap.get(new Double(initialTime)).containsKey(riskClass)){
        	for(String curveName : this.product.getClassifiedProduct().getCurveIndexNames()) setMap(riskClass, curveName, initialTime, model, sensiMap);
         }
         if(!sensiMap.containsKey(new Double(finalTime)) || !sensiMap.get(new Double(finalTime)).containsKey(riskClass)){
         	for(String curveName : this.product.getClassifiedProduct().getCurveIndexNames()) setMap(riskClass, curveName, finalTime, model, sensiMap);
         }
         
         // Get Sensitivities from sensiMap
         RandomVariableInterface[] initialSensitivities = sensiMap.get(new Double(initialTime)).get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().get().get(curveIndexName);
         RandomVariableInterface[] finalSensitivities   = sensiMap.get(new Double(finalTime)).get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().get().get(curveIndexName);

         // Multiply sensitivities with exercise indicator
         RandomVariableInterface indicator = new RandomVariable(exerciseIndicator);
         initialSensitivities = Arrays.stream(initialSensitivities).map(n -> n.mult(indicator)).toArray(RandomVariableInterface[]::new);
         finalSensitivities   = Arrays.stream(finalSensitivities).map(n -> n.mult(indicator)).toArray(RandomVariableInterface[]::new);
         
         // Perform linear interpolation
         double deltaT = finalTime-initialTime;
         double deltaTEval = evaluationTime-initialTime;
         
         if(deltaT==0) return finalSensitivities;
         
         RandomVariableInterface[] interpolatedSensis = new RandomVariable[finalSensitivities.length];
                  
         for(int bucketIndex=0; bucketIndex<interpolatedSensis.length; bucketIndex++){
        	 RandomVariableInterface slope = finalSensitivities[bucketIndex].sub(initialSensitivities[bucketIndex]).div(deltaT);
        	 interpolatedSensis[bucketIndex] = initialSensitivities[bucketIndex].add(slope.mult(deltaTEval));
         }
		 
         return interpolatedSensis;
		 
	  }

	  /** Checks if the sensitivity map of the class <code> SensitivityInterpolation <code> is empty
	   * 
	   * @param time  The time of the sensitivity
	   * @param riskClass The riskClass corresponding to the considered sensitivity
	   * @param curveIndexName The curve index name corresponding to the considered sensitivity (e.g. OIS, Libor6m)
	   * @return <code> true <code> if the map is empty
	   */
	  public boolean isEmptyMap(double time, String riskClass, String curveIndexName, HashMap<Double /*time*/,HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> sensiMap){
	      return (!sensiMap.containsKey(new Double(time))  ||  !sensiMap.get(new Double(time)).containsKey(riskClass) || !sensiMap.get(new Double(time)).get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().isPresent()); 
	  }
}
