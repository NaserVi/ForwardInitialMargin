package initialmargin.isdasimm.sensitivity;
import java.util.Arrays;
import java.util.stream.IntStream;

import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.products.AbstractSIMMProduct;
import initialmargin.isdasimm.products.SIMMBermudanSwaption;
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
public class SIMMSensitivityCalculation extends AbstractSIMMSensitivityCalculation{

    private double interpolationStep;
    private LIBORModelMonteCarloSimulationInterface model;
    
    /** Construct a SIMM sensitivity calculation scheme
     * 
     * @param sensitivityMode
     * @param liborWeightMode
     * @param interpolationStep
     * @param model
     * @param isUseTimeGridAdjustment
     * @param isUseAnalyticSwapSensitivities
     * @param isConsiderOISSensitivities
     */
    public SIMMSensitivityCalculation(SensitivityMode sensitivityMode, WeightMode liborWeightMode, 
            double interpolationStep, LIBORModelMonteCarloSimulationInterface model, boolean isUseTimeGridAdjustment, boolean isUseAnalyticSwapSensitivities, boolean isConsiderOISSensitivities){
            
    	super(sensitivityMode, liborWeightMode, isUseTimeGridAdjustment, isUseAnalyticSwapSensitivities,isConsiderOISSensitivities);
    	this.interpolationStep = interpolationStep;
    	this.model = model;
    }
    
    
    /** Construct a SIMM sensitivity calculation scheme
     * 
     * @param sensitivityMode
     * @param liborWeightMode
     * @param interpolationStep
     * @param model
     */
    public SIMMSensitivityCalculation(SensitivityMode sensitivityMode, WeightMode liborWeightMode, double interpolationStep, LIBORModelMonteCarloSimulationInterface model){
            
    	this(sensitivityMode, liborWeightMode, interpolationStep, model, true, true, true);
    	
    }
    
    
	@Override
    public RandomVariableInterface[] getDeltaSensitivitiesIR(AbstractSIMMProduct product,
			 												 String riskClass, 
			 												 String curveIndexName,
			 												 double evaluationTime, 
			 												 LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{

             RandomVariableInterface[] maturityBucketSensis = null;

             switch(sensitivityMode){

                 case Exact:  

                     maturityBucketSensis = doCalculateDeltaSensitivitiesIR(product, curveIndexName, evaluationTime, model);
                     
                     //for(int i=0; i<maturityBucketSensis.length;i++) System.out.println("mbsensis " + maturityBucketSensis[i].getAverage());
                     break;

                 case LinearMelting:
                	 
                	 // The time of the sensitivities used for melting     
        		     double initialMeltingTime = evaluationTime < product.getMeltingResetTime() ? 0 : product.getMeltingResetTime();
                	 
                	 maturityBucketSensis = getMeltedSensitivities(product, null /*given sensitivities*/, initialMeltingTime, evaluationTime, curveIndexName, riskClass);
                	 
                	 //for(int i=0;i<maturityBucketSensis.length;i++) System.out.println("before " + "\t" + evaluationTime + "\t" + maturityBucketSensis[i].getAverage());
                	 
                	 if(product instanceof SIMMBermudanSwaption) {
                		                		 
                		 maturityBucketSensis = ((SIMMBermudanSwaption)product).changeMeltedSensisOnExercisedPaths(evaluationTime, curveIndexName, maturityBucketSensis);
                	 }
                	 
                	 //for(int i=0;i<maturityBucketSensis.length;i++) System.out.println("after " + "\t" + evaluationTime + "\t" + maturityBucketSensis[i].getAverage());
                	 
                	 
                	 break;
                	 
                 case Interpolation:

                	 maturityBucketSensis = getInterpolatedSensitivities(product, riskClass, curveIndexName, evaluationTime, model);

                	 break;

                 default:
                	 break;

             }	

             return maturityBucketSensis;
	}
	
	@Override
	public RandomVariableInterface[] getExactDeltaSensitivitiesIR(AbstractSIMMProduct product, String riskClass, String curveIndexName,
			 													  double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{

			return doCalculateDeltaSensitivitiesIR(product, curveIndexName, evaluationTime, model);

    }
	
	
	/**Linear melting of the sensitivities given on the SIMM Buckets. The melting is perfromed linearly such that
	 * the time zero sensitivities on bucket with maturity N years have vanished after N years. After N years, half of
	 * the sensitivities which are originally on the 2N year bucket will have moved onto the N year bucket.
     * 
	 * @param initialMeltingTime The time at which the melting should start, i.e. time zero
	 * @param evaluationTime The time at which the melted sensitivites are calculated
	 * @param sensitivities The sensitivities on SIMM buckets or LiborPeriodDiscretization to be melted
	 * @param riskClass The SIMM risk class of the product whose sensitivities we consider
	 * @return The melted sensitivities
	 * @throws CalculationException 
	 * @throws CloneNotSupportedException 
	 * @throws SolverException 
     */
	// Add riskType as parameter if vega risk should be considered
	@Override
    public RandomVariableInterface[] getMeltedSensitivities(AbstractSIMMProduct product, 
    														RandomVariableInterface[] sensitivities, double meltingZeroTime,
    														double evaluationTime, String curveIndexName, String riskClass) throws SolverException, CloneNotSupportedException, CalculationException{
		  // also perform melting such that sensitivities are zero at final maturity  
		  
		  if(evaluationTime > product.getFinalMaturity()){
		    	   
		     return mapSensitivitiesOnBuckets(new RandomVariableInterface[]{new RandomVariable(0.0)}, riskClass, new int[]{12/*random*/}, model);
		    	   
		  }
		  
		  // Get sensitivities to melt if not provided as input to the function
		  if(sensitivities == null) {
			  sensitivities = product.getExactDeltaFromCache(meltingZeroTime, riskClass, curveIndexName);		 
		      //for(int i=0; i<sensitivities.length; i++) System.out.println("cacheSensis " + sensitivities[i].getAverage());
		  }
		  
		  double[] initialMeltingTime = new double[]{meltingZeroTime};
		  
	      int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};	
			   
	      // If sensitivities are given on LiborPeriodDiscretization, map them on SIMM Buckets 
		  if(sensitivities.length!=riskFactorsSIMM.length) sensitivities = mapSensitivitiesOnBuckets(sensitivities, riskClass, null, model); 
			   
		  // Get new riskFactor times
		  int[] riskFactorDays = Arrays.stream(riskFactorsSIMM).filter(n -> n > (int)Math.round(365*(evaluationTime-initialMeltingTime[0]))).map(n -> n-(int)Math.round(365*(evaluationTime-initialMeltingTime[0]))).toArray();
		       
		  // Find first bucket later than evaluationTime
		  int firstIndex = IntStream.range(0, riskFactorsSIMM.length)
			                          .filter(i -> riskFactorsSIMM[i]>(int)Math.round(365*(evaluationTime-initialMeltingTime[0]))).findFirst().getAsInt();
			   
		  //Calculate melted sensitivities
		  RandomVariableInterface[] meltedSensis = new RandomVariableInterface[sensitivities.length-firstIndex];
			   
		  for(int i=0;i<meltedSensis.length;i++){
			  meltedSensis[i]=sensitivities[i+firstIndex].mult(1.0-(double)Math.round(365*(evaluationTime-initialMeltingTime[0]))/(double)riskFactorsSIMM[i+firstIndex]);
		  }
			   
		  return mapSensitivitiesOnBuckets(meltedSensis, riskClass, riskFactorDays, null); 
			  
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
    private RandomVariableInterface[] getInterpolatedSensitivities(AbstractSIMMProduct product,
    															   String riskClass, 
			 													   String curveIndexName,
			 													   double evaluationTime, 
			 													   LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{
		 
		 
		 // time of initial and final sensitivities
		 TimeDiscretizationInterface exactSensiTimes = new TimeDiscretization(0,50,interpolationStep);
		 int initialIndex = exactSensiTimes.getTimeIndexNearestLessOrEqual(evaluationTime);
		 double initialTime = (exactSensiTimes.getTime(initialIndex) <= product.getMeltingResetTime()) && (exactSensiTimes.getTime(initialIndex+1) > product.getMeltingResetTime()) ? product.getMeltingResetTime() : exactSensiTimes.getTime(initialIndex);		 
		 double finalTime   = initialTime < product.getMeltingResetTime() ? Math.min(product.getMeltingResetTime(),exactSensiTimes.getTime(initialIndex+1)) : exactSensiTimes.getTime(initialIndex+1);
         
         // Get Sensitivities from exactDeltaCache 
         RandomVariableInterface[] initialSensitivities = product.getExactDeltaFromCache(initialTime, riskClass, curveIndexName);
         RandomVariableInterface[] finalSensitivities   = product.getExactDeltaFromCache(finalTime, riskClass, curveIndexName);

         // Perform linear interpolation
         double deltaT = finalTime-initialTime;
         double deltaTEval = evaluationTime-initialTime;
         
         if(deltaT==0) return finalSensitivities;
         
         //ArrayList<RandomVariableInterface> interpolatedSensis = new ArrayList<>();
         RandomVariableInterface[] interpolatedSensis = new RandomVariable[initialSensitivities.length];   
         for(int bucketIndex=0; bucketIndex<initialSensitivities.length; bucketIndex++){
        	 RandomVariableInterface slope = finalSensitivities[bucketIndex].sub(initialSensitivities[bucketIndex]).div(deltaT);
        	 interpolatedSensis[bucketIndex] = initialSensitivities[bucketIndex].add(slope.mult(deltaTEval));
        	 //interpolatedSensis.add(initialSensitivities[bucketIndex].add(slope.mult(deltaTEval)));
         }
		 //RandomVariableInterface[] result = interpolatedSensis.toArray(new RandomVariableInterface[interpolatedSensis.size()]);
         
		 return interpolatedSensis;
		 
	  }
	  
	  public double getInterpolationStep(){
			return this.interpolationStep;
	  }
	  
	  
}

