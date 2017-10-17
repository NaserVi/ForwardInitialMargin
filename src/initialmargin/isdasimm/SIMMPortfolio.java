package initialmargin.isdasimm;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import initialmargin.simm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import net.finmath.analytic.model.curves.DiscountCurve;
import net.finmath.exception.CalculationException;
//import net.finmath.marketdata.model.curves.DiscountCurve;
//import net.finmath.montecarlo.AbstractMonteCarloProduct;
import net.finmath.montecarlo.MonteCarloSimulationInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
//import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
//import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

//used to have LIBOR Market Model with analytic discount curve
import initialmargin.simm.changedfinmath.*;
import initialmargin.simm.changedfinmath.products.*;

public class SIMMPortfolio extends AbstractLIBORMonteCarloProduct{
	
	// Class variables of SIMMProductCollection
	private PortfolioInstrument[] portfolioProducts;
    private RandomVariableInterface[][] riskWeightToLiborAdjustments;
    private boolean isUseTimeGridAdjustment=true; // can be discarded later.. just to check how discount curve influences IM
	private double[] discountCurvePillars = {0.5 , 1.0, 2.0, 5.0, 30.0};// shouldn't the discount curve know its pillars ?
    private ConditionalExpectationEstimatorInterface conditionalExpectationOperator;// conditionalExpectationOperator;
    
    public enum WeightToLiborAdjustmentMethod{
		Constant,  //Sets dL/dS(t=0) for all forward IM times, i.e. leave the weight adjustment dL/dS constant
		Stochastic //Calculate dL/dS(t) for all forward IM times, i.e. (weakly) stochastic weight adjustment 
	}
    
    private WeightToLiborAdjustmentMethod liborWeightMethod;
    private LIBORModelMonteCarloSimulationInterface model = null;
    private SIMMSchemeMain SIMMScheme;
    private String calculationCCY; 
    
    public enum SensitivityMode{
    	LinearMelting,
    	Stochastic
    }
    private SensitivityMode sensitivityMode;
    
	public class PortfolioInstrument {
		private SIMMClassifiedProduct classifiedProduct;
	    private Map<Long, RandomVariableInterface> gradientOfProduct; // Same for all evaluationTimes; Is reset for different products
	    private double lastEvaluationTime;
	    
	    final private String[]  CreditMaturityBuckets = {"1y","2y","3y","5y","10y"};
        final private String[]  IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};

	    private HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,
	                HashMap<String/*maturityBucket*/,RandomVariableInterface>>>> deltaSensitivities = new HashMap<String,List<HashMap<String,HashMap<String,RandomVariableInterface>>>>(); // currently only for InterestRate riskClass
	    private RandomVariableInterface vegaSensitivity=null;
	    private HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>> meltingMap = new HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>();
	    private PortfolioInstrument(SIMMClassifiedProduct product){
		   this.classifiedProduct=product;
	    }
	
	    /**Calculates the sensitivity of this product if the requested sensitivity is available. 
	     * 
	     * @param productClass The SIMM product class of this product (RatesFx etc.)
	     * @param riskClass The SIMM risk class of this product (InterestRate etc.)
	     * @param maturityBucket The SIMM maturity bucket for the case of Interest or Credit risk class
	     * @param curveIndexNames The name of the relevant curves for this product (OIS, Libor6m etc.)
	     * @param bucketKey The SIMM bucket key of this product (it is the currency for InterestRate risk class)
	     * @param riskType The risk type: delta, vega or curvature
	     * @param evaluationTime The evaluation time
	     */
	    public RandomVariableInterface getSensitivity(String productClass, 
			                                      String riskClass, 
			                                      String maturityBucket, // only for IR and Credit risk class, null otherwise
			                                      String curveIndexName, // null if riskClass not IR
			                                      String bucketKey,      // currency for IR otherwise bucket nr.
			                                      String riskType, double evaluationTime) throws CalculationException, SolverException, CloneNotSupportedException{

		   RandomVariableInterface result=null;
		   if(!classifiedProduct.getHasOptionality() && riskType!="delta") return new RandomVariable(0.0);
		   if(evaluationTime!=lastEvaluationTime) {
			   clear();
			   this.lastEvaluationTime = evaluationTime;
		   }
		   if(gradientOfProduct==null) setGradient(model); // needs to be set only once
		
		   if(productClass==classifiedProduct.getProductClass() && Arrays.asList(classifiedProduct.getRiskClasses()).contains(riskClass)){
		   
		   switch(riskType){
		      case("delta"): 
			      switch(riskClass){
			          case("InterestRate"):
			 
			              if(Arrays.asList(classifiedProduct.getCurveIndexNames()).contains(curveIndexName) & bucketKey==classifiedProduct.getCurrency()){
				             // There exists a sensitivity. Check if the sensitivities (on all maturityBuckets) have already been calculated for given riskClass and riskType)
			            	  
			            	  if(!deltaSensitivities.containsKey(riskClass) || !deltaSensitivities.get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().isPresent()){
					            // The sensitivities need to be calculated for the given riskClass and riskType
			            		RandomVariableInterface[] maturityBucketSensis;
			            		if(sensitivityMode == SensitivityMode.LinearMelting & (!meltingMap.containsKey(riskClass) || !meltingMap.get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().isPresent())) {
			            			maturityBucketSensis = doCalculateDeltaSensitivitiesIR(curveIndexName, this, 0.0);			            			
			            			// Create a new element of the curveIndex List for given risk class		         
						            HashMap<String,RandomVariableInterface[]> curveIndexNameSensiMap = new HashMap<String,RandomVariableInterface[]>();
						            curveIndexNameSensiMap.put(curveIndexName,maturityBucketSensis);
						            // Check if list already exist
						            if(meltingMap.containsKey(riskClass)){
						            	meltingMap.get(riskClass).add(curveIndexNameSensiMap);
						            } else {
						            	List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>> list = new ArrayList<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>();
						            	list.add(curveIndexNameSensiMap);
						            	meltingMap.put(riskClass, list);
						            }	
			            		}  
			            		if(sensitivityMode == SensitivityMode.LinearMelting){
			            			maturityBucketSensis = meltingMap.get(riskClass).stream().filter(n -> n.containsKey(curveIndexName)).findFirst().get().get(curveIndexName);
			            		    maturityBucketSensis = getMeltedSensitivities(evaluationTime, maturityBucketSensis, riskClass);
			            		} else {
					                maturityBucketSensis = doCalculateDeltaSensitivitiesIR(curveIndexName, this, evaluationTime); // currently only for riskClass IR
			            		}
					            // Create a new element of the curveIndex List for given risk class		         
					            HashMap<String,HashMap<String,RandomVariableInterface>> curveIndexNameSensiMap = new HashMap<String,HashMap<String,RandomVariableInterface>>();
					            HashMap<String,RandomVariableInterface> bucketSensitivities = new HashMap<String,RandomVariableInterface>();
					           
					            for(int i=0;i<IRMaturityBuckets.length;i++) bucketSensitivities.put(IRMaturityBuckets[i], maturityBucketSensis[i]);
					            curveIndexNameSensiMap.put(curveIndexName,bucketSensitivities);
					            // Check if list already exist
					            if(deltaSensitivities.containsKey(riskClass)){
					            	deltaSensitivities.get(riskClass).add(curveIndexNameSensiMap);
					            } else {
					            	List<HashMap<String/*curveIndexName*/,HashMap<String/*maturityBucket*/,RandomVariableInterface>>> list = new ArrayList<HashMap<String/*curveIndexName*/,HashMap<String/*maturityBucket*/,RandomVariableInterface>>>();
					            	list.add(curveIndexNameSensiMap);
					            	deltaSensitivities.put(riskClass, list);
					            }					    					        					          
				             }
				             result = deltaSensitivities.get(riskClass).stream().filter(n -> n.containsKey(curveIndexName)).findFirst().get().get(curveIndexName).get(maturityBucket);			                  
			              } else result = new RandomVariable(0.0); // There exists no delta Sensi for risk Class InterestRate
			           break;
		               case("CreditQ"):
		               case("CreditNonQ"):
		               case("FX"):
		               case("Commodity"):
		               case("Equity"): result = null;
			       } break;
		      
		      case("vega"): 
			  case("curvature"):
			       switch(riskClass){
			           case("InterestRate"): //if(vegaSensitivity!=null) vegaSensitivity = getVegaSensitivityIR(curveIndexNames, product, evaluationTime, model);
			           case("CreditQ"):
			           case("CreditNonQ"):
			           case("FX"):
			           case("Commodity"):
			           case("Equity"): result=null;
			       }
		    }
		 }
	    return result;
	  } // end getSensitivity()
	
	  private void setGradient(LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		  RandomVariableDifferentiableInterface productValue = (RandomVariableDifferentiableInterface) getProduct().getValue(model.getTime(0), model);
		  Map<Long, RandomVariableInterface> gradientOfProduct = productValue.getGradient();
		  this.gradientOfProduct = gradientOfProduct;
	    }
	   
	  public SIMMClassifiedProduct getClassifiedProduct(){
		  return this.classifiedProduct;
	    }
	  
	  public AbstractLIBORMonteCarloProduct getProduct(){
	    	return getClassifiedProduct().getProduct();
	    }
	
	  public void clear(){
		  if(this.deltaSensitivities!=null) this.deltaSensitivities.clear();
		  this.vegaSensitivity = null;
		  
	    }

	  public Map<Long, RandomVariableInterface> getGradient() {
		   return this.gradientOfProduct;
	    }
  }// end class PortfolioInstrument
	
	/**Construct a <code> SIMMPortfolio </code> 
	 * 
	 * @param classifiedProducts The portfolio products 
	 * @param calculationCurrency The calculation currency
	 * @param method The libor weight adjustment method
	 */
	public SIMMPortfolio(SIMMClassifiedProduct[] classifiedProducts,
			             String calculationCurrency,
			             SensitivityMode sensiMode,
			             WeightToLiborAdjustmentMethod method){

		  this.portfolioProducts = createPortfolioInstruments(classifiedProducts);
		  this.calculationCCY = calculationCurrency;
		  this.sensitivityMode = sensiMode;
		  this.liborWeightMethod = method;
	}
	
	@Override
	public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		if(this.model==null || !model.equals(this.model)) { // At inception (t=0) or if the model is reset
			this.model = model; //...the (new) model must be set
			clearProducts();    //...the variables of classifiedSIMMProducts are reset to null (they have to be recalculated under the new model)
	        if(liborWeightMethod == WeightToLiborAdjustmentMethod.Constant){
	     			this.riskWeightToLiborAdjustments = getLiborSwapSensitivities(0.0 /*evaluationTime*/);
	        } 
			this.SIMMScheme= new SIMMSchemeMain(this,this.calculationCCY);
		}
		setConditionalExpectationOperator(evaluationTime);
		return SIMMScheme.getValue(evaluationTime);
	}
	

	/**Calculate the sensitivities dV/dS with respect to all swap rates for given product and curve. This applies to the risk class Interest Rates only.
	 * 
	 * @param curveIndexName The name of the curve to be considered (OIS, LiborXm)
	 * @param product The product whose sensitivity is to be considered.
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return The sensitivities dV/dS i.e. with respect to swap rates.
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	private RandomVariableInterface[] doCalculateDeltaSensitivitiesIR(String curveIndexName, // include inflation risk, ccybasis and other tenors 
                                                                      PortfolioInstrument product, 
                                                                      double evaluationTime) throws SolverException, CloneNotSupportedException, CalculationException{
		RandomVariableInterface[] delta;
		if(curveIndexName!="OIS"){ 
		   RandomVariableInterface[] dVdL = getValueLiborSensitivities(product, evaluationTime);
		   // the following line will be removed later. Just checking how timeGridAdjustment affects the result
		   int timeGridIndicator = 0; if(!isUseTimeGridAdjustment && !onLiborPeriodDiscretization(evaluationTime)) timeGridIndicator = 1;
		
		   delta = new RandomVariableInterface[dVdL.length-timeGridIndicator];
		   RandomVariableInterface[][] dLdS;
		   if(this.liborWeightMethod == WeightToLiborAdjustmentMethod.Stochastic){
			   dLdS = getLiborSwapSensitivities(evaluationTime);
		   } else dLdS = this.riskWeightToLiborAdjustments;
		   // Calculate Sensitivities wrt Swaps
		   for(int swapIndex = 0; swapIndex<dVdL.length-timeGridIndicator; swapIndex++){
			   RandomVariableInterface dVdS  =new RandomVariable(0.0);
			   RandomVariableInterface factor;
			   for(int liborIndex=0;liborIndex<dVdL.length-timeGridIndicator;liborIndex++){
				      factor = dLdS[liborIndex][swapIndex]==null ?  new RandomVariable(0.0) : dLdS[liborIndex][swapIndex];
					  dVdS = dVdS.addProduct(dVdL[liborIndex+timeGridIndicator], factor);
		       }
			   delta[swapIndex]=dVdS;
		   }
		 } else { // CurveIndexName == OIS
			delta = getDiscountCurveSensitivities(product,evaluationTime);
		 }
		 return getSensitivitiesOnBuckets(delta, "InterestRate", null);
	}
	
	/**Performs rebucketing of sensitivities to the SIMM buckets by linear interpolation.
	 * 
	 * @param sensitivities The sensitivities wrt swap rates dV/dS
	 * @param riskClass The risk class
	 * @return The sensitivities on the SIMM maturity buckets
	 */
	private RandomVariableInterface[] getSensitivitiesOnBuckets(RandomVariableInterface[] sensitivities, String riskClass, int[] riskFactorDays){
		//rebucketing to SIMM structure(buckets: 2w, 1m, 3m, 6m, 1y, 2y, 3y, 5y, 10y, 15y, 20y, 30y)	
		int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};	
		RandomVariableInterface[] deltaSIMM = new RandomVariableInterface[riskFactorsSIMM.length];
		for(int i = 0;i<deltaSIMM.length;i++) deltaSIMM[i] = new RandomVariable(0.0);
		if(riskFactorDays==null){// in case of sensitivities dV/dS at each time of the LiborPeriodDiscretization
		   riskFactorDays = new int[sensitivities.length];
		   // act/365 as default daycount convention
		   for(int i=0;i<sensitivities.length;i++) riskFactorDays[i] = (int)Math.round(365 * model.getLiborPeriodDiscretization().getTime(i+1));	
		}
		int counter = 0;
		for(int simmFactor =0; simmFactor<riskFactorsSIMM.length;simmFactor++){
			for(int i = counter; i<sensitivities.length; i++){
				
							    
					if(riskFactorDays[i] < riskFactorsSIMM[0]){
						deltaSIMM[0] = deltaSIMM[0].add(sensitivities[i]);
						counter++;
					}
					else{
						if(riskFactorDays[i] >= riskFactorsSIMM[riskFactorsSIMM.length-1]){
							deltaSIMM[deltaSIMM.length-1] = deltaSIMM[deltaSIMM.length-1].add(sensitivities[i]);
						}
					
						else{
							if(riskFactorDays[i] >= riskFactorsSIMM[simmFactor] && riskFactorDays[i] < riskFactorsSIMM[simmFactor+1]){
					
							deltaSIMM[simmFactor] = deltaSIMM[simmFactor].addProduct(sensitivities[i],((double)(riskFactorsSIMM[simmFactor+1] - riskFactorDays[i]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							deltaSIMM[simmFactor+1] = deltaSIMM[simmFactor+1].addProduct(sensitivities[i],((double)(riskFactorDays[i]-riskFactorsSIMM[simmFactor]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							counter++;
							}							
							else{
							break;
							}
						}
					
					}
			}
			
		}
		
	return deltaSIMM;		
			
	}
	
	
	/**Calculates dL/dS 
	 * 
	 * @param evaluationTime The time at which the sensitivity is calculated
	 * @return The matrix dL/dS 
	 * @throws CalculationException
	 */
	public RandomVariableInterface[][] getLiborSwapSensitivities(double evaluationTime) throws CalculationException{
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		RandomVariableInterface[][] dLdS=null;
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		
	    // Get index of first Libor starting >= evaluationTime
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int numberOfRemainingLibors = model.getNumberOfLibors()-nextLiborIndex;
		dLdS = new RandomVariableInterface [numberOfRemainingLibors][numberOfRemainingLibors];
					
		// Calculate dLdS directly  
		dLdS[0][0]=new RandomVariable(1.0);
		double discountTime = evaluationTime+liborPeriodLength;
		RandomVariableInterface sumDf = model.getNumeraire(discountTime).invert();
		for(int liborIndex = 1; liborIndex<dLdS.length;liborIndex++){
		    discountTime +=model.getLiborPeriodDiscretization().getTimeStep(0);
		    RandomVariableInterface df = model.getNumeraire(discountTime).invert();
		    RandomVariableInterface denominator = df.getConditionalExpectation(conditionalExpectationOperator);
		    dLdS[liborIndex][liborIndex-1]=sumDf.getConditionalExpectation(conditionalExpectationOperator).div(denominator).mult(-1.0);//dLdS[liborIndex][liborIndex-1]=-sumDf.getConditionalExpectation(conditionalExpectationOperator).getAverage()/denominator;
		    sumDf = sumDf.add(df);
		    dLdS[liborIndex][liborIndex] = sumDf.getConditionalExpectation(conditionalExpectationOperator).div(denominator);
		}
		
		return dLdS;
	}
	
	
	
	/**Since dV/dL is wrt the incorrect Libor times this function provides a matrix dL/dL to be multiplied with dV/dL in order to 
	 * have the correct libor times starting at evaluationTime. 
	 * @param evaluationTime The time at which the adjustment should be calculated.
	 * @return Pseudo Inverse of derivative band matrix; Identity matrix in case of evaluationTime on LiborPeriodDiscretization; 
	 * @throws CalculationException
	 */
	private RandomVariableInterface[][] getLiborTimeGridAdjustment(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime);
		
		// If evaluationTime lies on Libor Time Grid - return identity matrix
		if (onLiborPeriodDiscretization(evaluationTime)) {
			RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
			for(int i=0;i<dLdL.length;i++) dLdL[i][i]=new RandomVariable(1.0);
		    return dLdL;
		}
		
		// Calculate dLdL. It is a (n-1)x n Matrix!
		RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors+1];
		double swapTenorLength = model.getLiborPeriodDiscretization().getTimeStep(0); // Model must have same tenor as swap!
		double timeOfFirstLiborPriorToEval = getPreviousLiborTime(evaluationTime);
		int timeIndexAtEvaluationTime = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int timeIndexAtFirstLiborPriorToEval = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(timeOfFirstLiborPriorToEval);
		
		for(int liborIndex = 0; liborIndex <numberOfRemainingLibors; liborIndex++){
			double liborTime = evaluationTime+liborIndex*swapTenorLength; // t+j*\Delta T
		    int    previousLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(liborTime);
		    double previousLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex);
		    double firstNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+1);
		    double secondNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+2);
		    double factor1 = (secondNextLiborTime-(liborTime+swapTenorLength))/(secondNextLiborTime-firstNextLiborTime);
		    double factor2 = (liborTime-previousLiborTime)/(firstNextLiborTime-previousLiborTime);
		    int    timeIndex = liborIndex==0 ? timeIndexAtFirstLiborPriorToEval : timeIndexAtEvaluationTime;
		    // Get Libors. We take the average here to avoid calculation of pseudo inverse of a matrix of RV. Should be fixed.
		    RandomVariableInterface previousLibor = model.getLIBOR(timeIndex, previousLiborIndex);     
		    RandomVariableInterface nextLibor     = model.getLIBOR(timeIndex, previousLiborIndex + 1); 
		    RandomVariableInterface logInterpol = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).log().mult(-factor1);
		                            logInterpol = logInterpol.add(previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).log().mult(-factor2)).exp();
		    // Set derivatives
		    dLdL[liborIndex][liborIndex]   = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).mult(logInterpol).mult(1-factor2);// dLdL_i-1
		    dLdL[liborIndex][liborIndex+1] = previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).mult(logInterpol).mult(1-factor1);
		}
		
		// dLdL is (n-1) x n matrix. Get PseudoInverse for all paths and then put it back together as RV
		return getPseudoInverse(dLdL);
	}
	
	
	/**Calculates the row vector dV/dL
	 * 
	 * @param evaluationTime The time at which the forward sensistivity dVdL is calculated
	 * @return The forward sensisivity dVdL (as a row vector)
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getValueLiborSensitivities(PortfolioInstrument product, 
			                                                    double evaluationTime) throws CalculationException{
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		RandomVariableDifferentiableInterface numeraire = (RandomVariableDifferentiableInterface) model.getNumeraire(evaluationTime);
		
		// Calculate forward sensitivities
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime);
		int numberOfSensis = evaluationTime == getNextLiborTime(evaluationTime) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
		RandomVariableInterface[] valueLiborSensitivities = new RandomVariableInterface[numberOfSensis];// exclude last libor
		int timeIndexAtEval = model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		
		// Set all entries of dVdL
		// Set dVdL for last libor which is already fixed (if applicable)
		int timeGridIndicator = 0;
		int lastLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		if(numberOfSensis!=numberOfRemainingLibors){
			timeGridIndicator = 1;
			double lastLiborTime = model.getLiborPeriodDiscretization().getTime(lastLiborIndex);
			RandomVariableInterface lastLibor = model.getLIBOR(model.getTimeDiscretization().getTimeIndex(lastLiborTime), lastLiborIndex);
			RandomVariableInterface dVdL = getProductValueDerivative(product,lastLibor);
			valueLiborSensitivities[0] = dVdL.mult(numeraire);
		}
		
		for(int liborIndex=lastLiborIndex+timeGridIndicator;liborIndex<model.getNumberOfLibors(); liborIndex++){
			RandomVariableInterface liborAtTimeIndex = model.getLIBOR(timeIndexAtEval, liborIndex);
		    RandomVariableInterface dVdL = getProductValueDerivative(product,liborAtTimeIndex);
		    valueLiborSensitivities[liborIndex-lastLiborIndex] = dVdL.mult(numeraire).getConditionalExpectation(conditionalExpectationOperator);
		}
		
		if(isUseTimeGridAdjustment){
		// Up to now dVdL is wrt the Libors on the LiborPeriodDiscretization. Adjust it such that we have dVdL wrt Libors starting at evaluationTime 
		RandomVariableInterface[][] dLdL = getLiborTimeGridAdjustment(evaluationTime, model);
		RandomVariableInterface[] dVdLAdjusted = multiply(valueLiborSensitivities,dLdL);
		return dVdLAdjusted; 
		} else return valueLiborSensitivities;
	}
	
	/**Calculates dV/dS where S are swap rates of the discount curve.
	 * 
	 * @param evaluationTime The time at which dVdS is calculated
	 * @return The row vector of sensitivities wrt swap rates from discount curve.
	 * @throws CalculationException 
	 */
	public RandomVariableInterface[] getDiscountCurveSensitivities(PortfolioInstrument product, double evaluationTime) throws CalculationException{
		// We calculate dV/dP * dP/dS. dV/dP is at t=0 since the curve starts at 0 !?!
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		final double shift = 0.0001;
		
		// Remove first entry from pillars if it is at time 0.
		int index = discountCurvePillars[0]==0 ? 1 : 0;
		double[] pillars = new double[discountCurvePillars.length-index];
		for(int i=0;i<pillars.length;i++) pillars[i]=discountCurvePillars[i+index];
		
		RandomVariableInterface  value = product.getProduct().getValue(evaluationTime, model); //.mult(model.getNumeraire(evaluationTime));
		Map<Long, RandomVariableInterface> gradientOfProduct = ((RandomVariableDifferentiableInterface) value).getGradient();
		int numberOfP = getNumberOfRemainingLibors(evaluationTime);
		
	
		int lastPillarIndex = evaluationTime>pillars[0] ? new TimeDiscretization(pillars).getTimeIndexNearestLessOrEqual(evaluationTime) : 0;
		
        double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		RandomVariableInterface[] dVdP = new RandomVariableInterface[pillars.length-lastPillarIndex];//numberOfP];
        DiscountCurve discountCurve = (DiscountCurve) model.getModel().getDiscountCurve();
        // Define new pillars. Another option is to calculate dV/dP wrt original Pillars, and do interpolation: dV/dP*dP/d\tilde{P}
//		TimeDiscretization curveTimes = new TimeDiscretization(evaluationTime, numberOfP, liborPeriodLength);
        TimeDiscretization curveTimes = new TimeDiscretization(pillars);
		RandomVariableInterface[] discountFactors = new RandomVariableInterface[pillars.length];//numberOfP+1];
		// get discount factors
		for(int i=0;i<pillars.length;i++) discountFactors[i]=discountCurve.getDiscountFactor(pillars[i]);//curveTimes.getTime(i));
		// dV(t)/dP(T_i;0)
		for(int i=lastPillarIndex;i<pillars.length;i++){
//			discountFactors[i]+=shift;
//			DiscountCurve newDiscountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve", pillars/*curveTimes.getAsDoubleArray()*/, discountFactors); 
//		    // get clone of LMM with shifted curve. 
//			Map<String,Object> dataModified = new HashMap<String,Object>();
//		    dataModified.put("discountCurve", newDiscountCurve);
//		    LIBORModelMonteCarloSimulationInterface newModel = (LIBORModelMonteCarloSimulation) model.getCloneWithModifiedData(dataModified);
//		    dVdP[i-lastPillarIndex] = currentProduct.getValue(evaluationTime, newModel).mult(newModel.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
//            dVdP[i-lastPillarIndex]=dVdP[i-lastPillarIndex].sub(value).div(shift).mult(discountCurve.getDiscountFactor(evaluationTime));
//		    discountFactors[i]-=shift;
			dVdP[i-lastPillarIndex] = getDerivative(gradientOfProduct, discountFactors[i]).getConditionalExpectation(conditionalExpectationOperator);
		}
		// Get dP(T_i;0)/dP(t+i\delta T;0): Linear interpolation on log value per time
		RandomVariableInterface[][] dPdP = new RandomVariableInterface[numberOfP][pillars.length-lastPillarIndex];
		for(int i=0;i<dPdP.length;i++){
			double discountTime = evaluationTime + (i+1) * liborPeriodLength;
			if(discountTime < pillars[0]) {
				double term = Math.pow(discountTime/pillars[0],2.0);
				dPdP[i][0]=discountFactors[0].invert().mult(term).mult(discountFactors[0].log().mult(term).exp());
				continue;
			}
			// Get upper and lower index
			int lowerIndex = curveTimes.getTimeIndexNearestLessOrEqual(discountTime); // as 0 is included in time discretization but nut in pillars
			lowerIndex = lowerIndex <0 ? 0 : lowerIndex;
			int upperIndex = lowerIndex+1;
			double delta = (discountTime-pillars[lowerIndex])/curveTimes.getTimeStep(lowerIndex);
			RandomVariableInterface summand1 = discountFactors[lowerIndex].log().mult((1-delta)/pillars[lowerIndex]);
			RandomVariableInterface summand2 = discountFactors[upperIndex].log().mult(delta/pillars[upperIndex]);
			RandomVariableInterface factor   = summand1.add(summand2).mult(discountTime).exp();
			//Math.exp(((1-delta)/pillars[lowerIndex]*Math.log(discountFactors[lowerIndex])+delta/pillars[upperIndex]*Math.log(discountFactors[upperIndex]))*discountTime);
			dPdP[i][lowerIndex-lastPillarIndex]=factor.div(discountFactors[lowerIndex]).mult((1-delta)/pillars[lowerIndex]*discountTime);
			dPdP[i][upperIndex-lastPillarIndex]=factor.div(discountFactors[upperIndex]).mult(delta/pillars[upperIndex]*discountTime);
		}
        dVdP = multiply(dVdP,getPseudoInverse(dPdP));
		RandomVariableInterface[][] dPdS = getBondSwapSensitivity(evaluationTime);
		return multiply(dVdP,dPdS);
	}
	
	
	/**Calculates dPdS in a single curve context. Used for calculating sensis with respect to discount curve.
	 * 
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return The sensitivity of the discount curve (bonds) wrt to swap rates of the same curve.
	 * @throws CalculationException 
	 */
	private RandomVariableInterface[][] getBondSwapSensitivity(double evaluationTime) throws CalculationException{
		int numberOfBonds = getNumberOfRemainingLibors(evaluationTime);
		RandomVariableInterface sum= new RandomVariable(0.0);
		RandomVariableInterface[][] dSdP = new RandomVariableInterface[numberOfBonds][numberOfBonds];
		for(int bondIndex=0;bondIndex<dSdP[0].length;bondIndex++){
			RandomVariableInterface bond = model.getNumeraire(evaluationTime+(bondIndex+1)*0.5).invert().mult(model.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
		    sum = sum.add(bond);
		    for(int swapIndex=0;swapIndex<dSdP.length;swapIndex++){
		    	if(swapIndex<bondIndex) dSdP[swapIndex][bondIndex] = new RandomVariable(0.0);
		    	else if(swapIndex==bondIndex) dSdP[swapIndex][bondIndex] = sum.add(1.0).sub(bond).mult(-1.0).div(sum.squared());
		    	else dSdP[swapIndex][bondIndex] = bond.sub(1.0).div(sum.squared());    	
		    }
		} 
		return getPseudoInverse(dSdP); // PseudoInverse == Inverse for n x n matrix.
	}
	
	private RandomVariableInterface[] getMeltedSensitivities(double evaluationTime, RandomVariableInterface[] sensitivities, String riskClass){
		int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};	
		// Get new riskFactor times
		int[] riskFactorDays = Arrays.stream(riskFactorsSIMM).filter(n -> n > (int)Math.round(365*evaluationTime)).map(n -> n-(int)Math.round(365*evaluationTime)).toArray();
	    // Find first bucket later than evaluationTime
		int firstIndex = IntStream.range(0, riskFactorsSIMM.length)
		                          .filter(i -> riskFactorsSIMM[i]>(int)Math.round(365*evaluationTime)).findFirst().getAsInt();
		//Calculate melted sensitivities
		RandomVariableInterface[] meltedSensis = new RandomVariableInterface[sensitivities.length-firstIndex];
		for(int i=0;i<meltedSensis.length;i++){
			meltedSensis[i]=sensitivities[i+firstIndex].mult(1.0-(double)Math.round(365*evaluationTime)/(double)riskFactorsSIMM[i+firstIndex]);
		}
		return getSensitivitiesOnBuckets(meltedSensis, riskClass, riskFactorDays);       
	}
	
	//----------------------------------------------------------------------------------------------------------------------------------
	// Some auxiliary functions
	//----------------------------------------------------------------------------------------------------------------------------------

	/**Calculate Pseudo Inverse of matrix of type RandomVariableInterface[][]
	 * 
	 * @param matrix The matrix for which the pseudo inverse is calculated
	 * @return The pseudo inverse of the matrix
	 */
    private RandomVariableInterface[][] getPseudoInverse(RandomVariableInterface[][] matrix){
    	double[][][] inv = new double[matrix[0].length][matrix.length][model.getNumberOfPaths()];
		double[][] matrixOnPath = new double[matrix.length][matrix[0].length];
		for(int pathIndex=0; pathIndex<model.getNumberOfPaths(); pathIndex++){
			// Get double[][] matrix on path
			for(int i=0;i<matrixOnPath.length;i++){
				for(int j=0;j<matrixOnPath[0].length;j++){
					matrixOnPath[i][j]=matrix[i][j]==null ? 0 : matrix[i][j].get(pathIndex);
				}
			}
		    // Get Pseudo Inverse 
		    RealMatrix pseudoInverse = new SingularValueDecomposition(MatrixUtils.createRealMatrix(matrixOnPath)).getSolver().getInverse();
		    for(int j=0;j<pseudoInverse.getColumnDimension();j++){
			    double[] columnValues = pseudoInverse.getColumn(j);
			    for(int i=0;i<pseudoInverse.getRowDimension();i++){
				    inv[i][j][pathIndex]= columnValues[i];
			    }
		    }
		}
		// Wrap to RandomVariableInterface[][]
		RandomVariableInterface[][] pseudoInverse = new RandomVariableInterface[matrix[0].length][matrix.length];
		for(int i=0;i<pseudoInverse.length; i++){
			for(int j=0;j<pseudoInverse[0].length; j++){
				pseudoInverse[i][j] = new RandomVariable(0.0 /*should be evaluationTime*/,inv[i][j]);
			}
		}
		return pseudoInverse;
    }
   

	public static RandomVariableInterface[][] multiply(RandomVariableInterface[][] A,RandomVariableInterface[][] B){
		RandomVariableInterface[][] AB = new RandomVariableInterface[A.length][B.length];
		RandomVariableInterface ABproduct;
		for(int i=0;i<A.length;i++){
			for(int j=0; j<B.length; j++){
				AB[i][j] = new RandomVariable(0.0);
				for(int k=0;k<B.length;k++) {
					if(A[i][k]==null || B[k][j]==null) {ABproduct = new RandomVariable(0.0);}
					else {ABproduct = A[i][k].mult(B[k][j]);}
					AB[i][j]=AB[i][j].add(ABproduct);
				}
			}
		}
		return AB;
	}
	
	public static RandomVariableInterface[] multiply(RandomVariableInterface[] A,RandomVariableInterface[][] B){
		RandomVariableInterface[] AB = new RandomVariableInterface[B[0].length];
		RandomVariableInterface ABproduct;
		for(int i=0;i<B[0].length;i++){
				AB[i] = new RandomVariable(0.0);
				for(int k=0;k<A.length;k++) {
					if(A[k]==null || B[k][i]==null) {ABproduct = new RandomVariable(0.0);}
					else {ABproduct = A[k].mult(B[k][i]);}
					AB[i]=AB[i].add(ABproduct);
				}
		}
		return AB;
	}
	
	/**Calculates the derivative of the current portfolio product with respect to the specified parameter, dV/dX 
	 * 
	 * @param parameter The parameter with respect to which the derivative is calculated
	 * @return dV/dX 
	 * @throws CalculationException
	 */
	private RandomVariableInterface getProductValueDerivative(PortfolioInstrument product, 
			                                                  RandomVariableInterface parameter) throws CalculationException{

		RandomVariableInterface derivative = product.getGradient().get(((RandomVariableDifferentiableInterface)parameter).getID());
		return derivative==null ? new RandomVariable(0.0) : derivative;
	}
	
	private RandomVariableInterface getDerivative(Map<Long, RandomVariableInterface> gradient, RandomVariableInterface parameter){
		RandomVariableInterface derivative = gradient.get(((RandomVariableDifferentiableInterface)parameter).getID());
		return derivative==null ? new RandomVariable(0.0) : derivative;
	}
	
	private int getNumberOfRemainingLibors(double evaluationTime){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getNumberOfLibors()-nextLiborIndex;
	}
	
	private double getNextLiborTime(double evaluationTime){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex);
	}
	
	private double getPreviousLiborTime(double evaluationTime){
		if(evaluationTime==0) return 0.0;
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex-1);
	}
	
		
	private static ArrayList<RandomVariableInterface> getRegressionBasisFunctions(RandomVariableInterface[] libors, int order) {
		ArrayList<RandomVariableInterface> basisFunctions = new ArrayList<RandomVariableInterface>();
		// Create basis functions - here: 1, S, S^2, S^3, S^4
		for(int liborIndex=0; liborIndex<libors.length;liborIndex++){
		  for(int powerOfRegressionMonomial=0; powerOfRegressionMonomial<=order; powerOfRegressionMonomial++) {
			  basisFunctions.add(libors[liborIndex].pow(powerOfRegressionMonomial));
		  }
		}
		return basisFunctions;
	}
	
	/*
	 *  Getters and Setters
	 */
	
	public AbstractLIBORMonteCarloProduct getPortfolioProduct(int index) {
		return portfolioProducts[index].getProduct();
	}

	
	public WeightToLiborAdjustmentMethod getLiborWeightMethod(){
		return this.liborWeightMethod;
	}
	
	
	private boolean onLiborPeriodDiscretization(double evaluationTime){
		return (evaluationTime == getNextLiborTime(evaluationTime));
	}
	
	private void setConditionalExpectationOperator(double evaluationTime) throws CalculationException{
		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
        RandomVariableInterface[] regressor = new RandomVariableInterface[2];
        regressor[0]= model.getLIBOR(evaluationTime, evaluationTime,evaluationTime+model.getLiborPeriodDiscretization().getTimeStep(0));
		regressor[1]= model.getLIBOR(evaluationTime, evaluationTime, model.getLiborPeriodDiscretization().getTime(model.getNumberOfLibors()-1));
       	ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(regressor, 2);
//		Alternative definition of regressors
//		RandomVariableInterface[] libors = getRemainingLibors(evaluationTime, model);
//		ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(libors, 1 /*polyNomialOrder*/);
//		ConditionalExpectationEstimatorInterface conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));
       	this.conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));

	}
	
	// Delete later..
	public void setUseTimeGridAdjustment(boolean method){ 
		this.isUseTimeGridAdjustment = method;
	}
	
	
	private RandomVariableInterface[][] getPseudoInverse(double[][] matrix){
		RealMatrix pseudoInverse = new SingularValueDecomposition(MatrixUtils.createRealMatrix(matrix)).getSolver().getInverse();
		RandomVariableInterface[][] inv = new RandomVariableInterface[matrix[0].length][matrix.length];
		for(int j=0;j<pseudoInverse.getColumnDimension();j++){
		    double[] columnValues = pseudoInverse.getColumn(j);
		    for(int i=0;i<pseudoInverse.getRowDimension();i++){
			    inv[i][j]= new RandomVariable(columnValues[i]);
		    }		    
	    }
		return inv;
	}
	
	private void clearProducts(){
		for(int i=0;i<portfolioProducts.length;i++){
			portfolioProducts[i].clear();
		}
	}
	
	public PortfolioInstrument[] getProducts(){
		return portfolioProducts;
	}
	
	private PortfolioInstrument[] createPortfolioInstruments(SIMMClassifiedProduct[] classifiedProducts){
		PortfolioInstrument[] products = new PortfolioInstrument[classifiedProducts.length];
		for(int i=0;i<products.length;i++) products[i]= new PortfolioInstrument(classifiedProducts[i]);
		return products;
	}
	
	
} 