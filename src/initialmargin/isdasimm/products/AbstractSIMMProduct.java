package initialmargin.isdasimm.products;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.ArrayUtils;

import initialmargin.isdasimm.aggregationscheme.CalculationSchemeInitialMarginISDA;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.SensitivityMode;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.WeightMode;
import initialmargin.isdasimm.sensitivity.SIMMSensitivityCalculation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

/** This class contains the functions and methods which are shared by all products to be considered for 
 *  initial margin (MVA) calculation by the SIMM. 
 *  
 * @author Mario Viehmann
 *
 */
public abstract class AbstractSIMMProduct implements SIMMProductInterface {
	
    // Product classification within ISDA SIMM
    private String   productClass;      // RatesFX, Credit, 
    private String[] riskClass;         // InterestRate, CreditQ, CreditNonQ, Equity, Commodity
    private String[] curveIndexNames;   // e.g. OIS & Libor6m
    private String   currency;
    private boolean  hasOptionality;    // determines the relevance of vega and curvature risk (e.g. Swap has no curvature risk)
    private String   bucketKey;         // can be null (e.g. in risk class InterestRate it is null because the bucket is given by the currency
    
    // Further variables
    protected Map<Long, RandomVariableInterface> gradient = null;// Same for all evaluationTimes; Is reset for different products
    protected boolean isGradientOfDeliveryProduct = false;
    protected RandomVariableInterface exerciseIndicator;
    /*
     * Model and Product are separate! When we call "getInitialMargin(time, model)", we set modelCache = model! 
     * Thus, we can check if the model has changed. If it has changed, we have to re-calculate the gradient and clear the sensitivity maps.
     */
    protected LIBORModelMonteCarloSimulationInterface modelCache; 
    protected double lastEvaluationTime = -1;
    protected ConditionalExpectationEstimatorInterface conditionalExpectationOperator;
    protected AbstractSIMMSensitivityCalculation sensitivityCalculationScheme;
    private   CalculationSchemeInitialMarginISDA simmScheme;
    
    
    public static final String[]  IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};
    public static final RandomVariableInterface[] zeroBucketsIR = IntStream.range(0, IRMaturityBuckets.length).mapToObj(i->new RandomVariable(0.0)).toArray(RandomVariableInterface[]::new);
    
    // Define the sensitivity maps.
    /**
     *  The map of delta sensitivities at a specific time. This map is filled once per evaluation time step and the
     *  function <code> getSensitivity <code> defined in class <code> AbstractSIMMProduct </code> which is called in
     *  <code> MarginSchemeIRDelta </code> picks the sensitivies for a specified riskClass, curveIndexName and maturityBucket
     *  from this map. This map may - in contrast to the second map "exactDeltaCache" - contain interpolated sensitivities.
     */
    private HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,
            HashMap<String/*maturityBucket*/,RandomVariableInterface>>>> deltaAtTime = new HashMap<String,List<HashMap<String,HashMap<String,RandomVariableInterface>>>>(); // currently only for InterestRate riskClass
    
    /**
     * The cache for the exact delta sensitivities as given by AAD (or analytic). Unlike the map
     * "deltaAtTime", this map is not cleared if evaluationTime differs from lastEvaluationTime
     */
    private HashMap<Double /*time*/,List<HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,
            RandomVariableInterface[]>>>>> exactDeltaCache = new HashMap<Double /*time*/,List<HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>>>();
    
	//private RandomVariableInterface vegaSensitivity=null; 
    
    /**
     * The cache of numeraire adjustments used in the evaluation of this product in the <code> LIBORMarketModel <code>.
     * This data is the basis of the OIS curve sensitivities, which we calculate by applying AAD to the numeraire adjustments
     */
    protected Map<Double,RandomVariableInterface> numeraireAdjustmentMap = new HashMap<>();
    public static boolean isPrintSensis = false;
    
    /**Wraps an <code> AbstractLIBORMonteCarloProduct </code> into a product classified according to the SIMM methodology requirement.
     * 
     * @param productClass The SIMM product class of this product (RatesFx etc.)
     * @param riskClass The SIMM risk class of this product (InterestRate etc.)
     * @param curveIndexNames The name of the relevant curves for this product (OIS, Libor6m etc.)
     * @param currency The currency of this product
     * @param bucketKey The SIMM bucket key of this product (null for risk class InterestRate)
     * @param hasOptionality True if this product is not linear
     */
     public AbstractSIMMProduct(String   productClass,
		                        String[] riskClass,     // One product may contribute to several risk Classes
		                        String[] curveIndexNames,
		                        String   currency,
		                        String   bucketKey,
		                        boolean  hasOptionality){
    	
	   this.productClass = productClass; 
	   this.riskClass = riskClass;
	   this.curveIndexNames = curveIndexNames;
	   this.currency=currency;
	   this.hasOptionality = hasOptionality;
	   this.bucketKey = bucketKey;
	   
    }
     
    
    @Override
    public RandomVariableInterface getInitialMargin(double evaluationTime, LIBORModelMonteCarloSimulationInterface model, String calculationCCY) throws CalculationException{
    	return getInitialMargin(evaluationTime, model, calculationCCY, SensitivityMode.Exact, WeightMode.Constant, 0, true, false, true);
    }
 	
 	public RandomVariableInterface getInitialMargin(double evaluationTime, 
 			                                        LIBORModelMonteCarloSimulationInterface model, 
 			                                        String calculationCCY,
 			                                        SensitivityMode sensitivityMode,
 			                                        WeightMode liborWeightMode,
 			                                        double interpolationStep, 			                                        
 			                                        boolean isUseTimeGridAdjustment, 
 			                                        boolean isUseAnalyticSwapSensis, 
 			                                        boolean isConsiderOISSensitivities) throws CalculationException{
 		
 		if(evaluationTime >= getFinalMaturity()) return new RandomVariable(0.0);
 		
 		if(this.modelCache==null || !model.equals(this.modelCache) || (sensitivityCalculationScheme!=null && (sensitivityMode !=sensitivityCalculationScheme.getSensitivityMode() || liborWeightMode !=sensitivityCalculationScheme.getWeightMode()))) { // At inception (t=0) or if the model is reset            
 			setGradient(model); // Set the (new) gradient. The method setModel also clears the sensitivity maps and sets the model as modelCache.
 	        this.exerciseIndicator = null;
 	        this.exactDeltaCache.clear();
 	        this.sensitivityCalculationScheme = new SIMMSensitivityCalculation(sensitivityMode, liborWeightMode, interpolationStep, model, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensitivities);
 			this.simmScheme= new CalculationSchemeInitialMarginISDA(this,calculationCCY);
 		}
 		
 		return simmScheme.getValue(evaluationTime); 
 	}
     
 	
    @Override
    public RandomVariableInterface getSensitivity(String productClass, 
             									  String riskClass, 
             									  String maturityBucket, // only for IR and Credit risk class, null otherwise
             									  String curveIndexName, // null if riskClass is not IR
             									  String bucketKey,      // currency for IR otherwise bucket number
             									  String riskType, double evaluationTime) throws SolverException, CloneNotSupportedException, CalculationException{

         RandomVariableInterface result = null;	RandomVariableInterface[] maturityBucketSensis; // Sensitivities mapped on the SIMM Buckets

         if(!hasOptionality && riskType!="delta") return new RandomVariable(0.0);		   

         if(evaluationTime!=lastEvaluationTime) clearMaps();  // Clear the deltaSensitivity Map. It needs to be reset at each time step.
           
         if(productClass==this.productClass && Arrays.asList(this.riskClass).contains(riskClass)){

            switch(riskType){
               case("delta"): 
                  switch(riskClass){
                      case("InterestRate"):

                      if(Arrays.asList(curveIndexNames).contains(curveIndexName) && bucketKey==this.currency){
                         // There exists a sensitivity. Check if the sensitivities (on all maturityBuckets) have already been calculated for given riskClass and riskType)

                         if(!deltaAtTime.containsKey(riskClass) || !deltaAtTime.get(riskClass).stream().filter(n-> n.containsKey(curveIndexName)).findAny().isPresent()){

                            // The sensitivities need to be calculated for the given riskClass and riskType                     	            		                     
                            maturityBucketSensis = sensitivityCalculationScheme.getDeltaSensitivities(this, riskClass, curveIndexName, evaluationTime, modelCache);
                                             
                            if(isPrintSensis && curveIndexName=="Libor6m") {                            	
                            	System.out.println(evaluationTime + "\t" + maturityBucketSensis[3].getAverage() + "\t" + maturityBucketSensis[4].getAverage() + "\t"+ maturityBucketSensis[5].getAverage() + "\t"+ maturityBucketSensis[6].getAverage() + "\t"+ maturityBucketSensis[7].getAverage() + "\t"+ maturityBucketSensis[8].getAverage() + "\t"+maturityBucketSensis[9].getAverage() + "\t"+maturityBucketSensis[10].getAverage() + "\t"+maturityBucketSensis[11].getAverage());
                            }
                            // Create a new element of the curveIndex List for given risk class		         
                            HashMap<String,HashMap<String,RandomVariableInterface>> curveIndexNameexactDeltaCache = new HashMap<String,HashMap<String,RandomVariableInterface>>();
                            HashMap<String,RandomVariableInterface> bucketSensitivities = new HashMap<String,RandomVariableInterface>();

                            for(int i=0;i<IRMaturityBuckets.length;i++) bucketSensitivities.put(IRMaturityBuckets[i], maturityBucketSensis[i]);
                            curveIndexNameexactDeltaCache.put(curveIndexName,bucketSensitivities);
                            
                            // Check if list already exist
                            if(deltaAtTime.containsKey(riskClass)) deltaAtTime.get(riskClass).add(curveIndexNameexactDeltaCache);
                            else {
                              List<HashMap<String/*curveIndexName*/,HashMap<String/*maturityBucket*/,RandomVariableInterface>>> list = new ArrayList<HashMap<String/*curveIndexName*/,HashMap<String/*maturityBucket*/,RandomVariableInterface>>>();
                              list.add(curveIndexNameexactDeltaCache);
                              deltaAtTime.put(riskClass, list);
                            }					    					        					          
                         }
                         result = deltaAtTime.get(riskClass).stream().filter(n -> n.containsKey(curveIndexName)).findFirst().get().get(curveIndexName).get(maturityBucket);			                  
                    } else result = new RandomVariable(0.0); // There exists no delta Sensi for risk Class InterestRate
                    break;
                    // @Todo Add sensitivity calculation for the subsequent cases
                    case("CreditQ"): 
                    case("CreditNonQ"):
                    case("FX"):
                    case("Commodity"):
                    case("Equity"): result = null;
                  } break;
             // @Todo Add sensitivity calculation for the subsequent cases
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

      this.lastEvaluationTime = evaluationTime;
      return result;
    } // end getSensitivity()
    
 	
    /** Returns the cache of numeraire adjustments used in the last valuation of this product. We need the numeraire adjustments
     *  to calculate the sensitivities w.r.t. the OIS curve.
     * 
     * @return The cache of numeraire adjustments from the Libor market model
     * @throws CalculationException
     */
 	public Map<Double, RandomVariableInterface> getNumeraireAdjustmentMap() throws CalculationException{
	   if(this.numeraireAdjustmentMap==null) {
		  modelCache.clearNumeraireAdjustmentCache();
		  getLIBORMonteCarloProduct().getValue(0.0,modelCache);
		  this.numeraireAdjustmentMap = modelCache.getNumeraireAdjustmentMap();
	   }
	   return this.numeraireAdjustmentMap;
	}
	  
 	@Override
 	public Map<Long, RandomVariableInterface> getGradient() throws CalculationException {
		  
 	   if(gradient==null) {
	      // Clear cache of numeraire adjustments of the model to capture the numeraire adjustments from the product valuation
		  modelCache.clearNumeraireAdjustmentCache();		
		  // Calculate the product value as of time 0.
		  RandomVariableDifferentiableInterface productValue = (RandomVariableDifferentiableInterface) getLIBORMonteCarloProduct().getValue(0.0, modelCache);		     
		  // Get the map of numeraire adjustments used specifically for this product
		  this.numeraireAdjustmentMap.putAll(modelCache.getNumeraireAdjustmentMap());			 
		  // Calculate the gradient
		  Map<Long, RandomVariableInterface> gradientOfProduct = productValue.getGradient();
		  this.gradient = gradientOfProduct;
	   }		  
	   return this.gradient;		  
	}
 	  
 	 
 	/** Set the cache of exact delta sensitivities (calculated by AAD or analytically for Swaps). This cache is used in the class 
 	 *  <code> SIMMSensitivityCalculation <code> to obtain the sensitivities used for melting and interpolation.
 	 * 
 	 * @param riskClass The risk class
 	 * @param curveIndexName The name of the curve (OIS or Libor6m)
 	 * @param time The time for which the forward sensitivity is calculated
 	 * @param model The Libor Market model
 	 * @throws SolverException
 	 * @throws CloneNotSupportedException
 	 * @throws CalculationException
 	 */
    private void setExactDeltaCache(String riskClass, String curveIndexName, 
 				                   double time, LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{
 			
       // Calculate the sensitivities 
 	   RandomVariableInterface[] deltaSensis = sensitivityCalculationScheme.getExactDeltaSensitivities(this, curveIndexName, riskClass, time, model);
 			
 	   // Create a new element of the curveIndex List for given risk class		         
 	   HashMap<String,RandomVariableInterface[]> curveIndexNameDeltaCache = new HashMap<String,RandomVariableInterface[]>();
 	   curveIndexNameDeltaCache.put(curveIndexName,deltaSensis);
 	        
 	   // Check if the list of riskClasses in the HashMap already exist
 	   if(exactDeltaCache.containsKey(new Double(time))){
 	      if(exactDeltaCache.get(new Double(time)).stream().filter(n-> n.containsKey(riskClass)).findAny().isPresent()){
 	         exactDeltaCache.get(new Double(time)).stream().filter(n-> n.containsKey(riskClass)).findAny().get().get(riskClass).add(curveIndexNameDeltaCache);
 	        	  
 	      	           
 	      } else { // there is no of risk classes. Set it.
 	   	     List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>> curveList = new ArrayList<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>();
 	 	     curveList.add(curveIndexNameDeltaCache);
 	 	     HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>> riskClassMap = new HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>();
 	 	     riskClassMap.put(riskClass, curveList);	 	          
 	 	     exactDeltaCache.get(time).add(riskClassMap);	        	   
 	      }
 	           	            
 	   } else { // no entry at time
 	      List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>> curveList = new ArrayList<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>();
 	      curveList.add(curveIndexNameDeltaCache);
 	      HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>> riskClassMap = new HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>();
 	      riskClassMap.put(riskClass, curveList);
 	      List<HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>> riskClassList = new ArrayList<HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,RandomVariableInterface[]>>>>();
 	      riskClassList.add(riskClassMap);
 	      exactDeltaCache.put(time, riskClassList);
 	   }
 				
 	}
    
    /** Calculate the forward derivatives of the product w.r.t. the Libors at a given evaluation time.
     *  These derivatives are not w.r.t. the libors on the Libor period discretization, but w.r.t. the general libors 
     *  L(t+i\Delta_T,t+(i+1)\Delta_T;t). This function is called by the subclasses in the overridden functions
     *  <code> getValueLiborSensitivities <code>.
     * 
     * @param evaluationTime The time for which the forward sensitivities are calculated
     * @param model The libor market model
     * @return The forward Libor sensitivities of this product 
     * @throws CalculationException
     */
    protected RandomVariableInterface[] getValueLiborSensitivitiesAAD(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		setConditionalExpectationOperator(evaluationTime);
		RandomVariableDifferentiableInterface numeraire = (RandomVariableDifferentiableInterface) model.getNumeraire(evaluationTime);
		
		// Calculate forward sensitivities
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime,model);
		int numberOfSensis = evaluationTime == getNextLiborTime(evaluationTime,model) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
		RandomVariableInterface[] valueLiborSensitivities = new RandomVariableInterface[numberOfSensis];
		int timeIndexAtEval = model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		
		// Set all entries of dVdL
		// Set dVdL for last libor which is already fixed (if applicable)
		int timeGridIndicator = 0;
		int lastLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		
		if(numberOfSensis!=numberOfRemainingLibors){
			timeGridIndicator = 1;
			double lastLiborTime = model.getLiborPeriodDiscretization().getTime(lastLiborIndex);
			RandomVariableInterface lastLibor = model.getLIBOR(model.getTimeDiscretization().getTimeIndex(lastLiborTime), lastLiborIndex);
			RandomVariableInterface dVdL = getDerivative(lastLibor);
			valueLiborSensitivities[0] = dVdL.mult(numeraire);
		}
		
		for(int liborIndex=lastLiborIndex+timeGridIndicator;liborIndex<model.getNumberOfLibors(); liborIndex++){
			RandomVariableInterface liborAtTimeIndex = model.getLIBOR(timeIndexAtEval, liborIndex);
		    RandomVariableInterface dVdL = getDerivative(liborAtTimeIndex);
		    valueLiborSensitivities[liborIndex-lastLiborIndex] = dVdL.mult(numeraire).getConditionalExpectation(conditionalExpectationOperator);
		}

		if(sensitivityCalculationScheme.isUseTimeGridAdjustment){
		    // Up to now dVdL is wrt the Libors on the LiborPeriodDiscretization. Adjust it such that we have dVdL wrt Libors starting at evaluationTime 
		    RandomVariableInterface[][] dLdL = AbstractSIMMSensitivityCalculation.getLiborTimeGridAdjustment(evaluationTime, model);
		    RandomVariableInterface[] dVdLAdjusted = AbstractSIMMSensitivityCalculation.multiply(valueLiborSensitivities,dLdL);
		
		    return dVdLAdjusted; 
		} else return valueLiborSensitivities;
		
	}
 	  
 	/** Calculate the sensitivities w.r.t. the OIS curve: dV/dS. These sensitivities are calculated using the numeraire adjustments at future 
 	 *  discount times (times at which the model has called the <code> getNumeraire <code> function of the <code> LIBORMarketModel <code>.
 	 *  The function calculates the AAD derivatives dV/dA i.e. w.r.t. the adjustments at the future discount times and converts them into 
 	 *  dV/dS (sensitivties w.r.t. the swap rates). 
 	 *  dV/dP may be provided as input to this function such that only the conversion to dV/dS is perfromed. This is useful if we have 
 	 *  already obtained the sensitivities dV/dP analytically. 
 	 *  This function is called inside the subclasses in the overridden function 
 	 *  <code> getDiscountCurveSensitivities(String riskClass, double evaluationTime) <code>.
 	 * 
 	 * @param evaluationTime The time as of which the forward sensitivities are considered
 	 * @param futureDiscountTimes (may be null) The times at which the the numeraire has been used in the last valuation
 	 * @param dVdP (may be null) The sensitivities w.r.t. the OIS bond. Only used if dV/dP is given analytically.
 	 * @param riskClass The risk class
 	 * @param model The Libor market model
 	 * @return The forward sensitivities w.r.t. the OIS curve 
 	 * @throws CalculationException
 	 */
    protected RandomVariableInterface[] getDiscountCurveSensitivities(double evaluationTime, double[] futureDiscountTimes, 
    		                                                             RandomVariableInterface[] dVdP, String riskClass, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
       
       if(dVdP == null || futureDiscountTimes == null){ //i.e. need to calculate it with AAD
       
    	   // Get map with all numeraire adjustments used for this product
	       Map<Double, RandomVariableInterface> adjustmentMap = getNumeraireAdjustmentMap();
	
	       // Return zero if evaluationTime is later than the last time where an adjustment is available (i.e. the last time where a cash flow occurred)
	       if(!adjustmentMap.keySet().stream().filter(time -> time > evaluationTime).findAny().isPresent()){
		  
		      RandomVariableInterface zero = new RandomVariable(0.0);
		      return AbstractSIMMSensitivityCalculation.mapSensitivitiesOnBuckets(new RandomVariableInterface[]{zero}, riskClass, new int[]{17},model);
	   
	       }
	
	       // Calculate adjustment at evaluationTime
	       RandomVariableInterface adjustmentAtEval = model.getNumeraireAdjustment(evaluationTime);
			
	       // Get all adjustments after evaluationTime
	       Set<Double> relevant = new HashSet<Double>();
	       relevant = adjustmentMap.keySet().stream().filter(entry -> entry>evaluationTime).collect(Collectors.toCollection(HashSet::new));
	       adjustmentMap.keySet().retainAll(relevant);
	       futureDiscountTimes = ArrayUtils.toPrimitive(Arrays.stream(adjustmentMap.keySet().toArray()).sorted().toArray(Double[]::new));
	
	       //Calculate derivative w.r.t. adjustment
	       int numberOfSwaps = adjustmentMap.size();
	       dVdP = new RandomVariableInterface[numberOfSwaps];
	
	       setConditionalExpectationOperator(evaluationTime);
	       RandomVariableInterface numeraireAtEval  = model.getNumeraire(evaluationTime);
	
	       for(int i=0;i<dVdP.length;i++){
		
		      // Calculate dVdA
		      RandomVariableInterface adjustment = adjustmentMap.get(futureDiscountTimes[i]);
		      RandomVariableInterface dVdA = getDerivative(adjustment).getConditionalExpectation(conditionalExpectationOperator).mult(numeraireAtEval);
		
		      // Calculate dV(t)/dP(t_cf;t) where t_cf are the cash flow times of this product
		      RandomVariableInterface bond = model.getForwardBondLibor(futureDiscountTimes[i],evaluationTime);
		      //RandomVariableInterface test = model.getNumeraire(evaluationTime).div(model.getNumeraire(futureDiscountTimes[i])).mult(adjustment).div(adjustmentAtEval).getConditionalExpectation(conditionalExpectationOperator);
		      dVdP[i] = dVdA.mult(adjustment.squared()).mult(-1.0).div(bond).div(adjustmentAtEval);
		
	       }
       }
	
	   // Perform a log-linear interpolation of the discount factors to obtain dP(t_cf;t)/dP(t+i\Delta_T;t).
	   int numberOfP = getNumberOfRemainingLibors(evaluationTime, model);
	   RandomVariableInterface[][] dPdP = new RandomVariableInterface[futureDiscountTimes.length][numberOfP];
	
	   double deltaT = model.getLiborPeriodDiscretization().getTimeStep(0);
	   TimeDiscretizationInterface timesP = new TimeDiscretization(evaluationTime, numberOfP, deltaT);
	   for(int cfIndex=0; cfIndex<dPdP.length; cfIndex++){
		   int lowerIndex = timesP.getTimeIndexNearestLessOrEqual(futureDiscountTimes[cfIndex]);
		   double alpha = (futureDiscountTimes[cfIndex]-timesP.getTime(lowerIndex))/deltaT;
		   Arrays.fill(dPdP[cfIndex], new RandomVariable(0.0));
		   RandomVariableInterface bondLowerIndex = lowerIndex==0 ? new RandomVariable(1.0) : model.getForwardBondOIS(timesP.getTime(lowerIndex),evaluationTime);
		   RandomVariableInterface bondUpperIndex = model.getForwardBondOIS(timesP.getTime(lowerIndex+1),evaluationTime);
		   RandomVariableInterface bondAtCF       = model.getForwardBondOIS(futureDiscountTimes[cfIndex],evaluationTime);		   
		   dPdP[cfIndex][lowerIndex] = bondAtCF.mult(1-alpha).div(bondLowerIndex);
		   dPdP[cfIndex][lowerIndex+1] = bondAtCF.mult(alpha).div(bondUpperIndex);
	   }
	
	   // Calulate dV(t)/dP(t+i\Delta_T;t)
	   dVdP = AbstractSIMMSensitivityCalculation.multiply(dVdP,dPdP);
	
	   // Calculate dP(t+i\Delta_T;t)/dS_i(t) and dV(t)/dS_i(t)
	   RandomVariableInterface[][] dPdS = AbstractSIMMSensitivityCalculation.getBondSwapSensitivity(evaluationTime,model);
	   RandomVariableInterface[]   dVdS = AbstractSIMMSensitivityCalculation.multiply(dVdP,dPdS);
	
	   return AbstractSIMMSensitivityCalculation.mapSensitivitiesOnBuckets(dVdS, riskClass, null, model);
   }
 	  
 	/** Clear the time dependent delta cache and the vega sensitivity. 
 	 *  This is performed always upon change of the evaluation time of initial margin.
 	 */
	public void clearMaps(){
		   if(this.deltaAtTime!=null) this.deltaAtTime.clear();
		   //this.vegaSensitivity = null;
	}
	
	@Override
    public RandomVariableInterface[] getExactDeltaFromCache(double time, String riskClass, String curveIndexName) throws SolverException, CloneNotSupportedException, CalculationException{
    
    	if(!exactDeltaCache.containsKey(time) || !exactDeltaCache.get(time).stream().filter(n->n.containsKey(riskClass)).findAny().isPresent()){
        	
    		for(String curveName : curveIndexNames) setExactDeltaCache(riskClass, curveName, time, modelCache);
        
    	}
    	
    	return exactDeltaCache.get(time).stream().filter(n->n.containsKey(riskClass)).findAny().get().get(riskClass).stream().filter(n->n.containsKey(curveIndexName)).findAny().get().get(curveIndexName);
    				
    }
	
 	
    /*
     * Getters and setters
     */
    
    /** Set the libor market model in the model cache, set the gradient of the product w.r.t. the new model.  
     * 
     * @param model The libor market model
     * @throws CalculationException
     */
    public void setGradient(LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
    	// If the model is set, we must clear all maps and set the gradient to null.
    	clearMaps(); //...the maps containing sensitivities are reset to null (they have to be recalculated under the new model)
    	this.modelCache = model;
    	this.gradient = null;
    	this.gradient = getGradient();
    	this.isGradientOfDeliveryProduct = false; // for (bermudan) swaptions
    	
    }
    
    public abstract AbstractLIBORMonteCarloProduct getLIBORMonteCarloProduct();
    
    public String getProductClass(){
    	return this.productClass;
    }
    public String[] getRiskClasses(){
    	return this.riskClass;
    }
    public String[] getCurveIndexNames(){
    	return this.curveIndexNames;
    }
    public String getCurrency(){
    	return this.currency;
    }
    public boolean getHasOptionality(){
    	return this.hasOptionality;
    }
    public String getBucketKey(){
    	return this.bucketKey;
    }
    
    public HashMap<Double /*time*/,List<HashMap<String/*RiskClass*/,List<HashMap<String/*curveIndexName*/,
    RandomVariableInterface[]>>>>> getExactDeltaCache(){
    	return this.exactDeltaCache;
    }
    
    public void clearDeltaCache(){
    	this.exactDeltaCache.clear();
    }
    
    public void setSIMMSensitivityCalculation(AbstractSIMMSensitivityCalculation sensitivityCalculation){
    	this.sensitivityCalculationScheme = sensitivityCalculation;
    }
    
    public void setNullExerciseIndicator(){
    	this.exerciseIndicator=null;
    }
   
    
	protected RandomVariableInterface getDerivative(RandomVariableInterface parameter) throws CalculationException{
		Map<Long, RandomVariableInterface> gradient = getGradient();
		RandomVariableInterface derivative = gradient.get(((RandomVariableDifferentiableInterface)parameter).getID());
		return derivative==null ? new RandomVariable(0.0) : derivative;
	}
	
	protected double getNextLiborTime(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex);
	}
	
	protected int getNumberOfRemainingLibors(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getNumberOfLibors()-nextLiborIndex;
	}
	
	protected double getPreviousLiborTime(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		if(evaluationTime==0) return 0.0;
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex-1);
	}
	
    
}
