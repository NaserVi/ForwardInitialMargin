package initialmargin.isdasimm.products;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.SimpleSwap;
import initialmargin.isdasimm.changedfinmath.products.Swaption;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.stochastic.RandomVariableInterface;

/** This class describes a Swaption for SIMM initial margin (MVA) calculation.
 *  Both delivery types "Physical" and "CashSettled" are supported.
 * 
 * @author Mario Viehmann
 *
 */
public class SIMMSwaption extends AbstractSIMMProduct{
	
	// SIMM classification
	static final String productClass = "RatesFX";
	static final String[] riskClass = new String[]{"InterestRate"};
	
	// Swap after exercise
	private SimpleSwap swap = null; 
	
	private Swaption swaption;
	public enum DeliveryType {Physical, CashSettled};
	private DeliveryType deliveryType;
			
		
	/** Construct a swaption as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param swaption
	 * @param deliveryType
	 * @param curveIndexNames
	 * @param currency
	 * @throws CalculationException
	 */
	public SIMMSwaption(Swaption swaption, DeliveryType deliveryType, String[] curveIndexNames, String currency) throws CalculationException {
	    super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, true /*hasOptionality*/);
	    this.swaption = swaption;
	    this.deliveryType = deliveryType;
	    if(deliveryType==DeliveryType.Physical){
			
	    	this.swap = new SimpleSwap(swaption.getFixingDates(), swaption.getPaymentDates(),swaption.getSwaprates(),true /*isPayFix*/, 
		                               swaption.getNotional());
			
	     }
	}
		
	/** Construct a swaption as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param exerciseDate
	 * @param fixingDates
	 * @param paymentDates
	 * @param swapRates
	 * @param notional
	 * @param deliveryType
	 * @param curveIndexNames
	 * @param currency
	 * @throws CalculationException
	 */
    public SIMMSwaption(double exerciseDate, double[] fixingDates, double[] paymentDates, double[] swapRates, double notional, 
				        DeliveryType deliveryType, String[] curveIndexNames, String currency) throws CalculationException {
			
    	super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, true /*hasOptionality*/);
			
		this.swaption = new Swaption(exerciseDate,fixingDates,paymentDates,swapRates,notional);
		this.deliveryType = deliveryType;
		if(deliveryType==DeliveryType.Physical){
				
			this.swap = new SimpleSwap(fixingDates, paymentDates, swapRates, true /*isPayFix*/, notional);
				
		}
			
	}
		
		
    @Override
    public AbstractLIBORMonteCarloProduct getLIBORMonteCarloProduct() {
		return this.swaption;
	}

	
	@Override
	public RandomVariableInterface[] getValueLiborSensitivities(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		
		if(deliveryType == DeliveryType.Physical && evaluationTime >= swaption.getExerciseDate()){
			
			if(sensitivityCalculationScheme.isUseAnalyticSwapSensitivities){
			   
			   // Calculate sensis analytically
			   RandomVariableInterface[] swapSensis = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime, swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "Libor");
			   RandomVariableInterface indicator = getExerciseIndicator(evaluationTime);
			   swapSensis = Arrays.stream(swapSensis).map(n->n.mult(indicator)).toArray(RandomVariableInterface[]::new);
			   // Get time grid adjustment
			   RandomVariableInterface[][] dLdL =  AbstractSIMMSensitivityCalculation.getLiborTimeGridAdjustment(evaluationTime, model);
			   
			   return AbstractSIMMSensitivityCalculation.multiply(swapSensis,dLdL);
			
			} else setSwapGradient();
					
		}
		
		return getValueLiborSensitivitiesAAD(evaluationTime, model);
	}
	

	@Override
	public RandomVariableInterface[] getDiscountCurveSensitivities(String riskClass, 
																   double evaluationTime,
																   LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		
		double[] futureDiscountTimes = null; // the times of the times after evaluation time at which the numeraire has been used for this product
		RandomVariableInterface[] dVdP = null;
		
		if(deliveryType == DeliveryType.Physical && evaluationTime >= swaption.getExerciseDate()){
			
			// Return zero if evaluationTime is later than the last time where an adjustment is available (i.e. the last time where a cash flow occurred)
			if(!Arrays.stream(swap.getPaymentDates()).filter(time -> time > evaluationTime).findAny().isPresent()){
				return zeroBucketsIR;		
			}
			
			if(sensitivityCalculationScheme.isUseAnalyticSwapSensitivities) {
						
			   dVdP = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime,swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "OIS");
			   RandomVariableInterface indicator = getExerciseIndicator(evaluationTime);
			   dVdP = Arrays.stream(dVdP).map(n->n.mult(indicator)).toArray(RandomVariableInterface[]::new);
			   futureDiscountTimes = Arrays.stream(swap.getPaymentDates()).filter(n -> n > evaluationTime).toArray();
			
			} else setSwapGradient();
			
		}

	    return getDiscountCurveSensitivities(evaluationTime, futureDiscountTimes, dVdP /* null => use AAD*/, riskClass, model);
	}
		  

	@Override
	public RandomVariableInterface getExerciseIndicator(double time) throws CalculationException {
		if(exerciseIndicator==null) exerciseIndicator = swaption.getExerciseIndicator(modelCache);
		return this.exerciseIndicator;
	}
	
	@Override
	public double getFinalMaturity(){
		return deliveryType==DeliveryType.Physical ? swap.getPaymentDates()[swap.getPaymentDates().length-1] : swaption.getExerciseDate();
	}
	
	@Override
	public double getMeltingResetTime(){
		return swaption.getExerciseDate();
	}
	
	@Override
	public void setConditionalExpectationOperator(double evaluationTime) throws CalculationException{
		
		// Swaption: Set paths on which we have not exercised to zero
		RandomVariableInterface indicator = new RandomVariable(1.0);
		if(evaluationTime>=swaption.getExerciseDate()) indicator = getExerciseIndicator(evaluationTime); // 1 if exercised on this path		   
		
		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
        RandomVariableInterface[] regressor = new RandomVariableInterface[2];
        regressor[0]= modelCache.getLIBOR(evaluationTime, evaluationTime,evaluationTime+modelCache.getLiborPeriodDiscretization().getTimeStep(0));
		regressor[1]= modelCache.getLIBOR(evaluationTime, evaluationTime, modelCache.getLiborPeriodDiscretization().getTime(modelCache.getNumberOfLibors()-1));
       	ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(regressor, 2, indicator);
       	this.conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));

	}
	
	private static ArrayList<RandomVariableInterface> getRegressionBasisFunctions(RandomVariableInterface[] libors, int order, RandomVariableInterface indicator) {
		ArrayList<RandomVariableInterface> basisFunctions = new ArrayList<RandomVariableInterface>();
		// Create basis functions - here: 1, S, S^2, S^3, S^4
		
		for(int liborIndex=0; liborIndex<libors.length;liborIndex++){
		  for(int powerOfRegressionMonomial=0; powerOfRegressionMonomial<=order; powerOfRegressionMonomial++) {
			  basisFunctions.add(libors[liborIndex].pow(powerOfRegressionMonomial).mult(indicator));
		  }
		  
		}
		return basisFunctions;
	}
	
	
	/** Set the gradient of the swap in case of physical exercise. 
	 * 
	 * @throws CalculationException
	 */
    private void setSwapGradient() throws CalculationException{
		 if(!super.isGradientOfDeliveryProduct){
		    // Clear cache of numeraire adjustments of the model to capture the numeraire adjustments from the product valuation
		    modelCache.clearNumeraireAdjustmentCache();
		    // Calculate the product value as of time 0.
		    RandomVariableInterface indicator = getExerciseIndicator(swaption.getExerciseDate()+0.0001);
		    RandomVariableDifferentiableInterface productValue = (RandomVariableDifferentiableInterface) swap.getValue(0.0, modelCache).mult(indicator);
		    // Get the map of numeraire adjustments used specifically for this product
	        super.numeraireAdjustmentMap.putAll(modelCache.getNumeraireAdjustmentMap());
		    // Calculate the gradient
		    Map<Long, RandomVariableInterface> gradientOfProduct = productValue.getGradient(); 
		    // Set the gradient
		    super.gradient = gradientOfProduct;
		    super.isGradientOfDeliveryProduct = true;
		 }	 
	}

    
	public DeliveryType getDeliveryType(){
			return this.deliveryType;
	}

}

