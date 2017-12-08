package initialmargin.isdasimm.products;

import java.util.ArrayList;
import java.util.Arrays;

import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.SimpleSwap;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

/** This class describes a Swap for SIMM initial margin (MVA) calculation.
 * 
 * @author Mario Viehmann
 *
 */
public class SIMMSimpleSwap extends AbstractSIMMProduct{

	// SIMM classification
	static final String productClass = "RatesFX";
	static final String[] riskClass = new String[]{"InterestRate"};
	private SimpleSwap swap = null;
	
	
	/** Construct a swap as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param swap The swap as AbstractLIBORMonteCarloProduct
	 * @param curveIndexNames The names of the curves for this swap
	 * @param currency The currency of this swap
	 */
	public SIMMSimpleSwap(SimpleSwap swap, String[] curveIndexNames, String currency) {
        super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, false /*hasOptionality*/);
        this.swap = swap;
    }
	
	
	/** Construct a swap as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param fixingDates
	 * @param paymentDates
	 * @param swapRates
	 * @param isPayFix
	 * @param notional
	 * @param curveIndexNames
	 * @param currency
	 */
	public SIMMSimpleSwap(double[] fixingDates, 
			              double[] paymentDates, 
			              double[] swapRates, 
			              boolean isPayFix, 
			              double  notional, 
			              String[] curveIndexNames, String currency) {
		super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, false /*hasOptionality*/);
		this.swap = new SimpleSwap(fixingDates,paymentDates,swapRates,isPayFix,notional);
		
	}
	
	
	/** Calculate Swap Sensitivities dV/dL (Libor) or dV/dP (OIS) analytically for a Swap 
	 * 
	 * @param evaluationTime The time of evaluation
	 * @param fixingDates The fixing times of the swap floating leg
	 * @param swapRates The swap rates. May be <code> null <code> (only relevant for derivative w.r.t. discount curve).
	 * @param periodLength The constant period length of this swap
	 * @param model The LIBOR model used for simulation of the Libors
	 * @param withRespectTo "Libor" or "OIS"
	 * @return
	 * @throws CalculationException
	 */
	public static RandomVariableInterface[] getAnalyticSensitivities(double evaluationTime, 
            													     double[] fixingDates,
            													     double[] swapRates, /* Only relevant for OIS derivative*/
            													     double   periodLength,
            													     double   notional,
            													     LIBORModelMonteCarloSimulationInterface model,
            													     String withRespectTo) throws CalculationException{

		    RandomVariableInterface[] sensis = null;

			// periodIndex: Index of the swap period at evaluationTime
			int periodIndex = new TimeDiscretization(fixingDates).getTimeIndexNearestLessOrEqual(evaluationTime); 
			periodIndex = periodIndex < 0 ? 0 : periodIndex;
			
			// firstLiborIndex: Index of the Libor on the first period of the swap
			int currentLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
			int firstLiborIndex   = fixingDates[0] > evaluationTime ? model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(fixingDates[0]):currentLiborIndex;
            
			switch(withRespectTo){
			
			   case("Libor"):
				   
				   int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
				   int numberOfRemainingLibors = model.getNumberOfLibors()- nextLiborIndex;
				
				   int numberOfSensis = evaluationTime == model.getLiborPeriodDiscretization().getTime(nextLiborIndex) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
				   sensis = new RandomVariableInterface[numberOfSensis]; 
				   Arrays.fill(sensis, new RandomVariable(0.0));
				   
				   // return zero if evaluationTime > last payment date
				   if(evaluationTime>=fixingDates[fixingDates.length-1]+periodLength) return sensis; 
				   
			       // Actual Sensitivity Calculation: dV/dL = P(T,t)*periodLength
			       for(int liborIndex=currentLiborIndex;liborIndex<numberOfSensis+currentLiborIndex;liborIndex++){
				       int i = liborIndex < firstLiborIndex ? 0 : liborIndex-firstLiborIndex+1;
				       if(!(i>fixingDates.length-periodIndex || i==0) ){ 		
					       double paymentTime = fixingDates[periodIndex+i-1]+periodLength;
					       sensis[liborIndex-currentLiborIndex]=model.getForwardBondOIS(paymentTime, evaluationTime).mult(periodLength);			
				       }
			       }
			       break;
			       
			   case("OIS"):
				   
				   // Actual Sensitivity Calculation: dV/dL = P(T,t)*periodLength
				   numberOfSensis = fixingDates.length-periodIndex;
			   
			       // return zero if evaluationTime > last payment date, i.e. if numberOfSensis <= 0
			       if(numberOfSensis <= 0) return new RandomVariableInterface[]{new RandomVariable(0.0)};
			       
			       sensis = new RandomVariableInterface[numberOfSensis];
			       int timeIndex = model.getTimeIndex(evaluationTime);
			       timeIndex = timeIndex < 0 ? -timeIndex-2 : timeIndex;
			       firstLiborIndex   = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(fixingDates[0]);
			       
			       for(int liborIndex=0;liborIndex<numberOfSensis;liborIndex++){
					       double swapRate = swapRates[periodIndex+liborIndex];
					       int timeIndexFixing = model.getTimeIndex(fixingDates[periodIndex+liborIndex]);
					       sensis[liborIndex]=model.getLIBOR(Math.min(timeIndex,timeIndexFixing), firstLiborIndex+periodIndex+liborIndex).sub(swapRate).mult(periodLength);				     			       
			       }	       
			       break;	  
			}
			sensis = Arrays.stream(sensis).map(n->n.mult(notional)).toArray(RandomVariableInterface[]::new);
			return sensis;
	}
	
	public RandomVariableInterface[] getAnalyticSensitivities(double evaluationTime,
															  double periodLength,
		     												  LIBORModelMonteCarloSimulationInterface model,
		     												  String withRespectTo) throws CalculationException{

		return getAnalyticSensitivities(evaluationTime, swap.getFixingDates(), swap.getSwapRates(), periodLength, swap.getNotional(), model, withRespectTo);
			    
	}
	
	@Override
	public AbstractLIBORMonteCarloProduct getLIBORMonteCarloProduct() {
		return this.swap;
	}

	
	@Override
	public RandomVariableInterface[] getValueLiborSensitivities(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		
		if(sensitivityCalculationScheme.isUseAnalyticSwapSensitivities) {
			
			RandomVariableInterface[] swapSensis = getAnalyticSensitivities(evaluationTime, model.getLiborPeriodDiscretization().getTimeStep(0), model, "Libor");
			RandomVariableInterface[][] dLdL =  AbstractSIMMSensitivityCalculation.getLiborTimeGridAdjustment(evaluationTime, model);
			return AbstractSIMMSensitivityCalculation.multiply(swapSensis,dLdL);
			
		}
		
		return getValueLiborSensitivitiesAAD(evaluationTime, model);
		
	}
		
		
	@Override
	public RandomVariableInterface[] getDiscountCurveSensitivities(String riskClass, 
																   double evaluationTime,
																   LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		
		double[] futureDiscountTimes = null; // the times of the times after evaluation time at which the numeraire has been used for this product
		RandomVariableInterface[] dVdP = null;
		
		if(sensitivityCalculationScheme.isUseAnalyticSwapSensitivities) {
			
			// Return zero if evaluationTime is later than the last time where an adjustment is available (i.e. the last time where a cash flow occurred)
			if(!Arrays.stream(swap.getPaymentDates()).filter(time -> time > evaluationTime).findAny().isPresent()){
				
				RandomVariableInterface zero = new RandomVariable(0.0);
				return AbstractSIMMSensitivityCalculation.mapSensitivitiesOnBuckets(new RandomVariableInterface[]{zero}, riskClass, new int[]{17},model);
			
			}
			
			dVdP = getAnalyticSensitivities(evaluationTime, model.getLiborPeriodDiscretization().getTimeStep(0), model, "OIS");

		    futureDiscountTimes = Arrays.stream(swap.getPaymentDates()).filter(n -> n > evaluationTime).toArray();
			
		} 
			
		return getDiscountCurveSensitivities(evaluationTime, futureDiscountTimes, dVdP /* null => use AAD*/, riskClass, model);
		   
    }

	@Override
	public RandomVariableInterface getExerciseIndicator(double time) {
		return new RandomVariable(1.0);
	}
	
	@Override
	public double getFinalMaturity(){
		return swap.getPaymentDates()[swap.getPaymentDates().length-1];
	}
	
	@Override
	public double getMeltingResetTime(){
		return 0;
	}
	
	@Override
	public void setConditionalExpectationOperator(double evaluationTime) throws CalculationException{
		
		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
        RandomVariableInterface[] regressor = new RandomVariableInterface[2];
        regressor[0]= modelCache.getLIBOR(evaluationTime, evaluationTime,evaluationTime+modelCache.getLiborPeriodDiscretization().getTimeStep(0));
		regressor[1]= modelCache.getLIBOR(evaluationTime, evaluationTime, modelCache.getLiborPeriodDiscretization().getTime(modelCache.getNumberOfLibors()-1));
       	ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(regressor, 2);
       	this.conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));

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


}
