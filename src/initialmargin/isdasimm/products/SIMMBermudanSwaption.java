package initialmargin.isdasimm.products;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.BermudanSwaption;
import initialmargin.isdasimm.changedfinmath.products.SimpleSwap;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;

/** This class describes a Bermudan Swaption for SIMM initial margin (MVA) calculation.
 *  Both exercise types "Callable" and "Cancelable" are supported.
 * 
 * @author Mario Viehmann
 *
 */
public class SIMMBermudanSwaption extends AbstractSIMMProduct{

	// SIMM classification
	static final String productClass = "RatesFX";
	static final String[] riskClass = new String[]{"InterestRate"};
	
	// Swap after exercise
	private SimpleSwap swap = null; 
	//private boolean isUseAnalyticSensitivities = true;
	
	private BermudanSwaption bermudan;
	
	public enum ExerciseType {Callable, Cancelable};
	private ExerciseType exerciseType;
	private Map<String, RandomVariableInterface[]> swapSensitivityMap = new HashMap<>();

			
	/** Construct a bermudan swaption as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param bermudan
	 * @param curveIndexNames
	 * @param currency
	 * @throws CalculationException
	 */
	public SIMMBermudanSwaption(BermudanSwaption bermudan, String[] curveIndexNames, String currency) throws CalculationException {
	    
		super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, true /*hasOptionality*/);
	    
	    this.bermudan = bermudan;	    
	    this.exerciseType = bermudan.getIsCallable() ? ExerciseType.Callable : ExerciseType.Cancelable;
	    if(exerciseType==ExerciseType.Callable) this.swap = bermudan.getSwap();
	}
		
	
	/** Construct a bermudan swaption as a product for the SIMM. Initial margin and MVA can be calculated for this product.
	 * 
	 * @param fixingDates
	 * @param periodLengths
	 * @param paymentDates
	 * @param periodNotionals
	 * @param swapRates
	 * @param isPeriodStartDateExerciseDate
	 * @param exerciseType
	 * @param curveIndexNames
	 * @param currency
	 * @throws CalculationException
	 */
    public SIMMBermudanSwaption(double[] fixingDates, double[] periodLengths, double[] paymentDates, double[] periodNotionals,
    		                    double[] swapRates, boolean[] isPeriodStartDateExerciseDate, ExerciseType exerciseType, String[] curveIndexNames, String currency) throws CalculationException {
			super(productClass, riskClass, curveIndexNames, currency, null /*bucketKey*/, false /*hasOptionality*/);
			boolean isCallable = exerciseType == ExerciseType.Callable ? true : false;
			this.bermudan = new BermudanSwaption(isPeriodStartDateExerciseDate,fixingDates,periodLengths, paymentDates,periodNotionals,swapRates, isCallable);
		    this.exerciseType = exerciseType;
		    if(exerciseType==ExerciseType.Callable) this.swap = bermudan.getSwap();
		    
	}
		
		
    @Override
    public AbstractLIBORMonteCarloProduct getLIBORMonteCarloProduct() {
		return this.bermudan;
	}

	
	@Override
	public RandomVariableInterface[] getValueLiborSensitivities(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		
		// Get Bermudan sensitivities
		RandomVariableInterface[] bermudanSensis = getValueLiborSensitivitiesAAD(evaluationTime, model);
		
		// Set exercise indicator
		RandomVariableInterface[] indicatorOrig = new RandomVariableInterface[]{getExerciseIndicator(evaluationTime)};
		RandomVariableInterface[] indicator = indicatorOrig.clone();
		indicator[0] = indicator[0].sub(0.5); // i.e. -0.5: not exercised, 0.5: exercised
		
		// Set sensitivities on exercised paths
		if(evaluationTime >= bermudan.getLastValuationExerciseTime().getMin()){
			
			switch(exerciseType){
			   
			   case Callable:
				   
			      // Calculate sensis analytically
				  RandomVariableInterface[] swapSensisAna = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime, swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "Libor");				  
				  // Get time grid adjustment
				  RandomVariableInterface[][] dLdL =  AbstractSIMMSensitivityCalculation.getLiborTimeGridAdjustment(evaluationTime, model);	   
				  RandomVariableInterface[] swapSensis = AbstractSIMMSensitivityCalculation.multiply(swapSensisAna,dLdL);
			   
				  //for(int i=0; i<swapSensis.length; i++) System.out.println("SwapsensisANA " + swapSensis[i].getAverage());
				  if(evaluationTime>=bermudan.getExerciseTimes()[bermudan.getExerciseTimes().length-1]) {
			    		
					    // Set sensis of not exercised paths to zero
			    	    for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0],swapSensis[i],new RandomVariable(0.0));
			 			    	    			    	    
				  } else {
			    	    			    	
				        // Set sensitivities on paths: Bermudan sensis if not exercised, swap sensis if exercised.
				        for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0], swapSensis[i], bermudanSensis[i]);
			   	
				  }
				  
				  break;
				  
			   case Cancelable:
			
				   // Set the sensitivities on exercised paths to zero
				   for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0], new RandomVariable(0.0), bermudanSensis[i]);

				   break;
				
			   default:
				   break;
			}		
		}
		
		//for(int i=0; i<bermudanSensis.length; i++) System.out.println("L " + evaluationTime + " "+ bermudanSensis[i].getAverage());
		
		return bermudanSensis;
	}
	

	@Override
	public RandomVariableInterface[] getDiscountCurveSensitivities(String riskClass, 
																   double evaluationTime,
																   LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		
		// Get Bermudan sensitivities
		double[] futureDiscountTimes = null; // the times of the times after evaluation time at which the numeraire has been used for this product
		RandomVariableInterface[] dVdP = null;
		RandomVariableInterface[] bermudanSensis = getDiscountCurveSensitivities(evaluationTime, futureDiscountTimes, dVdP /* null => use AAD*/, riskClass, model);

		// Set exercise indicator
		RandomVariableInterface[] indicatorOrig = new RandomVariableInterface[]{getExerciseIndicator(evaluationTime)};
		RandomVariableInterface[] indicator = indicatorOrig.clone();
		indicator[0] = indicator[0].sub(0.5); // i.e. -0.5: not exercised, 0.5: exercised
		
		// Set sensitivities on exercised paths 
		if(evaluationTime >= bermudan.getLastValuationExerciseTime().getMin()){
			
			switch(exerciseType){
			
			   case Callable:
				   
				    // Return zero if evaluationTime is later than the last time where an adjustment is available (i.e. the last time where a cash flow occurred)
					if(!Arrays.stream(swap.getPaymentDates()).filter(time -> time > evaluationTime).findAny().isPresent()){
												
								   RandomVariableInterface zero = new RandomVariable(0.0);
								   return AbstractSIMMSensitivityCalculation.mapSensitivitiesOnBuckets(new RandomVariableInterface[]{zero}, riskClass, new int[]{17},model);
											
					}
			
				   // Get Swap Sensitivities analytically		
				   dVdP = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime,swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "OIS");
				   futureDiscountTimes = Arrays.stream(swap.getPaymentDates()).filter(n -> n > evaluationTime).toArray();
				   RandomVariableInterface[] swapSensis = getDiscountCurveSensitivities(evaluationTime, futureDiscountTimes, dVdP, riskClass, model);		    
				   					   
				   if(evaluationTime>=bermudan.getExerciseTimes()[bermudan.getExerciseTimes().length-1]) {
			    		
					    // Set sensis of not exercised paths to zero
			    	    for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0],swapSensis[i],new RandomVariable(0.0));
			 			    	    
			    	    
				   } else {
			    	    			    	
				        // Set sensitivities on paths: Bermudan sensis if not exercised, swap sensis if exercised.
				        for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0], swapSensis[i], bermudanSensis[i]);
			   	
				   }
				   
				   break;
				   
			   case Cancelable:
				   
				   // Set the sensitivities on exercised paths to zero
				   for(int i=0; i<bermudanSensis.length;i++) bermudanSensis[i] = bermudanSensis[i].barrier(indicator[0], new RandomVariable(0.0), bermudanSensis[i]);
				   
				   break;
				   
			   default:
				   break;
			}
			
		}

		//for(int i=0; i<bermudanSensis.length; i++) System.out.println("OIS " + evaluationTime + " "+ bermudanSensis[i].getAverage());
		
	    return bermudanSensis;
	    
	}
	
	
	/** This function sets the sensitivities on the exercised paths to either swap sensitivities (for a Callable) 
	 *  or to zero (for a Cancelable). This is necessary as the melting function in class <code> SIMMSensitivityCalculation <code>
	 *  does not know on which paths we have already exercised. 
	 * 
	 * @param evaluationTime The time of evaluation
	 * @param curveIndexName The name of the curve
	 * @param meltedBermudanSensis The bermudan sensitivities obtained by melting at evaluation time.
	 * @return The true melted bermudan sensitivities considering the paths on which we have exercised.
	 * @throws CalculationException
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 */
	public RandomVariableInterface[] changeMeltedSensisOnExercisedPaths(double evaluationTime, String curveIndexName, RandomVariableInterface[] meltedBermudanSensis) throws CalculationException, SolverException, CloneNotSupportedException{
				
		if(evaluationTime >= bermudan.getLastValuationExerciseTime().getMin()){
		   
		   RandomVariableInterface indicator = getExerciseIndicator(evaluationTime).sub(0.5);
		   
		   switch(exerciseType){
		      
		      case Callable:
		         // Get swap sensis
		         RandomVariableInterface[] swapSensis = getSwapSensitivitiesFromCache(curveIndexName);
		   
		         //for(int i=0; i<swapSensis.length; i++) System.out.println("Swapsensis " + swapSensis[i].getAverage());
		         // Melt swap sensis
		         double initialMeltingTime = swap.getStartTime();
		         RandomVariableInterface[] meltedSwapSensis = sensitivityCalculationScheme.getMeltedSensitivities(this, swapSensis, initialMeltingTime, evaluationTime, curveIndexName, "InterestRate");
		   
		         // Put swap sensis on exercised paths
		         //for(int i=0; i<meltedSwapSensis.length; i++) System.out.println("SwapsensisMelted " + meltedSwapSensis[i].getAverage());
		         //for(int i=0; i<meltedBermudanSensis.length; i++) System.out.println("Bersensis1 " + meltedBermudanSensis[i].getAverage());
		         
		         for(int i=0;i <meltedBermudanSensis.length;i++) meltedBermudanSensis[i] = meltedBermudanSensis[i].barrier(indicator, meltedSwapSensis[i], meltedBermudanSensis[i]);
		      
		         //for(int i=0; i<meltedBermudanSensis.length; i++) System.out.println("Bersensis2 " + meltedBermudanSensis[i].getAverage());
				 if(evaluationTime>=bermudan.getExerciseTimes()[bermudan.getExerciseTimes().length-1]) {
			    		
					  // Set sensis of not exercised paths to zero
			    	  for(int i=0; i<meltedBermudanSensis.length;i++) meltedBermudanSensis[i] = meltedBermudanSensis[i].barrier(indicator,meltedSwapSensis[i],new RandomVariable(0.0));
			 			    	    	
			    	  //for(int i=0; i<meltedBermudanSensis.length; i++) System.out.println("Bersensis3 " + meltedBermudanSensis[i].getAverage());
				 } else {
			    	    			    	
				      // Set sensitivities on paths: Bermudan sensis if not exercised, swap sensis if exercised.
				      for(int i=0; i<meltedBermudanSensis.length;i++) meltedBermudanSensis[i] = meltedBermudanSensis[i].barrier(indicator, meltedSwapSensis[i], meltedBermudanSensis[i]);
			   	
				 }
				 
				 //for(int i=0; i<meltedBermudanSensis.length; i++) System.out.println("Swapsensis " + meltedBermudanSensis[i].getAverage());
				 
		         break;
		      
		      case Cancelable:
		    	 // Set sensis on exercised paths to zero
		    	 for(int i=0;i <meltedBermudanSensis.length;i++) meltedBermudanSensis[i] = meltedBermudanSensis[i].barrier(indicator, new RandomVariable(0.0), meltedBermudanSensis[i]);
		   
		    	 break;
		    	 
		      default:
		    	  break;
		   }
		   
		}
		
		return meltedBermudanSensis;
			
	}
	
	/** Get the sensitivities of the underlying swap. This function is called when melting sensitivities at each time step.
	 * 
	 * @param curveIndexName The name of the curve
	 * @return The sensitivities of the underlying swap from the cache "swapSensitivityMap".
	 * @throws CalculationException
	 */
	private RandomVariableInterface[] getSwapSensitivitiesFromCache(String curveIndexName) throws CalculationException{
		
		if(!swapSensitivityMap.containsKey(curveIndexName)){
			double evaluationTime = swap.getStartTime();
			// Get OIS Sensis	
			RandomVariableInterface[] dVdP = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime,swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "OIS");
			double[] futureDiscountTimes = Arrays.stream(swap.getPaymentDates()).filter(n -> n > evaluationTime).toArray();
			RandomVariableInterface[] swapSensisOIS = getDiscountCurveSensitivities(evaluationTime, futureDiscountTimes, dVdP, "InterestRate", model);		    
            swapSensitivityMap.put("OIS", swapSensisOIS);
            
            // Get forward curve sensis
			RandomVariableInterface[] swapSensisAna = SIMMSimpleSwap.getAnalyticSensitivities(evaluationTime, swap.getFixingDates(), swap.getSwapRates(), model.getLiborPeriodDiscretization().getTimeStep(0), swap.getNotional(), model, "Libor");				  
			// Get time grid adjustment
			RandomVariableInterface[][] dLdL =  AbstractSIMMSensitivityCalculation.getLiborTimeGridAdjustment(evaluationTime, model);	   
			RandomVariableInterface[] dVdL = AbstractSIMMSensitivityCalculation.multiply(swapSensisAna,dLdL);
            
			// Calculate dV/dS = dV/dL * dL/dS
		    RandomVariableInterface[] swapSensisLibor = sensitivityCalculationScheme.getValueSwapSensitivities(evaluationTime, dVdL, model);  
		        //for(int i=0;i<swapSensisLibor.length;i++) System.out.println("dVdL " + swapSensisLibor[i].getAverage());
		        // Map Sensitivities on SIMM Buckets
		    swapSensisLibor = AbstractSIMMSensitivityCalculation.mapSensitivitiesOnBuckets(swapSensisLibor, "InterestRate" /*riskClass*/, null, model);
			
			
			swapSensitivityMap.put("Libor6m", swapSensisLibor);
			
			
		} 
		
		return swapSensitivityMap.get(curveIndexName);
	}
			  

	@Override
	public RandomVariableInterface getExerciseIndicator(double time) throws CalculationException {
		RandomVariableInterface indicator = bermudan.getLastValuationExerciseTime();
		indicator = indicator.barrier(new RandomVariable(indicator.sub(time+0.00001)), new RandomVariable(0.0), new RandomVariable(1.0));
	    return indicator;
	}
	
	@Override
	public double getFinalMaturity(){
		return bermudan.getPaymentDates()[bermudan.getPaymentDates().length-1];
	}
	
	@Override
	public double getMeltingResetTime(){
		return bermudan.getLastValuationExerciseTime().getMin();
	}
	
	@Override
	public void setConditionalExpectationOperator(double evaluationTime) throws CalculationException{
		
		// Bermudan Swaption: Set paths on which we have already exercised to zero.
	    RandomVariableInterface indicator = getExerciseIndicator(evaluationTime).barrier(new RandomVariable(getExerciseIndicator(evaluationTime).sub(0.5)), new RandomVariable(0.0), new RandomVariable(1.0));	    
		
		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
        RandomVariableInterface[] regressor = new RandomVariableInterface[2];
        regressor[0]= model.getLIBOR(evaluationTime, evaluationTime,evaluationTime+model.getLiborPeriodDiscretization().getTimeStep(0)).mult(indicator);
		regressor[1]= model.getLIBOR(evaluationTime, evaluationTime, model.getLiborPeriodDiscretization().getTime(model.getNumberOfLibors()-1)).mult(indicator);
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
	
	
	public ExerciseType getExerciseType(){
			return this.exerciseType;
	}
	
	
}
