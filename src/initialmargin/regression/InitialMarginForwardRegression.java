package initialmargin.regression;

import java.util.ArrayList;

//import com.marioviehmann.interestrate.products.Portfolio;
import initialmargin.regression.changedfinmath.products.*;

import net.finmath.exception.CalculationException;
import net.finmath.functions.NormalDistribution;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
//import com.marioviehmann.interestrate.*;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * This class implements the Dynamic Initial Margin by Regression as described in 
 * https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2911167
 * 
 * @author Mario Viehmann
 *
 */
public class InitialMarginForwardRegression {
	  private final double confidenceLevel = 0.99;
	  private final double MPOR = 0.03;
      private enum Method {SIMPLE, LSQREGRESSION, DELTAGAMMA};
      private Method method;
      private int polynomialOrder;
      
      private LIBORModelMonteCarloSimulationInterface model;
      public Portfolio portfolio;
      
      
      public InitialMarginForwardRegression(Portfolio portfolio,
    		               LIBORModelMonteCarloSimulationInterface model,
    		               int polynomialOrder,
    		               String method){
    	  this.model=model;
    	  this.portfolio = portfolio;
    	  this.polynomialOrder=polynomialOrder;
    	  this.method = Method.valueOf(method.toUpperCase());
      }
      
      
      public double getInitialMargin(double evaluationTime) throws CalculationException{
     	 
     	 RandomVariableInterface variance;
     	 double initialMargin = 0;
     	 switch (method){	   
 		   case DELTAGAMMA: // @ToDo delta and gamma as coefficients of the regression function
 			   break;
 		   
 		   case LSQREGRESSION:
 			   variance = getVarianceForecast(evaluationTime, model);
 			   double normalQuantile = NormalDistribution.inverseCumulativeDistribution(confidenceLevel);
 			   RandomVariableInterface initialMarginPathwise = variance.sqrt().mult(normalQuantile);
 			   initialMargin = initialMarginPathwise.getAverage();
 		       break;

 		   case SIMPLE:
 			   //RandomVariableInterface change = simple.getCleanPortfolioValueChange(evaluationTime);
 			   //for(int i=0;i<simple.getCleanPortfolioValueChange(evaluationTime).size(); i++) System.out.println(change.get(i));
 			   initialMargin = -getCleanPortfolioValueChange(evaluationTime).getQuantile(confidenceLevel);
 		       break;
 		       
 		   default:
 			   break;
     	 }
         return initialMargin;
      }
      
      public double[] getInitialMargin(TimeDiscretization initialMarginTimes) throws CalculationException{
     	 double[] initialMargin = new double[initialMarginTimes.getNumberOfTimes()];
     	 for(int timeIndex = 0; timeIndex<initialMarginTimes.getNumberOfTimes(); timeIndex++){
     		 initialMargin[timeIndex]=getInitialMargin(initialMarginTimes.getTime(timeIndex));
     	 }
     	 return initialMargin;
      }
     
      
      /**Calculates the clean portfolio value change, i.e. V(t+MPOR)-V(t)+CF({t,t+MPOR})
       * 
       * @param time the time t at which the marginal period of risk (MPOR) starts
       * @return 
       * @throws CalculationException
       */
      public RandomVariableInterface getCleanPortfolioValueChange(double time) throws CalculationException{
    	  
    	  double lastFixingTime = model.getLiborPeriodDiscretization().getTime(model.getLiborPeriodDiscretization().getTimeIndex(portfolio.getInitialLifeTime())-1);
    	  RandomVariableInterface cashFlows = portfolio.getCF(time, time+MPOR, model);
    	  
    	  RandomVariableInterface initialValue = portfolio.getValue(time, model);
    	  initialValue = initialValue.sub(cashFlows);
    	  
    	  if(time>0 && time < lastFixingTime && !(portfolio.getProducts()[0] instanceof Swap)) { // For swap we go forward along paths
    		  ConditionalExpectationEstimatorInterface condExpOperatorInitial = getConditionalExpectationEstimatorLibor(time, model);
    		  initialValue = initialValue.getConditionalExpectation(condExpOperatorInitial);
    	  }
    
    	  RandomVariableInterface finalValue = portfolio.getValue(time+MPOR, model);
    	  if(time+MPOR<lastFixingTime && !(portfolio.getProducts()[0] instanceof Swap)) {
    		  ConditionalExpectationEstimatorInterface condExpOperatorFinal = getConditionalExpectationEstimatorLibor(time+MPOR, model);
    		  finalValue = finalValue.getConditionalExpectation(condExpOperatorFinal);
    	  }
    	  
    	  //RandomVariableInterface finalValue = portfolio.getValue(time+MPOR, model);
    	  return finalValue.sub(initialValue); //.add(cashFlows);
      }
      
      /** Calculates the forecast of the variance of the clean portfolio value change over the marginal period of risk for all time points on a time discretization.
       * 
       * @param forwardVaRTimes The times for which the variance of the clean portfolio value change over the marginal period of risk is calculated
       * @param model Interface implementing the Libor Market Model 
       * @return
       * @throws CalculationException
       */
      public RandomVariableInterface[] getVarianceForecast(TimeDiscretizationInterface forwardVaRTimes, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
    	RandomVariableInterface[] VaRForecast = new RandomVariableInterface[forwardVaRTimes.getNumberOfTimes()];
  		for(int timeIndex = 0; timeIndex < forwardVaRTimes.getNumberOfTimes(); timeIndex++){
  			VaRForecast[timeIndex] = getVarianceForecast(forwardVaRTimes.getTime(timeIndex), model);
  		}
  		return VaRForecast;
  	}
  	
      

  	public RandomVariableInterface getVarianceForecast(double forwardVaRTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException
  	{
  		RandomVariableInterface cleanValueChange = getCleanPortfolioValueChange(forwardVaRTime);
  		//double test = cleanValueChange.squared().getAverage();
  		//for(int i=0; i<cleanValueChange.size(); i++) System.out.println(cleanValueChange.squared().get(i));
  		ConditionalExpectationEstimatorInterface condExpEstimator = getConditionalExpectationEstimator(forwardVaRTime, model);

  		RandomVariableInterface variance = cleanValueChange.squared().getConditionalExpectation(condExpEstimator).floor(0.0);
  		//for(int i=0; i<variance.size(); i++) System.out.println(variance.get(i));
  	    return variance;
  }
  	

   /**
    * Return the conditional expectation estimator suitable for this product.
    * 
    * @param forwardVaRTime The condition time.
    * @param model The model
    * @return The conditional expectation estimator suitable for this product
    * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
    */
   private ConditionalExpectationEstimatorInterface getConditionalExpectationEstimator(double forwardVaRTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
    	MonteCarloConditionalExpectationRegression condExpEstimator = new MonteCarloConditionalExpectationRegression(
	   		getRegressionBasisFunctions(forwardVaRTime, model)
		);
	    return condExpEstimator;
   }

  	/**
  	 * Provides basis funtions for the calculation of the forward variance which is the 
  	 * conditional expectation of the squared portfolio value change over the marginal period of risk 
  	 * conditional on the NPV
  	 * 
  	 * @param forwardVaRTime the time when the value at risk regression is performed
  	 * @param model The model
  	 * @return The basis functions based on the NPV of the portfolio under consideration in this class
  	 * @throws CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method. 
  	 */
  	private RandomVariableInterface[] getRegressionBasisFunctions(double forwardVaRTime,
  			                                                      LIBORModelMonteCarloSimulationInterface model
  			                                                      ) throws CalculationException {
  		// If Libor for last CF is not yet fixed
  		RandomVariableInterface NPV = portfolio.getValue(forwardVaRTime, model); 
  		double lastFixingTime = model.getLiborPeriodDiscretization().getTime(model.getLiborPeriodDiscretization().getTimeIndex(portfolio.getInitialLifeTime())-1);
    	if(forwardVaRTime < lastFixingTime && !(portfolio.getProducts()[0] instanceof Swap)){ 
  		   // to get NPV at time t
    		ConditionalExpectationEstimatorInterface condExpEstimatorLibor = getConditionalExpectationEstimatorLibor(forwardVaRTime, model);
  	       // State Variables: NPV of portfolio
  		   NPV = NPV.getConditionalExpectation(condExpEstimatorLibor);
    	}
    	
  		//RandomVariableInterface NPV = portfolio.getValue(forwardVaRTime, model);
  		ArrayList<RandomVariableInterface> basisFunctions = new ArrayList<RandomVariableInterface>();
  		
  		//for(int i=0; i<NPV.size(); i++) System.out.println(NPV.get(i));
  		
  		// Basis Functions
  		for(int orderIndex = 0; orderIndex <=polynomialOrder; orderIndex ++){
  		       basisFunctions.add(NPV.pow(orderIndex));
  	    }
  	
  		RandomVariableInterface[] finalBasisFunctions = basisFunctions.toArray(new RandomVariableInterface[basisFunctions.size()]);
  	    return finalBasisFunctions;
  	}
  	
  	
  	/**Provides a conditional expectation estimator for the calculation of the future portfolio value V(t)
  	 * 
  	 * @param forwardVaRTime the time at which the value at risk regression is performed
  	 * @param model
  	 * @return
  	 * @throws CalculationException
  	 */
    private ConditionalExpectationEstimatorInterface getConditionalExpectationEstimatorLibor(double forwardVaRTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
    	MonteCarloConditionalExpectationRegression condExpEstimator = new MonteCarloConditionalExpectationRegression(
	   		getRegressionBasisFunctionsLibor(forwardVaRTime, model)
		);
	    return condExpEstimator;
   }
    

    /** Provides basis functions based on Libor rates for the calculation of the future portfolio value V(t)
     * 
     * @param forwardVaRTime the time at which the value at risk regression is performed
     * @param model 
     * @return
     * @throws CalculationException
     */
  	private RandomVariableInterface[] getRegressionBasisFunctionsLibor(double forwardVaRTime,
              LIBORModelMonteCarloSimulationInterface model
              ) throws CalculationException {

       ArrayList<RandomVariableInterface> basisFunctions = new ArrayList<RandomVariableInterface>();
       ArrayList<RandomVariableInterface> libors = new ArrayList<RandomVariableInterface>();
       // State Variables: Libors
       int timeIndex = model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(forwardVaRTime);
       int firstLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(forwardVaRTime);
       int lastLiborIndex  = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(portfolio.getInitialLifeTime());
       for(int liborIndex = firstLiborIndex; liborIndex<lastLiborIndex; liborIndex++){
           libors.add(model.getLIBOR(timeIndex,liborIndex));
       }
       RandomVariableInterface[] finalLibors = libors.toArray(new RandomVariableInterface[libors.size()]);
       // Basis Functions
       //basisFunctions.add(new RandomVariable(1.0)); // order zero
       for(int liborIndex = 0; liborIndex <finalLibors.length; liborIndex ++){
    	   for(int orderIndex = 0; orderIndex<=2; orderIndex++){
               basisFunctions.add(finalLibors[liborIndex].pow(orderIndex));
    	   }
       }

       RandomVariableInterface[] finalBasisFunctions = basisFunctions.toArray(new RandomVariableInterface[basisFunctions.size()]);
       return finalBasisFunctions;
    }
  	

}
