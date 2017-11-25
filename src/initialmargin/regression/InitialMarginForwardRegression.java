package initialmargin.regression;

import java.util.ArrayList;

import initialmargin.regression.changedfinmath.products.*;

import net.finmath.exception.CalculationException;
import net.finmath.functions.NormalDistribution;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
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
	private final double MPR = 10.0/250.0; // The Marginal Period of Risk: 10 Days
	private enum Method {SIMPLE, LSQREGRESSION}; // The Method to calculate Initial Margin
	private Method method;
	private int polynomialOrder; // The order of the regression polynomial

	private LIBORModelMonteCarloSimulationInterface model;
	private Portfolio portfolio;
    //...

	public InitialMarginForwardRegression(Portfolio portfolio,
										  LIBORModelMonteCarloSimulationInterface model,
										  int polynomialOrder,
										  String method){
		this.model=model;
		this.portfolio = portfolio;
		this.polynomialOrder=polynomialOrder;
		this.method = Method.valueOf(method.toUpperCase());
	}


	/** Calculate initial margin at for a given time. 
	 * 
	 * @param evaluationTime The time at which the initial margin is calculaed.
	 * @return The initial margin of the portfolio.
	 * @throws CalculationException
	 */
	public double getInitialMargin(double evaluationTime) throws CalculationException{

		RandomVariableInterface variance;
		double initialMargin = 0;

		switch (method){	   

		case LSQREGRESSION: // Least Square Regression

			variance = getVarianceForecast(evaluationTime, model);

			double normalQuantile = NormalDistribution.inverseCumulativeDistribution(confidenceLevel);

			RandomVariableInterface initialMarginPathwise = variance.sqrt().mult(normalQuantile);

			initialMargin = initialMarginPathwise.getAverage();

			break;

		case SIMPLE: // Simple Dynamic Initial Margin

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
	
	/** Calculates the forecast of the variance of the clean portfolio value change over the marginal period of risk for a given time point.
	 * 
	 * @param forwardVaRTime The time for which the variance of the clean portfolio value change over the marginal period of risk is calculated
	 * @param model Interface implementing the Libor Market Model 
	 * @return The variance of the clean portfolio value change over the MPR
	 * @throws CalculationException
	 */
	public RandomVariableInterface getVarianceForecast(double forwardVaRTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException
	{
		RandomVariableInterface cleanValueChange = getCleanPortfolioValueChange(forwardVaRTime);

		ConditionalExpectationEstimatorInterface condExpEstimator = getConditionalExpectationEstimator(forwardVaRTime, model);

		RandomVariableInterface variance = cleanValueChange.squared().getConditionalExpectation(condExpEstimator).floor(0.0);

		return variance;
	}


	/**Calculates the clean portfolio value change, i.e. V(t+MPR)-V(t)+CF({t,t+MPR})
	 * 
	 * @param time the time t at which the marginal period of risk (MPR) starts
	 * @return The clean portfolio value change over the MPR
	 * @throws CalculationException
	 */
	public RandomVariableInterface getCleanPortfolioValueChange(double time) throws CalculationException{

		double lastFixingTime = model.getLiborPeriodDiscretization().getTime(model.getLiborPeriodDiscretization().getTimeIndex(portfolio.getInitialLifeTime())-1);

		RandomVariableInterface cashFlows = portfolio.getCF(time, time+MPR, model);

		RandomVariableInterface initialValue = portfolio.getValue(time, model);
		initialValue = initialValue.sub(cashFlows);

		if(time>0 && time < lastFixingTime && !(portfolio.getProducts()[0] instanceof Swap)) { // For swap we go forward along paths

			ConditionalExpectationEstimatorInterface condExpOperatorInitial = getConditionalExpectationEstimatorLibor(time, model);

			initialValue = initialValue.getConditionalExpectation(condExpOperatorInitial);

		}

		RandomVariableInterface finalValue = portfolio.getValue(time+MPR, model);

		if(time+MPR<lastFixingTime && !(portfolio.getProducts()[0] instanceof Swap)) {

			ConditionalExpectationEstimatorInterface condExpOperatorFinal = getConditionalExpectationEstimatorLibor(time+MPR, model);

			finalValue = finalValue.getConditionalExpectation(condExpOperatorFinal);

		}

		return finalValue.sub(initialValue);

	}

	
	/** Calculates the forecast of the variance of the clean portfolio value change over the marginal period of risk for all time points on a time discretization.
	 * 
	 * @param forwardVaRTimes The times for which the variance of the clean portfolio value change over the marginal period of risk is calculated
	 * @param model Interface implementing the Libor Market Model 
	 * @return The variance of the clean portfolio value change over the MPR
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getVarianceForecast(TimeDiscretizationInterface forwardVaRTimes, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{

		RandomVariableInterface[] VaRForecast = new RandomVariableInterface[forwardVaRTimes.getNumberOfTimes()];
		for(int timeIndex = 0; timeIndex < forwardVaRTimes.getNumberOfTimes(); timeIndex++){
			VaRForecast[timeIndex] = getVarianceForecast(forwardVaRTimes.getTime(timeIndex), model);
		}

		return VaRForecast;

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
