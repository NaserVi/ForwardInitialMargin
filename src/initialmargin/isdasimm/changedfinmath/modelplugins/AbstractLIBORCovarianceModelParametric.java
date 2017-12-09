package initialmargin.isdasimm.changedfinmath.modelplugins;
/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;
import java.util.logging.Level;
import java.util.logging.Logger;

import initialmargin.isdasimm.changedfinmath.LIBORMarketModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.MonteCarloSimulationInterface;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.OptimizerFactoryInterface;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.OptimizerInterface;
import net.finmath.optimizer.OptimizerInterface.ObjectiveFunction;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * Base class for parametric covariance models, see also {@link AbstractLIBORCovarianceModel}.
 * 
 * Parametric models feature a parameter vector which can be inspected
 * and modified for calibration purposes.
 * 
 * The parameter vector may have zero length, which indicated that the model
 * is not calibrateable.
 * 
 * This class includes the implementation of a generic calibration algorithm.
 * If you provide an arbitrary list of calibration products, the class can return
 * a new instance where the parameters are chosen such that the (weighted) root-mean-square 
 * error of the difference of the value of the calibration products and given target
 * values is minimized.
 * 
 * @author Christian Fries
 * @date 20.05.2006
 * @date 23.02.2014
 * @version 1.1
 */
public abstract class AbstractLIBORCovarianceModelParametric extends AbstractLIBORCovarianceModel {

	private static final Logger logger = Logger.getLogger("net.finmath");

	/**
	 * Constructor consuming time discretizations, which are handled by the super class.
	 * 
	 * @param timeDiscretization The vector of simulation time discretization points.
	 * @param liborPeriodDiscretization The vector of tenor discretization points.
	 * @param numberOfFactors The number of factors to use (a factor reduction is performed)
	 */
	public AbstractLIBORCovarianceModelParametric(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors) {
		super(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
	}
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model. The parameters are usually free parameters
	 * which may be used in calibration.
	 * 
	 * @return Parameter in {@link RandomVariableInterface}-array.
	 */
	public abstract RandomVariableInterface[] getParameterAsRandomVariable();
	
	/**
	 * Get the parameters of determining this parametric
	 * covariance model. The parameters are usually free parameters
	 * which may be used in calibration.
	 * 
	 * @return Parameter in double-array.
	 */
	public double[]	getParameter() {
		// get parameters
		RandomVariableInterface[] parameterAsRandomVariable = getParameterAsRandomVariable();

		// cover case of not calibrateable models
		if(parameterAsRandomVariable == null) return null;

		// get values of deterministic random variables
		return Arrays.stream(parameterAsRandomVariable).mapToDouble(param -> param.doubleValue()).toArray();
	}
	
	/**
	 * Get the parameter identifiers of determining this parametric
	 * covariance model, in case parameters are of 
	 * instance {@link RandomVariableDifferentiableInterface}.
	 * 
	 * @return Array of parameter identifiers, null if no internal 
	 * model is calibratable or random variables are not of 
	 * instance {@link RandomVariableDifferentiableInterface}
	 * */
	public long[] getParameterID() {
		RandomVariableInterface[] parameterAsRandomVariable = getParameterAsRandomVariable();
		
		if(parameterAsRandomVariable == null || !(parameterAsRandomVariable[0] instanceof RandomVariableDifferentiableInterface)) return null;
		
		return Arrays.stream(parameterAsRandomVariable).mapToLong(param -> ((RandomVariableDifferentiableInterface) param).getID()).toArray();
	}

	@Override
	public abstract Object clone();

	/**
	 * Return an instance of this model using a new set of parameters.
	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
	 * 
	 * @param parameters The new set of parameters.
	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
	 */
	public abstract AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(double[] parameters);
	
	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, double[] calibrationTargetValues, double[] calibrationWeights) throws CalculationException {
		return getCloneCalibrated(calibrationModel, calibrationProducts, calibrationTargetValues, calibrationWeights, null);
	}

	/**
	 * Performs a generic calibration of the parametric model by
	 * trying to match a given vector of calibration product to a given vector of target values
	 * using a given vector of weights.
	 * 
	 * Optional calibration parameters may be passed using the map calibrationParameters. The keys are (<code>String</code>s):
	 * <ul>
	 * 	<li><tt>brownianMotion</tt>: Under this key an object implementing {@link net.finmath.montecarlo.BrownianMotionInterface} may be provided. If so, this Brownian motion is used to build the valuation model.</li>
	 * 	<li><tt>maxIterations</tt>: Under this key an object of type Integer may be provided specifying the maximum number of iterations.</li>
	 * 	<li><tt>accuracy</tt>: Under this key an object of type Double may be provided specifying the desired accuracy. Note that this is understood in the sense that the solver will stop if the iteration does not improve by more than this number.</li>
	 * </ul>
	 * 
	 * @param calibrationModel The LIBOR market model to be used for calibrations (specifies forward curve and tenor discretization).
	 * @param calibrationProducts The array of calibration products.
	 * @param calibrationTargetValues The array of target values.
	 * @param calibrationWeights The array of weights.
	 * @param calibrationParameters A map of type Map&lt;String, Object&gt; specifying some (optional) calibration parameters.
	 * @return A new parametric model of the same type than <code>this</code> one, but with calibrated parameters.
	 * @throws CalculationException Thrown if calibration has failed.
	 */
	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, final double[] calibrationTargetValues, double[] calibrationWeights, Map<String,Object> calibrationParameters) throws CalculationException {

		if(calibrationParameters == null) calibrationParameters = new HashMap<String,Object>();
		Integer numberOfPathsParameter	= (Integer)calibrationParameters.get("numberOfPaths");
		Integer seedParameter			= (Integer)calibrationParameters.get("seed");
		Integer maxIterationsParameter	= (Integer)calibrationParameters.get("maxIterations");
		Double	parameterStepParameter	= (Double)calibrationParameters.get("parameterStep");
		Double	accuracyParameter		= (Double)calibrationParameters.get("accuracy");
		BrownianMotionInterface brownianMotionParameter	= (BrownianMotionInterface)calibrationParameters.get("brownianMotion");

		double[] initialParameters = this.getParameter();
		double[] lowerBound = new double[initialParameters.length];
		double[] upperBound = new double[initialParameters.length];
		double[] parameterStep = new double[initialParameters.length];
		double[] zero = new double[calibrationTargetValues.length];
		Arrays.fill(lowerBound, Double.NEGATIVE_INFINITY);
		Arrays.fill(upperBound, Double.POSITIVE_INFINITY);
		Arrays.fill(parameterStep, parameterStepParameter != null ? parameterStepParameter.doubleValue() : 1E-4);
		Arrays.fill(zero, 0);

		/*
		 * We allow for 2 simultaneous calibration models.
		 * Note: In the case of a Monte-Carlo calibration, the memory requirement is that of
		 * one model with 2 times the number of paths. In the case of an analytic calibration
		 * memory requirement is not the limiting factor.
		 */
		int numberOfThreads = 2;
		OptimizerFactoryInterface optimizerFactoryParameter = (OptimizerFactoryInterface)calibrationParameters.get("optimizerFactory");

		int numberOfPaths	= numberOfPathsParameter != null ? numberOfPathsParameter.intValue() : 2000;
		int seed			= seedParameter != null ? seedParameter.intValue() : 31415;
		int maxIterations	= maxIterationsParameter != null ? maxIterationsParameter.intValue() : 400;
		double accuracy		= accuracyParameter != null ? accuracyParameter.doubleValue() : 1E-7;
		final BrownianMotionInterface brownianMotion = brownianMotionParameter != null ? brownianMotionParameter : new BrownianMotion(getTimeDiscretization(), getNumberOfFactors(), numberOfPaths, seed);
		OptimizerFactoryInterface optimizerFactory = optimizerFactoryParameter != null ? optimizerFactoryParameter : new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		int numberOfThreadsForProductValuation = 2 * Math.max(2, Runtime.getRuntime().availableProcessors());
		final ExecutorService executor = null;//Executors.newFixedThreadPool(numberOfThreadsForProductValuation);

		ObjectiveFunction calibrationError = new ObjectiveFunction() {			
			// Calculate model values for given parameters
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {

				AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = AbstractLIBORCovarianceModelParametric.this.getCloneWithModifiedParameters(parameters);

				// Create a LIBOR market model with the new covariance structure.
				LIBORMarketModelInterface model = calibrationModel.getCloneWithModifiedCovarianceModel(calibrationCovarianceModel);
				ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion);
				final LIBORModelMonteCarloSimulation liborMarketModelMonteCarloSimulation =  new LIBORModelMonteCarloSimulation(model, process);

				ArrayList<Future<Double>> valueFutures = new ArrayList<Future<Double>>(calibrationProducts.length);
				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					final int workerCalibrationProductIndex = calibrationProductIndex;
					Callable<Double> worker = new  Callable<Double>() {
						public Double call() throws SolverException {
							try {
								return calibrationProducts[workerCalibrationProductIndex].getValue((MonteCarloSimulationInterface) liborMarketModelMonteCarloSimulation);
							} catch (CalculationException e) {
								// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
								return calibrationTargetValues[workerCalibrationProductIndex];
							} catch (Exception e) {
								// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
								return calibrationTargetValues[workerCalibrationProductIndex];
							}
						}
					};
					if(executor != null) {
						Future<Double> valueFuture = executor.submit(worker);
						valueFutures.add(calibrationProductIndex, valueFuture);
					}
					else {
						FutureTask<Double> valueFutureTask = new FutureTask<Double>(worker);
						valueFutureTask.run();
						valueFutures.add(calibrationProductIndex, valueFutureTask);
					}
				}
				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					try {
						double value = valueFutures.get(calibrationProductIndex).get();
						values[calibrationProductIndex] = value;
					}
					catch (InterruptedException e) {
						throw new SolverException(e);
					} catch (ExecutionException e) {
						throw new SolverException(e);
					}
				}
			}
		};

		OptimizerInterface optimizer = optimizerFactory.getOptimizer(calibrationError, initialParameters, lowerBound, upperBound, parameterStep, calibrationTargetValues);
		try {
			optimizer.run();
		}
		catch(SolverException e) {
			throw new CalculationException(e);
		}
		finally {
			if(executor != null) {
				executor.shutdown();
			}
		}

		// Get covariance model corresponding to the best parameter set.
		double[] bestParameters = optimizer.getBestFitParameters();
		AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = this.getCloneWithModifiedParameters(bestParameters);

		// Diagnostic output
		if (logger.isLoggable(Level.FINE)) {
			logger.fine("The solver required " + optimizer.getIterations() + " iterations. The best fit parameters are:");

			String logString = "Best parameters:";
			for(int i=0; i<bestParameters.length; i++) {
				logString += "\tparameter["+i+"]: " + bestParameters[i];
			}
			logger.fine(logString);
		}

		return calibrationCovarianceModel;    	
	}

	@Override
	public String toString() {
		return "AbstractLIBORCovarianceModelParametric [getParameter()="
				+ Arrays.toString(getParameter()) + "]";
	}
}


///*
// * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
// *
// * Created on 20.05.2006
// */
//
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.HashMap;
//import java.util.Map;
//import java.util.concurrent.Callable;
//import java.util.concurrent.ExecutionException;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Executors;
//import java.util.concurrent.Future;
//import java.util.concurrent.FutureTask;
//import java.util.logging.Level;
//import java.util.logging.Logger;
//import java.util.stream.IntStream;
//
//import net.finmath.exception.CalculationException;
//import net.finmath.montecarlo.BrownianMotion;
//import net.finmath.montecarlo.BrownianMotionInterface;
//import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
//import initialmargin.isdasimm.changedfinmath.LIBORMarketModelInterface;
//import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
//import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
//import net.finmath.montecarlo.process.ProcessEulerScheme;
//import net.finmath.montecarlo.process.ProcessEulerScheme.Scheme;
//import net.finmath.optimizer.OptimizerFactoryInterface;
//import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
//import net.finmath.optimizer.OptimizerInterface;
//import net.finmath.optimizer.OptimizerInterface.ObjectiveFunction;
//import net.finmath.optimizer.OptimizerInterfaceAAD.DerivativeFunction;
//import net.finmath.optimizer.SolverException;
//import net.finmath.stochastic.RandomVariableInterface;
//import net.finmath.time.TimeDiscretizationInterface;
//
///**
// * Base class for parametric covariance models, see also {@link AbstractLIBORCovarianceModel}.
// * 
// * Parametric models feature a parameter vector which can be inspected
// * and modified for calibration purposes.
// * 
// * The parameter vector may be null, which indicated that the model
// * is not calibrateable.
// * 
// * This class includes the implementation of a generic calibration algorithm.
// * If you provide an arbitrary list of calibration products, the class can return
// * a new instance where the parameters are chosen such that the (weighted) root-mean-square 
// * error of the difference of the value of the calibration products and given target
// * values is minimized.
// * 
// * @author Christian Fries
// * @author Stefan Sedlmair
// * @date 20.10.2017
// * @version 0.1
// */
//public abstract class AbstractLIBORCovarianceModelParametric extends AbstractLIBORCovarianceModel {
//
//	public enum OptimizerDerivativeType{
//		FiniteDifferences, AdjointAlgorithmicDifferentiation, AlgorithmicDifferentiation
//	}
//	
//	public enum OptimizerSolverType{
//		Vector, Scalar
//	}
//
//	private static final Logger logger = Logger.getLogger("net.finmath");
//	private OptimizerInterface optimizer = null;
//
//	/**
//	 * Constructor consuming time discretizations, which are handled by the super class.
//	 * 
//	 * @param timeDiscretization The vector of simulation time discretization points.
//	 * @param liborPeriodDiscretization The vector of tenor discretization points.
//	 * @param numberOfFactors The number of factors to use (a factor reduction is performed)
//	 */
//	public AbstractLIBORCovarianceModelParametric(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors) {
//		super(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
//	}
//
//	/**
//	 * Get the parameters of determining this parametric
//	 * covariance model. The parameters are usually free parameters
//	 * which may be used in calibration.
//	 * 
//	 * @return Parameter in {@link RandomVariableInterface}-array.
//	 */
//	public abstract RandomVariableInterface[] getParameterAsRandomVariable();
//	
//	/**
//	 * Get the parameters of determining this parametric
//	 * covariance model. The parameters are usually free parameters
//	 * which may be used in calibration.
//	 * 
//	 * @return Parameter in double-array.
//	 */
//	public double[]	getParameter() {
//		// get parameters
//		RandomVariableInterface[] parameterAsRandomVariable = getParameterAsRandomVariable();
//
//		// cover case of not calibrateable models
//		if(parameterAsRandomVariable == null) return null;
//
//		// get values of deterministic random variables
//		return Arrays.stream(parameterAsRandomVariable).mapToDouble(param -> param.doubleValue()).toArray();
//	}
//	
//	/**
//	 * Get the parameter identifiers of determining this parametric
//	 * covariance model, in case parameters are of 
//	 * instance {@link RandomVariableDifferentiableInterface}.
//	 * 
//	 * @return Array of parameter identifiers, null if no internal 
//	 * model is calibratable or random variables are not of 
//	 * instance {@link RandomVariableDifferentiableInterface}
//	 * */
//	public long[] getParameterID() {
//		RandomVariableInterface[] parameterAsRandomVariable = getParameterAsRandomVariable();
//		
//		if(parameterAsRandomVariable == null || !(parameterAsRandomVariable[0] instanceof RandomVariableDifferentiableInterface)) return null;
//		
//		return Arrays.stream(parameterAsRandomVariable).mapToLong(param -> ((RandomVariableDifferentiableInterface) param).getID()).toArray();
//	}
//
//	@Override
//	public abstract Object clone();
//
//	/**
//	 * Return an instance of this model using a new set of parameters.
//	 * Note: To improve performance it is admissible to return the same instance of the object given that the parameters have not changed. Models should be immutable.
//	 * 
//	 * @param parameters The new set of parameters.
//	 * @return An instance of AbstractLIBORCovarianceModelParametric with modified parameters.
//	 */
//	public abstract AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(double[] parameters);
//
//	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, double[] calibrationTargetValues, double[] calibrationWeights) throws CalculationException {
//		return getCloneCalibrated(calibrationModel, calibrationProducts, calibrationTargetValues, calibrationWeights, null);
//	}
//
//	/**
//	 * Performs a generic calibration of the parametric model by
//	 * trying to match a given vector of calibration product to a given vector of target values
//	 * using a given vector of weights.
//	 * 
//	 * Optional calibration parameters may be passed using the map calibrationParameters. The keys are (<code>String</code>s):
//	 * <ul>
//	 * 	<li><tt>brownianMotion</tt>: Under this key an object implementing {@link net.finmath.montecarlo.BrownianMotionInterface} may be provided. If so, this Brownian motion is used to build the valuation model.</li>
//	 * 	<li><tt>maxIterations</tt>: Under this key an object of type Integer may be provided specifying the maximum number of iterations.</li>
//	 * 	<li><tt>accuracy</tt>: Under this key an object of type Double may be provided specifying the desired accuracy. Note that this is understood in the sense that the solver will stop if the iteration does not improve by more than this number.</li>
//	 * </ul>
//	 * 
//	 * @param calibrationModel The LIBOR market model to be used for calibrations (specifies forward curve and tenor discretization).
//	 * @param calibrationProducts The array of calibration products.
//	 * @param calibrationTargetValues The array of target values.
//	 * @param calibrationWeights The array of weights.
//	 * @param calibrationParameters A map of type Map&lt;String, Object&gt; specifying some (optional) calibration parameters.
//	 * @return A new parametric model of the same type than <code>this</code> one, but with calibrated parameters.
//	 * @throws CalculationException Thrown if calibration has failed.
//	 */
//	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, final double[] calibrationTargetValues, double[] calibrationWeights, Map<String,Object> calibrationParameters) throws CalculationException {
//		
//		double[] initialParameters = this.getParameter();
//				
//		// if nothing to calibrate return the same model
//		if(initialParameters == null) return this;
//		
//		final int numberOfCalibrationProducts 	= calibrationProducts.length;
//		final int numberOfParameters 			= initialParameters.length;
//		
//		if(numberOfCalibrationProducts != calibrationTargetValues.length) throw new IllegalArgumentException("Each calibration product has to have a target value!");
//		
//		// get calibration parameters
//		if(calibrationParameters == null) calibrationParameters = new HashMap<String,Object>();
//		
//		final int 						numberOfPaths		= (int)calibrationParameters.getOrDefault(		"numberOfPaths", 		2000);
//		final int 						seed				= (int)calibrationParameters.getOrDefault(		"seed", 				31415);
//		final int 						maxIterations		= (int)calibrationParameters.getOrDefault(		"maxIterations", 		400);
//		final int 						numberOfThreads 	= (int)calibrationParameters.getOrDefault(		"numberOfThreads", 		2);
//		final double					parameterStepValue	= (double)calibrationParameters.getOrDefault(	"parameterStep", 		1E-4);
//		final double					accuracy			= (double)calibrationParameters.getOrDefault(	"accuracy", 			1E-7);
//		final Scheme 					processScheme		= (Scheme)calibrationParameters.getOrDefault(	"scheme",  				Scheme.EULER_FUNCTIONAL);
//		final double[] 					lowerBound 			= (double[])calibrationParameters.getOrDefault(	"parameterLowerBound", 	initialzeDoubleArray(Double.NEGATIVE_INFINITY, numberOfParameters));
//		final double[] 					upperBound 			= (double[])calibrationParameters.getOrDefault(	"parameterUpperBound", 	initialzeDoubleArray(Double.POSITIVE_INFINITY, numberOfParameters));
//		final OptimizerSolverType 		solverType 			= (OptimizerSolverType) calibrationParameters.getOrDefault(		"solverType", 			OptimizerSolverType.Vector);
//		final OptimizerDerivativeType	derivativeType 		= (OptimizerDerivativeType) calibrationParameters.getOrDefault(	"derivativeType", 		OptimizerDerivativeType.FiniteDifferences);
//		final BrownianMotionInterface 	brownianMotion		= (BrownianMotionInterface)calibrationParameters.getOrDefault(	"brownianMotion", 	new BrownianMotion(getTimeDiscretization(), getNumberOfFactors(), numberOfPaths, seed));
//		final OptimizerFactoryInterface optimizerFactory 	= (OptimizerFactoryInterface)calibrationParameters.getOrDefault("optimizerFactory", new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads));
//
//		final double[] parameterStep = initialzeDoubleArray(parameterStepValue, numberOfParameters);
//
//		final ExecutorService executor = (ExecutorService)calibrationParameters.getOrDefault("executor", (numberOfThreads > 1) ? Executors.newFixedThreadPool(numberOfThreads) : null);
//		
//		DerivativeFunction calibrationErrorExtended = new DerivativeFunction() {
//			
//			// avoid calculating the products twice with the same parameters
//			double[]								currentParameters = null;
//			RandomVariableInterface[] 				currentCalibratedPrices = null;
//			AbstractLIBORCovarianceModelParametric 	currentCalibrationCovarianceModel = null;
//			
//			private void updateCalibratedPriceStorage(double[] parameters){
//				// if parameters are the same do not calculate the prices again
//				if(!Arrays.equals(currentParameters, parameters)){
//						
//					currentParameters = parameters.clone();
//					currentCalibrationCovarianceModel = AbstractLIBORCovarianceModelParametric.this.getCloneWithModifiedParameters(currentParameters);
//				
//					currentCalibratedPrices = 
//						getFutureValuesFromParameters(currentCalibrationCovarianceModel, calibrationModel, brownianMotion, calibrationProducts, 
//														numberOfCalibrationProducts, executor, calibrationTargetValues, processScheme);
//				}
//			}
//			
//			@Override
//			public void setValues(double[] parameters, double[] values) throws SolverException {
//
//				updateCalibratedPriceStorage(parameters);
//								
//				switch(solverType) {
//				case Vector:
//					// copy the prices under the current calibration directly to the values
//					double[] calibratedPrices = Arrays.stream(currentCalibratedPrices).mapToDouble(price -> price.getAverage()).toArray();
//					System.arraycopy(calibratedPrices, 0, values, 0, numberOfCalibrationProducts);
//					break;
//				case Scalar:
//					// calculate the mean-square-error
//					double errorRMS = calculateRMSError(currentCalibratedPrices, calibrationTargetValues).doubleValue();
//					System.arraycopy(new double[] {errorRMS} , 0, values, 0, 1);
//					break;
//				}		
//			}
//			
//			@Override
//			public void setDerivatives(double[] parameters, double[][] derivatives) throws SolverException {				
//				
//				updateCalibratedPriceStorage(parameters);
//				
//				// for convenience
//				RandomVariableInterface zero = calibrationModel.getRandomVariableForConstant(0.0);
//
//				switch(solverType) {
//				case Vector:
//					switch(derivativeType) {
//					case AdjointAlgorithmicDifferentiation:
//						// get parameter keys
//						long[] keys = currentCalibrationCovarianceModel.getParameterID();					
//
//						ArrayList<Future<Map<Long, RandomVariableInterface>>> derivativeFutureAAD = new ArrayList<>(numberOfCalibrationProducts);
//						
//						// loop over all leaf nodes of the operator tree
//						for(int productIndex=0; productIndex < numberOfCalibrationProducts; productIndex++) {
//							final RandomVariableInterface calibratedPrice = currentCalibratedPrices[productIndex];
//
//							// perform AAD for every leaf node
//							Callable<Map<Long, RandomVariableInterface>> worker = new Callable<Map<Long,RandomVariableInterface>>() {
//								
//								@Override
//								public Map<Long, RandomVariableInterface> call() throws Exception {
//									Map<Long, RandomVariableInterface> gradient = ((RandomVariableDifferentiableInterface) calibratedPrice).getGradient();
//									// only the averages of the derivatives are needed!
//									gradient.replaceAll((k, v) -> v.average());
//									return gradient;
//								}
//							};
//							
//							if(executor != null){
//								Future<Map<Long, RandomVariableInterface>> gradientFuture = executor.submit(worker);
//								derivativeFutureAAD.add(gradientFuture);
//							} else {
//								FutureTask<Map<Long, RandomVariableInterface>> gradientFutureTask = new FutureTask<>(worker);
//								gradientFutureTask.run();
//								derivativeFutureAAD.add(gradientFutureTask);
//							}
//						}
//						
//						// request the ids of the parameters from the calibrated model
//						for(int productIndex=0; productIndex < numberOfCalibrationProducts; productIndex++) {
//							Map<Long, RandomVariableInterface> gradient = null;
//							try { gradient = derivativeFutureAAD.get(productIndex).get();} catch (InterruptedException | ExecutionException e) {e.printStackTrace();}
//							
//							for(int parameterIndex = 0; parameterIndex < parameters.length; parameterIndex++) 
//								// do not stop the optimizer when derivative is not found. Set default to zero.
//								derivatives[parameterIndex][productIndex] = gradient.getOrDefault(keys[parameterIndex], zero).doubleValue();
//						}	
//						break;
//						default: 
//							break;
//					}
//					
//				case Scalar:				
//					RandomVariableInterface errorRMS = calculateRMSError(currentCalibratedPrices, calibrationTargetValues);
//					
//					// take gradient of the mean-square-error (here AAD should bring the most improvement!)
//					Map<Long, RandomVariableInterface> gradient = ((RandomVariableDifferentiableInterface) errorRMS).getGradient();
//										
//					// request the ids of the parameters from the calibrated model
//					long[] keys = currentCalibrationCovarianceModel.getParameterID();
//					
//					// fill in the calculated gradient in the derivative matrix
//					for(int parameterIndex = 0; parameterIndex < parameters.length; parameterIndex++) 
//						// do not stop the optimizer when derivative is not found. Set default to zero.
//						derivatives[parameterIndex][0] = gradient.getOrDefault(keys[parameterIndex], zero).getAverage();
//					break;
//				}
//			}
//		};
//
//		// define evaluating functions
//		ObjectiveFunction calibrationError = new ObjectiveFunction() {		
//
//			// Calculate model values for given parameters
//			@Override
//			public void setValues(double[] parameters, double[] values) throws SolverException {
//				calibrationErrorExtended.setValues(parameters, values);
//			}		
//		};
//		
//		// in case of a scalar solving the target value will be zero
//		double[] targetValues = null;
//		switch(solverType) {
//		case Vector:
//			targetValues = calibrationTargetValues.clone();
//			break;
//		case Scalar:
//			targetValues = new double[] {0.0};
//			break;
//		}
//
//		switch(derivativeType) {
//		case FiniteDifferences:
//			optimizer = optimizerFactory.getOptimizer(calibrationError, initialParameters, lowerBound, upperBound, parameterStep, targetValues);
//			break;
//		case AdjointAlgorithmicDifferentiation:
//		case AlgorithmicDifferentiation:
//			optimizer = optimizerFactory.getOptimizer(calibrationErrorExtended, initialParameters, lowerBound, upperBound, parameterStep, targetValues);
//			break;
//		}
//					
//		try {
//			optimizer.run();
//		}
//		catch(SolverException e) {
//			throw new CalculationException(e);
//		}
//		finally {
//			if(executor != null) {
//				executor.shutdown();
//			}
//		}
//
//		// Get covariance model corresponding to the best parameter set.
//		final double[] bestParameters = optimizer.getBestFitParameters();
//		final AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = this.getCloneWithModifiedParameters(bestParameters);
//
//		// Diagnostic output
//		if (logger.isLoggable(Level.FINE)) {
//			logger.fine("The solver required " + optimizer.getIterations() + " iterations. The best fit parameters are:");
//
//			String logString = "Best parameters:";
//			for(int i=0; i<bestParameters.length; i++) {
//				logString += "\tparameter["+i+"]: " + bestParameters[i];
//			}
//			logger.fine(logString);
//		}
//
//		return calibrationCovarianceModel;    	
//	}
//
//	
//	
//	private RandomVariableInterface[] getFutureValuesFromParameters(AbstractLIBORCovarianceModel calibrationCovarianceModel, LIBORMarketModelInterface calibrationModel, BrownianMotionInterface brownianMotion, final AbstractLIBORMonteCarloProduct[] calibrationProducts, int numberOfCalibrationProducts, ExecutorService executor, double[] calibrationTargetValues, Scheme processScheme){
//				
//		// Create a LIBOR market model with the new covariance structure.
//		final LIBORMarketModelInterface 		model 									= calibrationModel.getCloneWithModifiedCovarianceModel(calibrationCovarianceModel);
//		final ProcessEulerScheme 				process 								= new ProcessEulerScheme(brownianMotion, processScheme);
//		final LIBORModelMonteCarloSimulation 	liborMarketModelMonteCarloSimulation 	= new LIBORModelMonteCarloSimulation(model, process);
//
//		ArrayList<Future<RandomVariableInterface>> valueFutures = new ArrayList<Future<RandomVariableInterface>>(numberOfCalibrationProducts);
//		for(int calibrationProductIndex=0; calibrationProductIndex<numberOfCalibrationProducts; calibrationProductIndex++) {
//			final int workerCalibrationProductIndex = calibrationProductIndex;
//			Callable<RandomVariableInterface> worker = new  Callable<RandomVariableInterface>() {
//				public RandomVariableInterface call() throws SolverException {
//					try {
//						return calibrationProducts[workerCalibrationProductIndex].getValue(0.0, liborMarketModelMonteCarloSimulation);
//					} catch (CalculationException e) {
//						// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
//						return  calibrationModel.getRandomVariableForConstant(calibrationTargetValues[workerCalibrationProductIndex]); 
//					} catch (Exception e) {
//						// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
//						return  calibrationModel.getRandomVariableForConstant(calibrationTargetValues[workerCalibrationProductIndex]); 
//					}
//				}
//			};
//			if(executor != null) {
//				Future<RandomVariableInterface> valueFuture = executor.submit(worker);
//				valueFutures.add(calibrationProductIndex, valueFuture);
//			}
//			else {
//				FutureTask<RandomVariableInterface> valueFutureTask = new FutureTask<RandomVariableInterface>(worker);
//				valueFutureTask.run();
//				valueFutures.add(calibrationProductIndex, valueFutureTask);
//			}
//		}
//		
//		// get calculated prices	
//		RandomVariableInterface[] calibratedPrices = valueFutures.stream().map(futureValue -> {
//				RandomVariableInterface calibratedPrice = null;
//				try {
//					calibratedPrice = futureValue.get().average();
//				} catch (InterruptedException | ExecutionException e) {
//					e.printStackTrace();
//				}
//				return calibratedPrice;
//		}).toArray(RandomVariableInterface[]::new);
//		
//		return calibratedPrices;
//	}
//	
//	public OptimizerInterface getCalibrationOptimizer() {
//		return optimizer; 
//	}
//	
//	
//	@Override
//	public String toString() {
//		return "AbstractLIBORCovarianceModelParametric [getParameter()="
//				+ Arrays.toString(getParameter()) + "]";
//	}
//	
//	private static double[] initialzeDoubleArray(double value, int length) {
//		double[] array = new double[length];
//		Arrays.fill(array, value);
//		return array;
//	}
//	
//	/**
//	 * calculate the root mean square error for a give {@link RandomVariableInterface} Array with corresponding target values
//	 * <center>errorRMS = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(E[x<sub>i</sub>]-y<sub>i</sub>)<up>2</up>}</center>
//	 * 
//	 * @param model array of {@link RandomVariableInterface}s of model values
//	 * @param target double array of target values
//	 * @return deterministic {@link RandomVariableInterface} of RMS error 
//	 * */
//	private static RandomVariableInterface calculateRMSError(RandomVariableInterface[] model, double[] target){
//		int numberOfValues = model.length;
//		RandomVariableInterface errorRMS = IntStream.range(0, numberOfValues).parallel().mapToObj(
//				i -> model[i].average().sub(target[i]).squared()
//					).reduce(RandomVariableInterface::add).get().div(numberOfValues).sqrt();
//		return errorRMS;
//	}
//}
