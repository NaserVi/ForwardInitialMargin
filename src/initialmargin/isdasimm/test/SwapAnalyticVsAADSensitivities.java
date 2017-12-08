package initialmargin.isdasimm.test;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import org.junit.Assert;
import org.junit.Test;

import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SimpleSwap;
import net.finmath.montecarlo.interestrate.products.components.AbstractNotional;
import net.finmath.montecarlo.interestrate.products.components.Notional;
import net.finmath.montecarlo.interestrate.products.indices.AbstractIndex;
import net.finmath.montecarlo.interestrate.products.indices.LIBORIndex;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.interestrate.products.Swap;
import net.finmath.montecarlo.interestrate.products.SwapLeg;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;

public class SwapAnalyticVsAADSensitivities {
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
	final static DecimalFormat formatterSensi 	= new DecimalFormat("0.000000000");

	public  Map<Long,RandomVariableInterface> gradient;

	@Test
	public void testAADVsAnalyticSensis() throws CalculationException{ 

		// Create a Libor market Model
		AbstractRandomVariableFactory randomVariableFactory = createRandomVariableFactoryAAD();
		DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
				new double[] {0.996 , 0.995, 0.994, 0.993, 0.98} /*discountFactors*/);

		ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,					                                                            
				new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
				0.5/* tenor / period length */);

		LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(randomVariableFactory,5000/*numberOfPaths*/, 1 /*numberOfFactors*/, 
				discountCurve,
				forwardCurve,0.0 /* Correlation */);


		// ------------------------------------------------------------------------------------------------------------
		// Create Swap
		// ------------------------------------------------------------------------------------------------------------

		double     startingTime     = 0.0;
		double     constantSwapRate = 0.02;
		int        numberOfPeriods  = 8;
		double     periodLength     = 0.5;
		//double     notional         = 100;

		double[]   fixingDates     = new double[numberOfPeriods];
		double[]   paymentDates    = new double[numberOfPeriods];
		double[]   swapRates       = new double[numberOfPeriods];

		// Fill Arrays
		fixingDates  = IntStream.range(0, fixingDates.length).mapToDouble(i->startingTime+i*periodLength).toArray();
		paymentDates = IntStream.range(0,paymentDates.length).mapToDouble(i->startingTime+(i+1)*periodLength).toArray();
		Arrays.fill(swapRates, constantSwapRate); 

		// Create Products
		AbstractLIBORMonteCarloProduct simpleSwap = new SimpleSwap(fixingDates,paymentDates,swapRates); //Notional 1


		// ------------------------------------------------------------------------------------------------------------
		// Compare Sensis
		// ------------------------------------------------------------------------------------------------------------

		// Calculate forward sensitivities
		double timeStep = 0.1;
		double relError = 0;
		int    counter=0;

		for(int timeIndex=0;timeIndex<(int)(paymentDates[paymentDates.length-1]/timeStep);timeIndex++){
			double time = timeStep*timeIndex;
			RandomVariableInterface[] sensisAAD = getAADSwapLiborSensitivities(time,simpleSwap,model);
			RandomVariableInterface[] sensisANA = getAnalyticSwapLiborSensitivities(time,periodLength,fixingDates,model);
			System.out.println("Time" + "\t" + "LiborIndex" + "\t" + "SensisAAD" + "\t" + "SensisAnalytic");
			for(int liborIndex=0;liborIndex<sensisANA.length;liborIndex++){
				System.out.println(formatterTime.format(timeIndex*timeStep) + "\t" + liborIndex + "\t" +
						formatterSensi.format(sensisAAD[liborIndex].getAverage()) + "\t" +
						formatterSensi.format(sensisANA[liborIndex].getAverage()));

				if(sensisANA[liborIndex].getAverage()!=0){
					//l2error += Math.sqrt(sensisAAD[liborIndex].sub(sensisANA[liborIndex]).squared().getAverage()*1000);
					relError +=Math.abs(sensisAAD[liborIndex].getAverage()-sensisANA[liborIndex].getAverage())/Math.abs(sensisANA[liborIndex].getAverage());
					counter++;
				}
			} 

		}
		System.out.println("Rel. Error " + relError/counter);
		Assert.assertTrue(relError/counter<0.01);
	}


	/** Calculate dV/dL
	 * 
	 * @param evaluationTime
	 * @param periodLength
	 * @param paymentDates
	 * @param model
	 * @return
	 * @throws CalculationException 
	 */
	public  RandomVariableInterface[] getAnalyticSwapLiborSensitivities(double evaluationTime, 
			double  periodLength, 
			double[] fixingDates,
			LIBORModelMonteCarloSimulationInterface model) throws CalculationException{

		MonteCarloConditionalExpectationRegression cOperator = getConditionalExpectationOperator(evaluationTime,model);
		// Calculate forward sensitivities
		int periodIndex = new TimeDiscretization(fixingDates).getTimeIndexNearestLessOrEqual(evaluationTime); 
		periodIndex = periodIndex < 0 ? 0 : periodIndex;
		int firstLiborIndex = fixingDates[0] > evaluationTime ? model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(fixingDates[0]):model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime,model);
		int numberOfSensis = evaluationTime == getNextLiborTime(evaluationTime,model) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
		RandomVariableInterface[] sensis = new RandomVariableInterface[numberOfSensis]; Arrays.fill(sensis, new RandomVariable(0.0));
		int currentLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		RandomVariableInterface numeraireAtEval = model.getNumeraire(evaluationTime);

		for(int liborIndex=currentLiborIndex;liborIndex<numberOfSensis+currentLiborIndex;liborIndex++){
			int i = liborIndex < firstLiborIndex ? 0 : liborIndex-firstLiborIndex+1;
			if(!(i>fixingDates.length-periodIndex || i==0) ){ //fixingDates[i-1]+periodLength<evaluationTime		
				// Actual Sensitivity Calculation: dV/dL = P(T,t)*periodLength
				RandomVariableInterface numeraireAtPayment = model.getNumeraire(fixingDates[periodIndex+i-1]+periodLength);
				sensis[liborIndex-currentLiborIndex]=numeraireAtEval.div(numeraireAtPayment).mult(periodLength).getConditionalExpectation(cOperator);			

			}
		}
		return sensis;
	}

	/**Calculates the row vector dV/dL
	 * 
	 * @param evaluationTime The time at which the forward sensistivity dVdL is calculated
	 * @return The forward sensisivity dVdL (as a row vector)
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getAADSwapLiborSensitivities(double evaluationTime,
			AbstractLIBORMonteCarloProduct product,
			LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		if(this.gradient==null) {
			RandomVariableDifferentiableInterface value = (RandomVariableDifferentiableInterface) product.getValue(0.0, model);
			this.gradient = value.getGradient();
		}

		MonteCarloConditionalExpectationRegression cOperator = getConditionalExpectationOperator(evaluationTime, model);
		RandomVariableInterface numeraireAtEval = model.getNumeraire(evaluationTime);

		// Calculate forward sensitivities
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime,model);
		int numberOfSensis = evaluationTime == getNextLiborTime(evaluationTime,model) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
		RandomVariableInterface[] valueLiborSensitivities = new RandomVariableInterface[numberOfSensis];// exclude last libor
		int timeIndexAtEval = model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);

		// Set all entries of dVdL
		// Set dVdL for last libor which is already fixed (if applicable)
		int lastLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		int indexAdjustment = 0;
		if(numberOfSensis!=numberOfRemainingLibors){
			indexAdjustment = 1;
			double lastLiborTime = model.getLiborPeriodDiscretization().getTime(lastLiborIndex);
			RandomVariableInterface lastLibor = model.getLIBOR(model.getTimeDiscretization().getTimeIndex(lastLiborTime), lastLiborIndex);
			RandomVariableInterface dVdL = gradient.get(((RandomVariableDifferentiableInterface)lastLibor).getID());
			valueLiborSensitivities[0] = dVdL.mult(numeraireAtEval);
		}

		for(int liborIndex=lastLiborIndex+indexAdjustment;liborIndex<model.getNumberOfLibors(); liborIndex++){
			RandomVariableInterface liborAtTimeIndex = model.getLIBOR(timeIndexAtEval, liborIndex);
			RandomVariableInterface dVdL = gradient.get(((RandomVariableDifferentiableInterface)liborAtTimeIndex).getID());
			valueLiborSensitivities[liborIndex-lastLiborIndex] = dVdL==null ? new RandomVariable(0.0) : dVdL.mult(numeraireAtEval).getConditionalExpectation(cOperator);
		}
		return valueLiborSensitivities;
	}


	private static int getNumberOfRemainingLibors(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getNumberOfLibors()-nextLiborIndex;
	}

	private static double getNextLiborTime(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex);
	}


	public static  LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			AbstractRandomVariableFactory randomVariableFactory,
			int numberOfPaths, int numberOfFactors, DiscountCurve discountCurve, ForwardCurve forwardCurve, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 12.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		DiscountCurveInterface appliedDiscountCurve;
		if(discountCurve==null) {
			appliedDiscountCurve = (DiscountCurveInterface) new DiscountCurveFromForwardCurve(forwardCurve);
		} else {
			appliedDiscountCurve = discountCurve;
		}
		/*
		 * Create a simulation time discretization
		 */
		double lastTime	= 12.0;
		double dt		= 0.1;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		//double a = 0.0 / 20.0, b = 0.0, c = 0.25, d = 0.3 / 20.0 / 2.0;


		double[][] volatility = new double[timeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
		for (int timeIndex = 0; timeIndex < volatility.length; timeIndex++) {
			for (int liborIndex = 0; liborIndex < volatility[timeIndex].length; liborIndex++) {
				// Create a very simple volatility model here
				double time = timeDiscretization.getTime(timeIndex);
				double maturity = liborPeriodDiscretization.getTime(liborIndex);
				double timeToMaturity = maturity - time;

				double instVolatility;
				if(timeToMaturity <= 0)
					instVolatility = 0;				// This forward rate is already fixed, no volatility
				else
					instVolatility = 0.3 + 0.2 * Math.exp(-0.25 * timeToMaturity);

				// Store
				volatility[timeIndex][liborIndex] = instVolatility;
			}
		}
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFromGivenMatrix(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, volatility);



		//___________________________________________________

		/*
		 * Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(
				timeDiscretization, liborPeriodDiscretization, numberOfFactors,
				correlationDecayParam);


		/*
		 * Combine volatility model and correlation model to a covariance model
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization,
						liborPeriodDiscretization, volatilityModel, correlationModel);

		// Set model properties
		Map<String, String> properties = new HashMap<String, String>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose log normal model
		properties.put("stateSpace", LIBORMarketModel.StateSpace.LOGNORMAL.name());

		// Empty array of calibration items - hence, model will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */

		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, null, forwardCurve, appliedDiscountCurve, randomVariableFactory, covarianceModel, calibrationItems, properties);

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER_FUNCTIONAL);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}


	public static AbstractLIBORMonteCarloProduct[] createSwaps(String[] maturities){
		AbstractLIBORMonteCarloProduct[] swaps = new AbstractLIBORMonteCarloProduct[maturities.length];
		// 1) Create Portfolio of swaps -------------------------------------------------------------------------------
		for(int swapIndex = 0; swapIndex < maturities.length; swapIndex++){
			// Floating Leg
			LocalDate	referenceDate = LocalDate.of(2017, 8, 12);
			int			spotOffsetDays = 0;
			String		forwardStartPeriod = "0D";
			String		maturity = maturities[swapIndex];
			String		frequency = "semiannual";
			String		daycountConvention = "30/360";

			/*
			 * Create Monte-Carlo leg
			 */
			AbstractNotional notional = new Notional(100.0);//*(1+Math.max(Math.random(), -0.7)));
			AbstractIndex index = new LIBORIndex(0.0, 0.5);
			double spread = 0.0;
			ScheduleInterface schedule = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturity, frequency, daycountConvention, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
			SwapLeg leg = new SwapLeg(schedule, notional, index, spread, false /* isNotionalExchanged */);

			// Fixed Leg
			LocalDate	referenceDateF = LocalDate.of(2017, 8, 12);
			int			spotOffsetDaysF = 0;
			String		forwardStartPeriodF = "0D";
			String		maturityF = maturities[swapIndex];
			String		frequencyF = "semiannual";
			String		daycountConventionF = "30/360";

			/*
			 * Create Monte-Carlo leg
			 */
			AbstractNotional notionalF = notional;
			AbstractIndex indexF = null;
			double spreadF = 0.00;
			ScheduleInterface scheduleF = ScheduleGenerator.createScheduleFromConventions(referenceDateF, spotOffsetDaysF, forwardStartPeriodF, maturityF, frequencyF, daycountConventionF, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
			SwapLeg legF = new SwapLeg(scheduleF, notionalF, indexF, spreadF, false /* isNotionalExchanged */);

			// Swap
			AbstractLIBORMonteCarloProduct swap = new Swap(leg,legF);
			swaps[swapIndex]=swap;
		}
		return swaps;
	}

	public static AbstractRandomVariableFactory createRandomVariableFactoryAAD(){
		Map<String, Object> properties = new HashMap<String, Object>();
		properties.put("isGradientRetainsLeafNodesOnly", new Boolean(false));
		return new RandomVariableDifferentiableAADFactory(new RandomVariableFactory(), properties);
	}

	private static MonteCarloConditionalExpectationRegression getConditionalExpectationOperator(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{

		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
		RandomVariableInterface[] regressor = new RandomVariableInterface[2];
		regressor[0]= model.getLIBOR(evaluationTime, evaluationTime,evaluationTime+model.getLiborPeriodDiscretization().getTimeStep(0));
		regressor[1]= model.getLIBOR(evaluationTime, evaluationTime, model.getLiborPeriodDiscretization().getTime(model.getNumberOfLibors()-1));
		ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(regressor, 2);
		return new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));

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
