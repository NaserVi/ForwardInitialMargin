package initialmargin.isdasimm.test;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

//import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModel.CalibrationItem;
import initialmargin.isdasimm.changedfinmath.LIBORModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
//import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import initialmargin.isdasimm.changedfinmath.modelplugins.AbstractLIBORCovarianceModelParametric;
import initialmargin.isdasimm.changedfinmath.modelplugins.BlendedLocalVolatilityModel;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialFormIntegrated;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.OptimizerFactoryInterface;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

public class LMMCalibrationTest{
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterParam		= new DecimalFormat(" #0.000;-#0.000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	

	private CalibrationItem createCalibrationItem(double weight, double exerciseDate, double swapPeriodLength, int numberOfPeriods, double moneyness, double targetVolatility, String targetVolatilityType, ForwardCurveInterface forwardCurve, DiscountCurveInterface discountCurve) throws CalculationException {

		double[]	fixingDates			= new double[numberOfPeriods];
		double[]	paymentDates		= new double[numberOfPeriods];
		double[]	swapTenor			= new double[numberOfPeriods + 1];

		for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
			fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
			paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex + 1) * swapPeriodLength;
			swapTenor[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
		}
		swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;

		// Swaptions swap rate
		double swaprate = moneyness + getParSwaprate(forwardCurve, discountCurve, swapTenor);

		// Set swap rates for each period
		double[] swaprates = new double[numberOfPeriods];
		Arrays.fill(swaprates, swaprate);

		/*
		 * We use Monte-Carlo calibration on implied volatility.
		 * Alternatively you may change here to Monte-Carlo valuation on price or
		 * use an analytic approximation formula, etc.
		 */
		SwaptionSimple swaptionMonteCarlo = new SwaptionSimple(swaprate, swapTenor, SwaptionSimple.ValueUnit.valueOf(targetVolatilityType));
//		double targetValuePrice = AnalyticFormulas.blackModelSwaptionValue(swaprate, targetVolatility, fixingDates[0], swaprate, getSwapAnnuity(discountCurve, swapTenor));
		return new CalibrationItem(swaptionMonteCarlo, targetVolatility, weight);
	}
/**
	 * Brute force Monte-Carlo calibration of swaptions.
	 * 
	 * @throws CalculationException
	 * @throws SolverException
	 */
	@Test
	public void testATMSwaptionCalibration() throws CalculationException, SolverException {

		final int numberOfPaths		= 1000;
		final int numberOfFactors	= 1;
		
		
		/*
		 * Calibration test
		 */
		System.out.println("Calibration to Swaptions.\n");

		DiscountCurveInterface discountCurve = (DiscountCurveInterface) DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
				// Times 
				new double[] {0,0.02739726,0.065753425,0.095890411,0.178082192,0.254794521,0.345205479,0.421917808,0.506849315,0.594520548,0.673972603,0.764383562,0.843835616,0.926027397,1.01369863,1.254794521,1.512328767,2.01369863,3.010958904,4.010958904,5.010958904,6.010958904,7.019178082,8.016438356,9.01369863,10.01369863,11.01643836,12.02191781,15.01917808,18.02465753,20.02191781,25.02739726,30.03287671,40.04109589,50.04109589},
                // Discount Factors
				new double[] {1,0.942220253,1.14628676,0.973644156,0.989291916,0.988947387,0.989030365,0.989540089,0.989760412,0.990003764,0.990397338,0.990628687,0.990878391,0.991165682,0.991574886,0.992229531,0.993347703,0.993022409,0.992927371,0.990353891,0.98534136,0.977964157,0.968209156,0.956438149,0.942562961,0.927724566,0.911915214,0.895097576,0.84499878,0.798562566,0.769568088,0.707863301,0.654037617,0.562380546,0.496026132}
				);

		ForwardCurveInterface  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
				// Fixings of the forward
				new double[] {0.504109589,1.504109589,2.509589041,3.506849315,4.506849315,5.506849315,6.509589041,7.515068493,8.512328767,9.509589041,10.51232877,11.51232877,12.51232877,13.51780822,14.51506849,15.51506849,16.51506849,17.51506849,18.52328767,19.52054795,20.51780822,21.51780822,22.52054795,23.52054795,24.5260274,25.52328767,26.52328767,27.52328767,28.52328767,29.52328767,34.52876712,39.53150685,44.53424658,49.5369863,54.54246575,59.54520548},
				// Forward Rates                                                         
				new double[] {-0.002630852,-6.82E-04,0.002757708,0.005260602,0.007848164,0.010749576,0.012628982,0.014583704,0.017103188,0.017791957,0.01917447,0.019788258,0.020269155,0.02327218,0.01577317,0.026503375,0.017980753,0.016047889,0.024898978,0.010798547,0.027070148,0.014816786,0.018220786,0.016549747,0.008028913,0.020022068,0.015134412,0.016604122,0.014386016,0.026732673,0.003643934,0.024595029,0.002432369,0.02233176,0.003397059,0.020576206},
				0.5/* tenor / period length */);

		
		/*
		 * Calibration of model volatilities
		 */
		System.out.println("Brute force Monte-Carlo calibration of model volatilities:");

		/*
		 * Create a set of calibration products.
		 */
		ArrayList<String>					calibrationItemNames	= new ArrayList<String>();
		final ArrayList<CalibrationItem>	calibrationItems		= new ArrayList<CalibrationItem>();

		double	swapPeriodLength	= 0.5;

		String[] atmExpiries = {"1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y"};
		String[] atmTenors = {"1Y",	"1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y      ",	" 2Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              "};
          				         
	    double[] atmNormalVolatilities = {0.085,0.095,0.096,0.127,0.161,	0.198,	0.266,	0.34,	0.476,	0.56,	0.611,	0.657,	0.667,	0.634,	0.592,	0.56,	0.53,	0.135,	0.143,	0.149,	0.186,	0.228,	0.269,	0.337,	0.4,	0.507,	0.575,	0.617,	0.659,	0.67,	0.635,	0.591,	0.56,	0.53,	0.183,	0.191,	0.204,	0.247,	0.28,	0.316,	0.382,	0.438,	0.526,	0.582,	0.62,	0.657,	0.666,	0.631,	0.59,	0.558,	0.525,	0.224,	0.24,	0.27,	0.304,	0.336,	0.361,	0.42,	0.467,	0.539,	0.588,	0.623,	0.654,	0.662,	0.627,	0.586,	0.553,	0.52,	0.242,	0.281,	0.299,	0.339,	0.375,	0.405,	0.449,	0.492,	0.555,	0.595,	0.625,	0.652,	0.661,	0.627,	0.585,	0.549,	0.516,	0.263,	0.303,	0.321,	0.364,	0.401,	0.425,	0.463,	0.502,	0.564,	0.598,	0.628,	0.651,	0.658,	0.622,	0.58,	0.54,	0.507,	0.276,	0.317,	0.336,	0.379,	0.415,	0.438,	0.476,	0.511,	0.567,	0.603,	0.629,	0.65,	0.655,	0.615,	0.572,	0.532,	0.497,	0.295,	0.329,	0.345,	0.392,	0.43,	0.452,	0.488,	0.522,	0.575,	0.607,	0.63,	0.65,	0.651,	0.611,	0.567,	0.524,	0.486,	0.298,	0.337,	0.351,	0.402,	0.436,	0.459,	0.496,	0.53,	0.579,	0.608,	0.629,	0.65,	0.647,	0.607,	0.562,	0.518,	0.477,	0.304,	0.342,	0.359,	0.407,	0.439,	0.466,	0.504,	0.535,	0.58,	0.609,	0.63,	0.646,	0.644,	0.605,	0.557,	0.512,	0.471,	0.329,	0.371,	0.385,	0.429,	0.457,	0.481,	0.514,	0.537,	0.568,	0.586,	0.598,	0.603,	0.596,	0.554,	0.509,	0.469,	0.433,	0.344,	0.387,	0.4,	0.44,	0.463,	0.488,	0.518,	0.541,	0.568,	0.581,	0.587,	0.585,	0.57,	0.52,	0.469,	0.428,	0.394,	0.36,	0.402,	0.412,	0.446,	0.466,	0.491,	0.518,	0.539,	0.561,	0.569,	0.572,	0.564,	0.548,	0.498,	0.446,	0.406,	0.375,	0.373,	0.41,	0.416,	0.449,	0.468,	0.491,	0.518,	0.535,	0.557,	0.564,	0.564,	0.555,	0.535,	0.48,	0.427,	0.388,	0.361};

	    		
	    

		LocalDate referenceDate = LocalDate.of(2017, Month.DECEMBER, 8); 
		BusinessdayCalendarExcludingTARGETHolidays cal = new BusinessdayCalendarExcludingTARGETHolidays();
		DayCountConvention_ACT_365 modelDC = new DayCountConvention_ACT_365();
		for(int i=0; i<atmNormalVolatilities.length; i++ ) {

			LocalDate exerciseDate = cal.createDateFromDateAndOffsetCode(referenceDate, atmExpiries[i]);
			LocalDate tenorEndDate = cal.createDateFromDateAndOffsetCode(exerciseDate, atmTenors[i]);
			double	exercise		= modelDC.getDaycountFraction(referenceDate, exerciseDate);
			double	tenor			= modelDC.getDaycountFraction(exerciseDate, tenorEndDate);

			// We consider an idealized tenor grid (alternative: adapt the model grid)
			exercise	= Math.round(exercise/0.25)*0.25;
			tenor		= Math.round(tenor/0.25)*0.25;

			if(exercise < 1.0) continue;

			int numberOfPeriods = (int)Math.round(tenor / swapPeriodLength);

			double	moneyness			= 0.0;
			double	targetVolatility	= atmNormalVolatilities[i];

			String	targetVolatilityType = "VOLATILITYNORMAL";

			double	weight = 1.0;

			calibrationItems.add(createCalibrationItem(weight, exercise, swapPeriodLength, numberOfPeriods, moneyness, targetVolatility, targetVolatilityType, forwardCurve, discountCurve));
			calibrationItemNames.add(atmExpiries[i]+"\t"+atmTenors[i]);
		}

		/*
		 * Create a simulation time discretization
		 */
		// If simulation time is below libor time, exceptions will be hard to track.
		double lastTime	= 40.0;
		double dt		= 0.1;
		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);
		final TimeDiscretizationInterface liborPeriodDiscretization = timeDiscretization;

		/*
		 * Create Brownian motions 
		 */
		final BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 31415 /* seed */);
		
		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(timeDiscretization, liborPeriodDiscretization, new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0), new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0), 0.50 / 100);
		double a = 0.0 / 20.0, b = 0.0, c = 0.25, d = 0.3 / 20.0 / 2.0;
	    LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormIntegrated(timeDiscretization, liborPeriodDiscretization, a, b, c, d, true);		
					
		LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriodDiscretization, numberOfFactors, 0.04, true);
		// Create a covariance model
		//AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretization, liborPeriodDiscretization, numberOfFactors, new double[] { 0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 0.10} );
		AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter 0.0 (that is "lognormal").
		AbstractLIBORCovarianceModelParametric covarianceModelBlended = new BlendedLocalVolatilityModel(covarianceModelParametric, 0.2, true);

		// Set model properties
		Map<String, Object> properties = new HashMap<String, Object>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose normal state space for the Euler scheme (the covariance model above carries a linear local volatility model, such that the resulting model is log-normal).
		properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).		
		Double accuracy = new Double(1E-4);	// Lower accuracy to reduce runtime of the unit test
		int maxIterations = 100;
		int numberOfThreads = 4;
		OptimizerFactoryInterface optimizerFactory = new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		double[] parameterStandardDeviation = new double[covarianceModelParametric.getParameter().length];
		double[] parameterLowerBound = new double[covarianceModelParametric.getParameter().length];
		double[] parameterUpperBound = new double[covarianceModelParametric.getParameter().length];
		Arrays.fill(parameterStandardDeviation, 0.20/100.0);
		Arrays.fill(parameterLowerBound, 0.0);
		Arrays.fill(parameterUpperBound, Double.POSITIVE_INFINITY);

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).
		Map<String, Object> calibrationParameters = new HashMap<String, Object>();
		calibrationParameters.put("accuracy", accuracy);
		calibrationParameters.put("brownianMotion", brownianMotion);
		calibrationParameters.put("optimizerFactory", optimizerFactory);
		calibrationParameters.put("parameterStep", new Double(1E-4));
		properties.put("calibrationParameters", calibrationParameters);

		long millisCalibrationStart = System.currentTimeMillis();

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModel.CalibrationItem[] calibrationItemsLMM = new LIBORMarketModel.CalibrationItem[calibrationItemNames.size()];
		for(int i=0; i<calibrationItemNames.size(); i++) calibrationItemsLMM[i] = new LIBORMarketModel.CalibrationItem(calibrationItems.get(i).calibrationProduct,calibrationItems.get(i).calibrationTargetValue,calibrationItems.get(i).calibrationWeight);
		LIBORModelInterface liborMarketModelCalibrated = new LIBORMarketModel(
				liborPeriodDiscretization,
				null,
				forwardCurve, discountCurve,
				covarianceModelBlended,
				calibrationItemsLMM,
				properties);

		long millisCalibrationEnd = System.currentTimeMillis();

		System.out.println("\nCalibrated parameters are:");
		double[] param = ((AbstractLIBORCovarianceModelParametric)((LIBORMarketModel) liborMarketModelCalibrated).getCovarianceModel()).getParameter();
		for (double p : param) System.out.println(p);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion);
		LIBORModelMonteCarloSimulationInterface simulationCalibrated = new LIBORModelMonteCarloSimulation(liborMarketModelCalibrated, process);

		System.out.println("\nValuation on calibrated model:");
		double deviationSum			= 0.0;
		double deviationSquaredSum	= 0.0;
		for (int i = 0; i < calibrationItems.size(); i++) {
			AbstractLIBORMonteCarloProduct calibrationProduct = calibrationItems.get(i).calibrationProduct;
			try {
				double valueModel = calibrationProduct.getValue(simulationCalibrated);
				double valueTarget = calibrationItems.get(i).calibrationTargetValue;
				double error = valueModel-valueTarget;
				deviationSum += error;
				deviationSquaredSum += error*error;
				System.out.println(calibrationItemNames.get(i) + "\t" + "Model: " + formatterValue.format(valueModel) + "\t Target: " + formatterValue.format(valueTarget) + "\t Deviation: " + formatterDeviation.format(valueModel-valueTarget));// + "\t" + calibrationProduct.toString());
			}
			catch(Exception e) {
			}
		}
		
		
		System.out.println("Calibration of volatilities..." + (millisCalibrationEnd-millisCalibrationStart)/1000.0);
		
		double averageDeviation = deviationSum/calibrationItems.size();
		System.out.println("Mean Deviation:" + formatterValue.format(averageDeviation));
		System.out.println("RMS Error.....:" + formatterValue.format(Math.sqrt(deviationSquaredSum/calibrationItems.size())));
		System.out.println("__________________________________________________________________________________________\n");

		Assert.assertTrue(Math.abs(averageDeviation) < 1E-2);
	}
	
	private static double getParSwaprate(ForwardCurveInterface forwardCurve, DiscountCurveInterface discountCurve, double[] swapTenor) throws CalculationException {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretization(swapTenor), new TimeDiscretization(swapTenor), forwardCurve, discountCurve);
	}
}
