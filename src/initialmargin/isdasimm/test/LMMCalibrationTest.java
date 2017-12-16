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

import initialmargin.isdasimm.changedfinmath.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModel.CalibrationItem;
import initialmargin.isdasimm.changedfinmath.LIBORModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.modelplugins.AbstractLIBORCovarianceModelParametric;
import initialmargin.isdasimm.changedfinmath.modelplugins.BlendedLocalVolatilityModel;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORVolatilityModel;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORVolatilityModelPiecewiseConstant;
import initialmargin.isdasimm.changedfinmath.products.SwaptionSimple;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.OptimizerFactoryInterface;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.businessdaycalendar.BusinessdayCalendar;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.daycount.DayCountConventionInterface;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

/*
 * Modified code. Source: net.finmath.montecarlo.interestrate.LIBORMarketModelCalibrationTest
 */
public class LMMCalibrationTest{
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));
	

	public static CalibrationItem createCalibrationItem(double weight, double exerciseDate, double swapPeriodLength, int numberOfPeriods, double moneyness, double targetVolatility, String targetVolatilityType, ForwardCurveInterface forwardCurve, DiscountCurveInterface discountCurve) throws CalculationException {

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
		
		// Curve Data as of December 8, 2017
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

		// Swaption ATM normal implied volatiliy as of December 8, 2017
		String[] atmExpiries = {"1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y","1M","2M","3M","6M","9M","1Y","18M","2Y","3Y","4Y","5Y","7Y","10Y","15Y","20Y","25Y","30Y"};
		String[] atmTenors =   {"1Y",	"1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","1Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y","2Y      ",	" 2Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 3Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 4Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 5Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 6Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 7Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 8Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 9Y               ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 10Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 15Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 20Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 25Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              ",	" 30Y              "};         				         
	    double[] atmNormalVolatilities = {0.00085,	0.00095,	0.00096,	0.00127,	0.00161,	0.00198,	0.00266,	0.0034,	0.00476,	0.0056,	0.00611,	0.00657,	0.00667,	0.00634,	0.00592,	0.0056,	0.0053,	0.00135,	0.00143,	0.00149,	0.00186,	0.00228,	0.00269,	0.00337,	0.004,	0.00507,	0.00575,	0.00617,	0.00659,	0.0067,	0.00635,	0.00591,	0.0056,	0.0053,	0.00183,	0.00191,	0.00204,	0.00247,	0.0028,	0.00316,	0.00382,	0.00438,	0.00526,	0.00582,	0.0062,	0.00657,	0.00666,	0.00631,	0.0059,	0.00558,	0.00525,	0.00224,	0.0024,	0.0027,	0.00304,	0.00336,	0.00361,	0.0042,	0.00467,	0.00539,	0.00588,	0.00623,	0.00654,	0.00662,	0.00627,	0.00586,	0.00553,	0.0052,	0.00242,	0.00281,	0.00299,	0.00339,	0.00375,	0.00405,	0.00449,	0.00492,	0.00555,	0.00595,	0.00625,	0.00652,	0.00661,	0.00627,	0.00585,	0.00549,	0.00516,	0.00263,	0.00303,	0.00321,	0.00364,	0.00401,	0.00425,	0.00463,	0.00502,	0.00564,	0.00598,	0.00628,	0.00651,	0.00658,	0.00622,	0.0058,	0.0054,	0.00507,	0.00276,	0.00317,	0.00336,	0.00379,	0.00415,	0.00438,	0.00476,	0.00511,	0.00567,	0.00603,	0.00629,	0.0065,	0.00655,	0.00615,	0.00572,	0.00532,	0.00497,	0.00295,	0.00329,	0.00345,	0.00392,	0.0043,	0.00452,	0.00488,	0.00522,	0.00575,	0.00607,	0.0063,	0.0065,	0.00651,	0.00611,	0.00567,	0.00524,	0.00486,	0.00298,	0.00337,	0.00351,	0.00402,	0.00436,	0.00459,	0.00496,	0.0053,	0.00579,	0.00608,	0.00629,	0.0065,	0.00647,	0.00607,	0.00562,	0.00518,	0.00477,	0.00304,	0.00342,	0.00359,	0.00407,	0.00439,	0.00466,	0.00504,	0.00535,	0.0058,	0.00609,	0.0063,	0.00646,	0.00644,	0.00605,	0.00557,	0.00512,	0.00471,	0.00329,	0.00371,	0.00385,	0.00429,	0.00457,	0.00481,	0.00514,	0.00537,	0.00568,	0.00586,	0.00598,	0.00603,	0.00596,	0.00554,	0.00509,	0.00469,	0.00433,	0.00344,	0.00387,	0.004,	0.0044,	0.00463,	0.00488,	0.00518,	0.00541,	0.00568,	0.00581,	0.00587,	0.00585,	0.0057,	0.0052,	0.00469,	0.00428,	0.00394,	0.0036,	0.00402,	0.00412,	0.00446,	0.00466,	0.00491,	0.00518,	0.00539,	0.00561,	0.00569,	0.00572,	0.00564,	0.00548,	0.00498,	0.00446,	0.00406,	0.00375,	0.00373,	0.0041,	0.00416,	0.00449,	0.00468,	0.00491,	0.00518,	0.00535,	0.00557,	0.00564,	0.00564,	0.00555,	0.00535,	0.0048,	0.00427,	0.00388,	0.00361};
	    		
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
		double lastTime	= 30.0;
		double dt		= 0.1;
		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);
		
		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 30.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		/*
		 * Create Brownian motions 
		 */
		final BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 31415 /* seed */);
		
		// Create a volatility model: Piecewise constant volatility
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(new RandomVariableFactory(),timeDiscretization, liborPeriodDiscretization, new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), new double[]{ 0.50 / 100},true);
					
		// Create a correlation model
		LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriodDiscretization, numberOfFactors, 0.04, true);
		
		// Create a covariance model
		AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter 0.0 (that is "lognormal").
		AbstractLIBORCovarianceModelParametric covarianceModelBlended = new BlendedLocalVolatilityModel(new RandomVariableFactory(),covarianceModelParametric, 0.5, true);

		// Set model properties
		Map<String, Object> properties = new HashMap<String, Object>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose normal state space for the Euler scheme (the covariance model above carries a linear local volatility model, such that the resulting model is log-normal).
		properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).		
		Double accuracy = new Double(1E-5);	// Lower accuracy to reduce runtime of the unit test
		int maxIterations = 200;
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
				new RandomVariableFactory(), // No AAD here
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

		Assert.assertTrue(Math.abs(averageDeviation) < 1E-3);
	}
	
	public static double getParSwaprate(ForwardCurveInterface forwardCurve, DiscountCurveInterface discountCurve, double[] swapTenor) throws CalculationException {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretization(swapTenor), new TimeDiscretization(swapTenor), forwardCurve, discountCurve);
	}

	public static LIBORMarketModel.CalibrationItem[] createCalibrationItems(ForwardCurveInterface forwardCurve, DiscountCurveInterface discountCurve, String[] atmExpiries, String[] atmTenors, double[] atmNormalVolatilities, LocalDate referenceDate, BusinessdayCalendar cal, DayCountConventionInterface modelDC, double swapPeriodLength) throws CalculationException{

		final ArrayList<CalibrationItem>	calibrationItems		= new ArrayList<CalibrationItem>();

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
		}
		LIBORMarketModel.CalibrationItem[] calibrationItemsLMM = new LIBORMarketModel.CalibrationItem[calibrationItems.size()];
		for(int i=0; i<calibrationItems.size(); i++) calibrationItemsLMM[i] = new LIBORMarketModel.CalibrationItem(calibrationItems.get(i).calibrationProduct,calibrationItems.get(i).calibrationTargetValue,calibrationItems.get(i).calibrationWeight);

		return calibrationItemsLMM;
	}

	public static Map<String, Object> getModelCalibrationPropertiesMap(double accuracy, double parameterStep, int maxIterations, int numberOfThreads, BrownianMotionInterface brownianMotion){	
		// Set model properties
		Map<String, Object> properties = new HashMap<String, Object>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose normal state space for the Euler scheme (the covariance model above carries a linear local volatility model, such that the resulting model is log-normal).
		properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).		
		OptimizerFactoryInterface optimizerFactory = new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).
		Map<String, Object> calibrationParameters = new HashMap<String, Object>();
		calibrationParameters.put("accuracy", new Double(accuracy));
		calibrationParameters.put("brownianMotion", brownianMotion);
		calibrationParameters.put("optimizerFactory", optimizerFactory);
		calibrationParameters.put("parameterStep", new Double(parameterStep));
		properties.put("calibrationParameters", calibrationParameters);

		return properties;
	}
	
	
	public static double[] getCalibratedParameters(LIBORModelInterface liborMarketModelCalibrated){
	   return ((AbstractLIBORCovarianceModelParametric)((LIBORMarketModel) liborMarketModelCalibrated).getCovarianceModel()).getParameter();
	}

	public static double[] getTargetValuesUnderCalibratedModel(LIBORModelInterface liborMarketModelCalibrated, BrownianMotionInterface brownianMotion, CalibrationItem[] calibrationItems){
		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion);
		LIBORModelMonteCarloSimulationInterface simulationCalibrated = new LIBORModelMonteCarloSimulation(liborMarketModelCalibrated, process);

		double[] valueModel = new double[calibrationItems.length];
		for (int i = 0; i < calibrationItems.length; i++) {
			AbstractLIBORMonteCarloProduct calibrationProduct = calibrationItems[i].calibrationProduct;
			try {
				valueModel[i] = calibrationProduct.getValue(simulationCalibrated);
			}
			catch(Exception e) {
			}
		}
		return valueModel;
	}
	
	public static double[] getCalibratedVolatilities(LIBORModelInterface liborMarketModelCalibrated){
		double[] calibratedParameters = getCalibratedParameters(liborMarketModelCalibrated);
		double[] calibratedVols = new double[calibratedParameters.length-2];
		for(int i=0;i<calibratedVols.length;i++) calibratedVols[i]=calibratedParameters[i];
		return calibratedVols;
	}
	
	public static double getCalibratedBlendingParameter(LIBORModelInterface liborMarketModelCalibrated){
		double[] calibratedParameters = getCalibratedParameters(liborMarketModelCalibrated);
		return calibratedParameters[calibratedParameters.length-1];
	}
}
