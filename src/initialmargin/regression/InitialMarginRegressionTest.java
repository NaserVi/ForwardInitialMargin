package initialmargin.regression;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;


import initialmargin.regression.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.regression.changedfinmath.products.Portfolio;
import initialmargin.regression.changedfinmath.products.SimpleSwap;
import initialmargin.regression.changedfinmath.products.Swap;
import initialmargin.regression.changedfinmath.products.SwapLeg;
import initialmargin.regression.changedfinmath.products.Swaption;
//import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
/*
import net.finmath.montecarlo.interestrate.products.Swap;
import net.finmath.montecarlo.interestrate.products.SwapLeg;
import net.finmath.montecarlo.interestrate.products.components.AbstractNotional;
import net.finmath.montecarlo.interestrate.products.components.Notional;
import net.finmath.montecarlo.interestrate.products.indices.AbstractIndex;
import net.finmath.montecarlo.interestrate.products.indices.LIBORIndex;
*/
import initialmargin.regression.changedfinmath.products.components.AbstractNotional;
import initialmargin.regression.changedfinmath.products.components.Notional;
import initialmargin.regression.changedfinmath.products.indices.AbstractIndex;
import initialmargin.regression.changedfinmath.products.indices.LIBORIndex;
import initialmargin.simm.changedfinmath.LIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.CurveInterface;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarInterface;

public class InitialMarginRegressionTest {
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
	final static DecimalFormat formatterIM  	= new DecimalFormat("0.00000000000");
	
	private final static int numberOfPaths		= 1000;
	private final static int numberOfFactors	= 1;
 
     public static void main(String[] args) throws CalculationException{
    	 // Create Libor market model 
     	 DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
     			                                                                                 new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
                                                                                                 new double[] {0.99 , 0.98, 0.97, 0.9, 0.7} /*discountFactors*/);
         ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
                                                                                       new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,
                                                                                       new double[] {0.05, 0.05, 0.05, 0.05, 0.05},
                                                                                       0.5/* tenor / period length */);

         LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(new RandomVariableFactory(),numberOfPaths, numberOfFactors, 
                                                                                 discountCurve,
                                                                                 forwardCurve, 0.0 /* Correlation */);
        // Another model with different volatility structure. 
        LIBORModelMonteCarloSimulationInterface model2 = createLIBORMarketModel2(5000, 2, 0.2);
     	
        // IM Portfolio Products. First test: Simple IR Swap
 		AbstractLIBORMonteCarloProduct[] products = new Swap[1];
 		products[0] = createSwap(0.0,5);
   		
 		double finalTime = 6.0;
 		double timeStep  = 1.0/20.0;
   	    // Create Portfolio of single 10y swap
  		Portfolio portfolio = new Portfolio(products, new double[]{1});
  		portfolio.setInitialLifeTime(5.0);
  		InitialMarginForwardRegression imModel = new InitialMarginForwardRegression(portfolio, model, 2 /*polynomialOrder*/, "LSQREGRESSION");
   		
  		System.out.println("Initial Margin of swap by Regression ");
   		System.out.println("Time " + "\t" + "Initial Margin");
   		for(int i = 1; i<(5.0/timeStep); i++){
   		    System.out.println(formatterTime.format(i*timeStep) + "\t " +
   		                       /*formatterIM.format(*/imModel.getInitialMargin(i*timeStep));
   		}
    
   		// Swaption
   		double     exerciseDate  = 2.0;	// Exercise date
   		double[]   fixingDates   = {2.0, 2.5, 3.0, 3.5, 4.0,4.5,5.0,5.5};   // Vector of fixing dates (must be sorted)
   		double[]   paymentDates  = {    2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};	// Vector of payment dates (same length as fixing dates)
   		double[]   swaprates     = {-0.01,-0.01,-0.01,-0.01,-0.01,-0.01,-0.01,-0.01};// Vector of strikes
   		
   		//AbstractLIBORMonteCarloProduct[] swapSimple = new AbstractLIBORMonteCarloProduct[] {new SimpleSwap(fixingDates,paymentDates,swaprates, 100.0)};
   		AbstractLIBORMonteCarloProduct[] swaption = new AbstractLIBORMonteCarloProduct[] {new Swaption(exerciseDate,fixingDates,paymentDates,swaprates,100.0)};
   		
   		Portfolio portfolio2 = new Portfolio(swaption, new double[] {1});
   		portfolio2.setInitialLifeTime(6.0);
   		InitialMarginForwardRegression IMRegression = new InitialMarginForwardRegression(portfolio2, model, 2 /*polynomialOrder*/,"LSQREGRESSION");
   		
   		System.out.println("Initial Margin of swaption by Regression ");
   		System.out.println("Time " + "\t" + "Initial Margin");
   		for(int i = 1; i<=(finalTime/timeStep); i++){
   		   System.out.println(formatterTime.format(i*timeStep) + "\t " +
   		                      formatterIM.format(IMRegression.getInitialMargin(i*timeStep)));
   		}
   		
   	}
     
     
 	public static  LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			AbstractRandomVariableFactory randomVariableFactory,
			int numberOfPaths, int numberOfFactors, DiscountCurve discountCurve, ForwardCurve forwardCurve, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 10.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

		DiscountCurveInterface appliedDiscountCurve;
		if(discountCurve==null) {
			appliedDiscountCurve = new DiscountCurveFromForwardCurve(forwardCurve);
		} else {
			appliedDiscountCurve = discountCurve;
		}
		/*
		 * Create a simulation time discretization
		 */
		double lastTime	= 10.0;
		double dt		= 0.0025;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);
        
		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		double a = 0.0 / 20.0, b = 0.0, c = 0.25, d = 0.3 / 20.0 / 2.0;
		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormIntegrated(timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);		
/*		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);
		double[][] volatilityMatrix = new double[timeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimeSteps(); timeIndex++) Arrays.fill(volatilityMatrix[timeIndex], d);
		volatilityModel = new LIBORVolatilityModelFromGivenMatrix(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, volatilityMatrix);
*/
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
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFromGivenMatrix(timeDiscretization, liborPeriodDiscretization, volatility);

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
 	
 	
 	
    // Another LMM with different volatility structure
	public static LIBORModelMonteCarloSimulationInterface createLIBORMarketModel2(int numberOfPaths, int numberOfFactors, double correlationDecayParam) throws CalculationException {

		// Create the forward curve (initial value of the LIBOR market model)
		ForwardCurve forwardCurve = ForwardCurve.createForwardCurveFromForwards(
				"forwardCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 6.0, 7.0}	/* fixings of the forward */,
				new double[] {0.05, 0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
				0.5 /*Tenor*/);

		// Create the discount curve
		DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromZeroRates(
				"discountCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 6.0, 7.0}	/* maturities */,
				new double[] {0.04, 0.04, 0.04, 0.04, 0.05, 0.05}	/* zero rates */
				);
		AnalyticModelInterface model = new AnalyticModel(new CurveInterface[] { forwardCurve , discountCurve });

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 10.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);		

		/*
		 * Create a simulation time discretization
		 */
		double lastTime	= 7.0;
		double dt		= 0.0025;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
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
		LIBORVolatilityModelFromGivenMatrix volatilityModel = new LIBORVolatilityModelFromGivenMatrix(timeDiscretization, liborPeriodDiscretization, volatility);

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

		// BlendedLocalVolatlityModel (future extension)
		//		AbstractLIBORCovarianceModel covarianceModel2 = new BlendedLocalVolatlityModel(covarianceModel, 0.00, false);

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
		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(
				liborPeriodDiscretization, model, forwardCurve, discountCurve, covarianceModel, calibrationItems, properties);

		
		ProcessEulerScheme process = new ProcessEulerScheme(
				new net.finmath.montecarlo.BrownianMotion(timeDiscretization,
						numberOfFactors, numberOfPaths, 3141 /* seed */));
		//		process.setScheme(ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
	
	
	
	
	// Product of which we calculate IM
	public static AbstractLIBORMonteCarloProduct createSwap(double evaluationTime, int maturityInYears){
	   
	    // Floating Leg
	    	// Reference Date: today
	    	Calendar calRef = Calendar.getInstance();
	    	calRef.set(Calendar.YEAR, 2017);
	    	calRef.set(Calendar.MONTH, Calendar.JANUARY);
	    	calRef.set(Calendar.DAY_OF_MONTH, 7);
	    	Date referenceDate = calRef.getTime();
	    	//Start Date
	    	Calendar calStart = calRef;
	    	int wholeYears = (int)Math.round(evaluationTime)<=evaluationTime ? (int)Math.round(evaluationTime):(int)Math.round(evaluationTime)-1;
	    	int wholeMonths= (int)Math.round((evaluationTime-wholeYears)*12.0)<=(evaluationTime-wholeYears)*12.0 ? (int)Math.round((evaluationTime-wholeYears)*12.0):(int)Math.round((evaluationTime-wholeYears)*12.0)-1;
	    	int wholeDays  = (int)((evaluationTime-wholeYears)*12.0-wholeMonths)*30;
	    	calStart.add(Calendar.YEAR, wholeYears);
	    	calStart.add(Calendar.MONTH, wholeMonths);
	    	calStart.add(Calendar.DAY_OF_MONTH, wholeDays);
			Date startDate = calStart.getTime();
			
	    	//Maturity Date
	    	Calendar calMat = calStart;
	    	calMat.add(Calendar.YEAR, maturityInYears);
	    	Date maturityDate = calMat.getTime();
		
		
		String		frequency = "semiannual";
		String		daycountConvention = "30/360";

		/*
		 * Create Monte-Carlo leg
		 */
		AbstractNotional notional = new Notional(100.0);//*(1+Math.max(Math.random(), -0.7)));
		AbstractIndex index = new LIBORIndex(0.0, 0.5);
		double spread = 0.0;
		
		ScheduleInterface schedule = ScheduleGenerator.createScheduleFromConventions(referenceDate, startDate, maturityDate, frequency, daycountConvention, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
		//ScheduleInterface schedule = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturity, frequency, daycountConvention, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
		SwapLeg leg = new SwapLeg(schedule, notional, index, spread, false /* isNotionalExchanged */);
	    
		// Fixed Leg
		
		
		String		frequencyF = "semiannual";
		String		daycountConventionF = "30/360";

		/*
		 * Create Monte-Carlo leg
		 */
		AbstractNotional notionalF = notional;
		AbstractIndex indexF = null;
		double spreadF = 0.00;
		ScheduleInterface scheduleF = ScheduleGenerator.createScheduleFromConventions(referenceDate, startDate, maturityDate, frequencyF, daycountConventionF, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), 0, 0);
		//ScheduleInterface scheduleF = ScheduleGenerator.createScheduleFromConventions(referenceDateF, spotOffsetDaysF, forwardStartPeriodF, maturityF, frequencyF, daycountConventionF, "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
		SwapLeg legF = new SwapLeg(scheduleF, notionalF, indexF, spreadF, false /* isNotionalExchanged */);

		// Swap
		return new Swap(leg,legF);

	    }
   	
    	 
    	 
     
     
}

