package initialmargin.isdasimm.test;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import initialmargin.isdasimm.changedfinmath.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModelWithTenorRefinement;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.TermStructureModelInterface;
import initialmargin.isdasimm.changedfinmath.TermStructureModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.Swap;
import initialmargin.isdasimm.changedfinmath.products.SwapLeg;
import initialmargin.isdasimm.changedfinmath.products.components.AbstractNotional;
import initialmargin.isdasimm.changedfinmath.products.components.Notional;
import initialmargin.isdasimm.changedfinmath.products.indices.AbstractIndex;
import initialmargin.isdasimm.changedfinmath.products.indices.LIBORIndex;
import initialmargin.isdasimm.products.AbstractSIMMProduct;
import initialmargin.isdasimm.products.SIMMBermudanSwaption;
import initialmargin.isdasimm.products.SIMMBermudanSwaption.ExerciseType;
import initialmargin.isdasimm.products.SIMMPortfolio;
import initialmargin.isdasimm.products.SIMMSimpleSwap;
import initialmargin.isdasimm.products.SIMMSwaption;
import initialmargin.isdasimm.products.SIMMSwaption.DeliveryType;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.SensitivityMode;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.WeightMode;
import net.finmath.analytic.model.curves.DiscountCurve;
import net.finmath.analytic.model.curves.DiscountCurveInterface;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
//import net.finmath.montecarlo.interestrate.TermStructureModelInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.interestrate.modelplugins.TermStructCovarianceModelFromLIBORCovarianceModelParametric;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretization.ShortPeriodLocation;
import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;


public class SIMMTest {
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
	final static double upperQuantile = 0.025;
	final static double lowerQuantile = 0.975;

	final static boolean isPrintAverage  = true;
	final static boolean isPrintQuantile = false;
	final static boolean isPrintPaths    = false;

	final static boolean isCalculatePortfolio = false;
	final static boolean isCalculateSwap      = false;
	final static boolean isCalculateSwaption  = false;
	final static boolean isCalculateBermudan  = true;

	public static void main(String[] args) throws CalculationException{

		/*
		 *  Create a Libor market Model
		 */

		AbstractRandomVariableFactory randomVariableFactory = createRandomVariableFactoryAAD();
		DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
				getRVAAD(new double[] {0.996 , 0.995, 0.994, 0.993, 0.98}) /*discountFactors*/);

		ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,					                                                            
				new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
				0.5/* tenor / period length */);

		LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(false,randomVariableFactory,200/*numberOfPaths*/, 1 /*numberOfFactors*/, 
				discountCurve,
				forwardCurve,0.0 /* Correlation */);


		/*
		 *  Create Products. Input for Swap
		 */
		double     startTime            = 0.0;	// Exercise date
		double     constantSwapRateSwap = 0.0;
		int        numberOfPeriodsSwap  = 10;
		double     notionalSwap        = 100;
		double[]   fixingDatesSwap     = new double[numberOfPeriodsSwap];
		double[]   paymentDatesSwap    = new double[numberOfPeriodsSwap];  	
		double[]   swapRatesSwap       = new double[numberOfPeriodsSwap];

		// Fill data
		fixingDatesSwap = IntStream.range(0, fixingDatesSwap.length).mapToDouble(i->startTime+i*0.5).toArray();
		paymentDatesSwap = IntStream.range(0, paymentDatesSwap.length).mapToDouble(i->startTime+(i+1)*0.5).toArray();
		Arrays.fill(swapRatesSwap, constantSwapRateSwap); 

		/*
		 *  Create Products. Input for (Bermudan) Swaption
		 */
		double     exerciseTime     = 5.0;	// Exercise date //5
		double     constantSwapRate = 0.025;
		int        numberOfPeriods  = 10;//20;
		double     notional         = 100;
		double[]   fixingDates     = new double[numberOfPeriods];
		double[]   paymentDates    = new double[numberOfPeriods];
		double[]   periodLength    = new double[paymentDates.length];
		double[]   periodNotionals = new double[periodLength.length];
		double[]   swapRates       = new double[numberOfPeriods];
		boolean[]  isPeriodStartDateExerciseDate = new boolean[periodLength.length]; // for Bermudan

		// Set values
		fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(i->exerciseTime+i*0.5).toArray();
		paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(i->exerciseTime+(i+1)*0.5).toArray();
		Arrays.fill(periodLength, 0.5);
		Arrays.fill(periodNotionals, notional);
		Arrays.fill(swapRates, constantSwapRate); 
		Arrays.fill(isPeriodStartDateExerciseDate, false);
		isPeriodStartDateExerciseDate[0]=true;
		isPeriodStartDateExerciseDate[2]=true;
		isPeriodStartDateExerciseDate[4]=true;
		isPeriodStartDateExerciseDate[8]=true;




		/*
		 *  Create SIMMProducts and a SIMMPortfolio
		 */	   
		AbstractSIMMProduct SIMMSwap = new SIMMSimpleSwap(fixingDatesSwap, paymentDatesSwap, swapRatesSwap, true /*isPayFix*/,notionalSwap, new String[]{"OIS", "Libor6m"}, "EUR");

		AbstractSIMMProduct SIMMSwaption = new SIMMSwaption(exerciseTime, fixingDates, paymentDates, swapRates, notional, 
				DeliveryType.Physical, new String[]{"OIS","Libor6m"}, "EUR");

		AbstractSIMMProduct SIMMBermudan = new SIMMBermudanSwaption(fixingDates, periodLength, paymentDates, periodNotionals,
				swapRates, isPeriodStartDateExerciseDate, ExerciseType.Cancelable, new String[]{"OIS", "Libor6m"}, "EUR");

		SIMMPortfolio SIMMPortfolio = new SIMMPortfolio(new AbstractSIMMProduct[]{SIMMSwap, SIMMSwaption},"EUR");

		/*
		 *  Set calculation parameters
		 */
		double finalIMTime=exerciseTime+model.getLiborPeriodDiscretization().getTimeStep(0)*numberOfPeriods;
		double timeStep = 0.1;
		double interpolationStep = 1.0;
		boolean isUseAnalyticSwapSensis = false;
		boolean isUseTimeGridAdjustment = true;
		boolean isConsiderOISSensis     = true;
        // time measurement variables
		long timeStart;
		long timeEnd;

		
		/*
		 * Perform calculations
		 */


		// Portfolio
		
		if(isCalculatePortfolio){
			RandomVariableInterface[][] valuesPortfolio = new RandomVariableInterface[4][(int)(finalIMTime/timeStep)+1];

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep;i++) valuesPortfolio[0][i] = SIMMPortfolio.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.LinearMelting, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for Portfolio, Melting: " + formatterTime.format((timeEnd-timeStart)/1000.0)+" s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep;i++) valuesPortfolio[1][i] = SIMMPortfolio.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Interpolation, WeightMode.Constant, interpolationStep, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			if(isPrintAverage){	   
				System.out.println("Expected Forward IM for Portfolio");
				System.out.println("Exact" + "\t" + "Melting " + "\t" + "Interpolation");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesPortfolio[0][i].getAverage() + "\t" + valuesPortfolio[1][i].getAverage() + "\t" + valuesPortfolio[2][i].getAverage());
				}
			}
			if(isPrintQuantile){
				System.out.println("Quantiles Forward IM for Portfolio");
				System.out.println("Upper Bound" + "\t" + "Lower Bound");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesPortfolio[0][i].getQuantile(upperQuantile) + "\t" + valuesPortfolio[0][i].getQuantile(lowerQuantile));
				}
			}
			if(isPrintPaths){
				System.out.println("Some paths of Forward IM for Portfolio");

				for(int i=0;i<finalIMTime/timeStep+1;i++){
					for(int j=0; j<10; j++){

						System.out.print(valuesPortfolio[0][i].get(j) + "\t"); 

					}
					System.out.println();
				}
			}
		}
		
		// Swap

		if(isCalculateSwap){
			RandomVariableInterface[][] valuesSwap = new RandomVariableInterface[4][(int)(finalIMTime/timeStep)+1];

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwap[0][i] = SIMMSwap.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Exact, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAP, AAD in every step, constant weights: " + formatterTime.format((timeEnd-timeStart)/1000.0)+" s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwap[1][i] = SIMMSwap.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.LinearMelting, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAP, Melting: " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwap[2][i] = SIMMSwap.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Interpolation, WeightMode.Constant, interpolationStep, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAP, Interpolation with step " + interpolationStep + ": " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwap[3][i] = SIMMSwap.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Exact, WeightMode.Stochastic, interpolationStep, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAP, AAD in every step, stochastic weights " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			if(isPrintAverage){	   
				System.out.println("Expected Forward IM for Swap");
				System.out.println("Exact, constant weights" + "\t" + "Exact, stochastic weights" + "\t" + "Melting " + "\t" + "Interpolation");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesSwap[0][i].getAverage() + "\t" + valuesSwap[3][i].getAverage() + "\t" + valuesSwap[1][i].getAverage());
				}
			}
			if(isPrintQuantile){
				System.out.println("Quantiles Forward IM for Swap");
				System.out.println("Upper Bound" + "\t" + "Lower Bound");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesSwap[0][i].getQuantile(upperQuantile) + "\t" + valuesSwap[0][i].getQuantile(lowerQuantile));
				}
			}
			if(isPrintPaths){
				System.out.println("Some paths of Forward IM for Swap");

				for(int i=0;i<finalIMTime/timeStep+1;i++){
					for(int j=0; j<10; j++){

						System.out.print(valuesSwap[0][i].get(j) + "\t"); 

					}
					System.out.println();
				}
			}
		}

		// Swaption

		if(isCalculateSwaption){
			RandomVariableInterface[][] valuesSwaption = new RandomVariableInterface[4][(int)(finalIMTime/timeStep)+1];

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[0][i] = SIMMSwaption.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Exact, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAPTION, AAD in every step, constant weights: " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[1][i] = SIMMSwaption.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.LinearMelting, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAPTION, Melting: " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[2][i] = SIMMSwaption.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Interpolation, WeightMode.Constant, interpolationStep, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for SWAPTION, Interpolation with step " + interpolationStep + ": " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			if(isPrintAverage){	   
				System.out.println("Expected Forward IM for Swaption");
				System.out.println("Exact, constant weights" + "\t" + "Melting " + "\t" + "Interpolation");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesSwaption[0][i].getAverage() + "\t" + valuesSwaption[1][i].getAverage() + "\t" + valuesSwaption[2][i].getAverage());
				}
			}
			if(isPrintQuantile){
				System.out.println("Quantiles Forward IM for Swaption");
				System.out.println("Upper Bound" + "\t" + "Lower Bound");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesSwaption[0][i].getQuantile(upperQuantile) + "\t" + valuesSwaption[0][i].getQuantile(lowerQuantile));
				}
			}
			if(isPrintPaths){
				System.out.println("Some paths of Forward IM for Swaption");

				for(int i=0;i<finalIMTime/timeStep+1;i++){
					for(int j=0; j<10; j++){

						System.out.print(valuesSwaption[0][i].get(j) + "\t"); 

					}
					System.out.println();
				}
			}

		}


		// Bermudan

		if(isCalculateBermudan){
			RandomVariableInterface[][] valuesBermudan = new RandomVariableInterface[4][(int)(finalIMTime/timeStep)+1];

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesBermudan[0][i] = SIMMBermudan.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Exact, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, true, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for BERMUDAN, AAD in every step, constant weights: " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesBermudan[1][i] = SIMMBermudan.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.LinearMelting, WeightMode.Constant, 1.0, isUseTimeGridAdjustment, true, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for BERMUDAN, Melting: " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			timeStart = System.currentTimeMillis();
			for(int i=0;i<finalIMTime/timeStep+1;i++) valuesBermudan[2][i] = SIMMBermudan.getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Interpolation, WeightMode.Constant, interpolationStep, isUseTimeGridAdjustment, true, isConsiderOISSensis);
			timeEnd = System.currentTimeMillis();

			System.out.println("Time for BERMUDAN, Interpolation with step " + interpolationStep + ": " + formatterTime.format((timeEnd-timeStart)/1000.0)+"s");

			if(isPrintAverage){	   
				System.out.println("Expected Forward IM for Bermudan");
				System.out.println("Exact" + "\t" + "Melting " + "\t" + "Interpolation");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesBermudan[0][i].getAverage() + "\t" + valuesBermudan[1][i].getAverage() + "\t" + valuesBermudan[2][i].getAverage());
				}
			}
			if(isPrintQuantile){
				System.out.println("Quantiles Forward IM for Bermudan");
				System.out.println("Upper Bound" + "\t" + "Lower Bound");
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesBermudan[0][i].getQuantile(upperQuantile) + "\t" + valuesBermudan[0][i].getQuantile(lowerQuantile));
				}
			}
			if(isPrintPaths){
				System.out.println("Some paths of Forward IM for Bermudan");

				for(int i=0;i<finalIMTime/timeStep+1;i++){
					for(int j=0; j<20; j++){

						System.out.print(valuesBermudan[0][i].get(j) + "\t"); 

					}
					System.out.println();
				}
			}


		}


	}


	public static  LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			boolean isUseTenorRefinement,
			AbstractRandomVariableFactory randomVariableFactory,
			int numberOfPaths, int numberOfFactors, DiscountCurve discountCurve, ForwardCurve forwardCurve, double correlationDecayParam) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 20.0;
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
		double lastTime	= 20.0;
		double dt		= 0.1;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

		/*
		 * Create a volatility structure v[i][j] = sigma_j(t_i)
		 */
		//double a = 0.0 / 20.0, b = 0.0, c = 0.25, d = 0.3 / 20.0 / 2.0;
		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormIntegrated(timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);		
		/*		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, a, b, c, d, false);
		double[][] volatilityMatrix = new double[timeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimeSteps(); timeIndex++) Arrays.fill(volatilityMatrix[timeIndex], d);
		volatilityModel = new LIBORVolatilityModelFromGivenMatrix(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, volatilityMatrix);
		 */
		//___________________________________________________

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
					instVolatility = 0.2 + 0.2 * Math.exp(-0.4 * timeToMaturity);

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

		BrownianMotionInterface brownianMotion = new net.finmath.montecarlo.BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, 3141 /* seed */);

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER_FUNCTIONAL);

		if(!isUseTenorRefinement){
			LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, null, forwardCurve, appliedDiscountCurve, randomVariableFactory, covarianceModel, calibrationItems, properties);

			return new LIBORModelMonteCarloSimulation(liborMarketModel, process);

		} else {

			TimeDiscretizationInterface liborPeriodDiscretizationFine = new TimeDiscretization(0.0, 40.0, 0.0625, ShortPeriodLocation.SHORT_PERIOD_AT_START);
			TimeDiscretizationInterface liborPeriodDiscretizationMedium = new TimeDiscretization(0.0, 40.0, 0.25, ShortPeriodLocation.SHORT_PERIOD_AT_START);
			TimeDiscretizationInterface liborPeriodDiscretizationCoarse = new TimeDiscretization(0.0, 40.0, 4.0, ShortPeriodLocation.SHORT_PERIOD_AT_START);
			TermStructureModelInterface liborMarketModel = new LIBORMarketModelWithTenorRefinement(
					new TimeDiscretizationInterface[] { liborPeriodDiscretizationFine, liborPeriodDiscretizationMedium, liborPeriodDiscretizationCoarse },
					new Integer[] { 4, 8, 200 },
					null,
					forwardCurve, appliedDiscountCurve, new TermStructCovarianceModelFromLIBORCovarianceModelParametric(null, covarianceModel),
					null /*calibrationItems*/, properties);
			return new TermStructureModelMonteCarloSimulation(liborMarketModel, process);
		}
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

	public static RandomVariableInterface[] getRVAAD(double[] rates){
		RandomVariableInterface[] rv = new RandomVariableInterface[rates.length];
		for(int i=0;i<rv.length;i++) rv[i]=new RandomVariableDifferentiableAAD(rates[i]);
		return rv;
	}







}

