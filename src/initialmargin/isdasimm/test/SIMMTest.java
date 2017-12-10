package initialmargin.isdasimm.test;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import initialmargin.isdasimm.changedfinmath.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.modelplugins.AbstractLIBORCovarianceModelParametric;
import initialmargin.isdasimm.changedfinmath.modelplugins.BlendedLocalVolatilityModel;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORVolatilityModel;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import initialmargin.isdasimm.changedfinmath.modelplugins.LIBORVolatilityModelPiecewiseConstant;
import initialmargin.isdasimm.products.AbstractSIMMProduct;
import initialmargin.isdasimm.products.SIMMBermudanSwaption;
import initialmargin.isdasimm.products.SIMMBermudanSwaption.ExerciseType;
import initialmargin.isdasimm.products.SIMMPortfolio;
import initialmargin.isdasimm.products.SIMMSimpleSwap;
import initialmargin.isdasimm.products.SIMMSwaption;
import initialmargin.isdasimm.products.SIMMSwaption.DeliveryType;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.SensitivityMode;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.WeightMode;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;


public class SIMMTest {
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
	final static double upperQuantile = 0.025;
	final static double lowerQuantile = 0.975;

	final static boolean isPrintAverage  = true;
	final static boolean isPrintQuantile = false;
	final static boolean isPrintPaths    = false;

	final static boolean isCalculatePortfolio = false;
	final static boolean isCalculateSwap      = true;
	final static boolean isCalculateSwaption  = false;
	final static boolean isCalculateBermudan  = false;

	public static void main(String[] args) throws CalculationException{

		/*
		 *  Create a Libor market Model
		 */

		AbstractRandomVariableFactory randomVariableFactory = createRandomVariableFactoryAAD();

		// Curve Data as of December 8, 2017
		DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
				// Times 
				new double[] {0,0.02739726,0.065753425,0.095890411,0.178082192,0.254794521,0.345205479,0.421917808,0.506849315,0.594520548,0.673972603,0.764383562,0.843835616,0.926027397,1.01369863,1.254794521,1.512328767,2.01369863,3.010958904,4.010958904,5.010958904,6.010958904,7.019178082,8.016438356,9.01369863,10.01369863,11.01643836,12.02191781,15.01917808,18.02465753,20.02191781,25.02739726,30.03287671,40.04109589,50.04109589},
                // Discount Factors
				new double[] {1,0.942220253,1.14628676,0.973644156,0.989291916,0.988947387,0.989030365,0.989540089,0.989760412,0.990003764,0.990397338,0.990628687,0.990878391,0.991165682,0.991574886,0.992229531,0.993347703,0.993022409,0.992927371,0.990353891,0.98534136,0.977964157,0.968209156,0.956438149,0.942562961,0.927724566,0.911915214,0.895097576,0.84499878,0.798562566,0.769568088,0.707863301,0.654037617,0.562380546,0.496026132}
				);

		ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
				// Fixings of the forward
				new double[] {0.504109589,1.504109589,2.509589041,3.506849315,4.506849315,5.506849315,6.509589041,7.515068493,8.512328767,9.509589041,10.51232877,11.51232877,12.51232877,13.51780822,14.51506849,15.51506849,16.51506849,17.51506849,18.52328767,19.52054795,20.51780822,21.51780822,22.52054795,23.52054795,24.5260274,25.52328767,26.52328767,27.52328767,28.52328767,29.52328767,34.52876712,39.53150685,44.53424658,49.5369863,54.54246575,59.54520548},
				// Forward Rates                                                         
				new double[] {-0.002630852,-6.82E-04,0.002757708,0.005260602,0.007848164,0.010749576,0.012628982,0.014583704,0.017103188,0.017791957,0.01917447,0.019788258,0.020269155,0.02327218,0.01577317,0.026503375,0.017980753,0.016047889,0.024898978,0.010798547,0.027070148,0.014816786,0.018220786,0.016549747,0.008028913,0.020022068,0.015134412,0.016604122,0.014386016,0.026732673,0.003643934,0.024595029,0.002432369,0.02233176,0.003397059,0.020576206},
				0.5/* tenor / period length */);
		
//		DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
//				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
//				new double[] {0.996 , 0.995, 0.994, 0.993, 0.98} /*discountFactors*/);
//
//		ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
//				new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,					                                                            
//				new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
//				0.5/* tenor / period length */);
		
	
		LIBORModelMonteCarloSimulationInterface modell = createLIBORMarketModel(false,randomVariableFactory,100/*numberOfPaths*/, 1 /*numberOfFactors*/, 
				discountCurve,
				forwardCurve);
		
		
		LIBORModelMonteCarloSimulationInterface model = getZeroVolatilityModel(modell);
		
		
		/*
		 *  Create Products. Input for Swap
		 */
		double     startTime            = 0.0;	// Exercise date
		double     constantSwapRateSwap = 0.013;
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
		double     constantSwapRate = 0.9;
		int        numberOfPeriods  = 10;
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

	
	public static  LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(boolean isUseTenorRefinement,
										AbstractRandomVariableFactory randomVariableFactory,
										int numberOfPaths, int numberOfFactors, DiscountCurve discountCurve, ForwardCurve forwardCurve) throws CalculationException {

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
		
		// Create a volatility model: Piecewise constant volatility calibrated to Swaption Normal implied volatility of December 8, 2017
		double[] volatility = new double[]{
		0.0035380523915104246,
		0.004191317634739811,
		0.008841173374561527,
		0.010367178689341235,
		0.009683514692001837,
		0.00881065062410322,
		0.005,
		0.005553622644567465,
		0.006240047498020553,
		0.008993528127078064,
		0.009894813615533201,
		0.009632854002962725,
		0.009785740680837035,
		0.0037906111648575865,
		0.014616121514330995,
		0.011590302354861921,
		0.012136753578600236,
		0.009878601075748226,
		0.008283683236194047,
		0.01158663971536579,
		0.011596322104781735,
		0.010557210170556731,
		0.011936780200631499,
		0.007661672888823457,
		0.003682948768971966,
		0.004960044546431093,
		0.01198892262469119,
		0.007086263635340348,
		0.01173222179819021,
		0.008696172864293873,
		0.005,
		0.007599201382279569,
		0.007148972951937137,
		0.006788680889933558,
		0.011140027803944855,
		0.005,
		0.005};
		
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(randomVariableFactory,timeDiscretization, liborPeriodDiscretization, new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), /*new double[]{0.01}*/volatility,false);
				
//		double[][] volatility = new double[timeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
//		for (int timeIndex = 0; timeIndex < volatility.length; timeIndex++) {
//			for (int liborIndex = 0; liborIndex < volatility[timeIndex].length; liborIndex++) {
//				// Create a very simple volatility model here
//				double time = timeDiscretization.getTime(timeIndex);
//				double maturity = liborPeriodDiscretization.getTime(liborIndex);
//				double timeToMaturity = maturity - time;
//
//				double instVolatility;
//				if(timeToMaturity <= 0)
//					instVolatility = 0;				// This forward rate is already fixed, no volatility
//				else
//					instVolatility = (0.2 + 0.2 * Math.exp(-0.4 * timeToMaturity))/25.0; // / 50 if statespace = normal
//
//				// Store
//				volatility[timeIndex][liborIndex] = instVolatility;
//			}
//		}
//		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFromGivenMatrix(randomVariableFactory, timeDiscretization, liborPeriodDiscretization, volatility);
//		
		// Create a correlation model
		LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriodDiscretization, numberOfFactors, 0.04 /*correlationDecayParameter*/, false);
		
		// Create a covariance model
		AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter 0.0 (that is "lognormal").
		double displacementParameter = 0.5880313623110442;
		AbstractLIBORCovarianceModelParametric covarianceModelBlended = new BlendedLocalVolatilityModel(randomVariableFactory,covarianceModelParametric, displacementParameter, false);

		// Set model properties
		Map<String, Object> properties = new HashMap<String, Object>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());

		// Choose normal state space for the Euler scheme (the covariance model above carries a linear local volatility model, such that the resulting model is log-normal).
		properties.put("stateSpace", LIBORMarketModel.StateSpace.NORMAL.name());

		// Empty array of calibration items - hence, model will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */

		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER_FUNCTIONAL);

		LIBORMarketModelInterface liborMarketModel = new LIBORMarketModel(liborPeriodDiscretization, null, forwardCurve, discountCurve, randomVariableFactory, covarianceModelBlended, calibrationItems, properties);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
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
	
	
	public static LIBORModelMonteCarloSimulationInterface getZeroVolatilityModel(LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		
		// Set brownian motion with one path
		BrownianMotionInterface originalBM = model.getBrownianMotion();
		BrownianMotionInterface brownianMotion = new BrownianMotion(originalBM.getTimeDiscretization(), originalBM.getNumberOfFactors(), 1 /* numberOfPaths */, 3141);
        
		// Get process
		ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion, ProcessEulerScheme.Scheme.EULER_FUNCTIONAL);
		
		// Create zero volatility model
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(createRandomVariableFactoryAAD(),model.getTimeDiscretization(), model.getLiborPeriodDiscretization(), new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), new TimeDiscretization(0.00, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0), new double[]{0.0}/*volatility*/,false);

		//double[][] volatility = new double[model.getTimeDiscretization().getNumberOfTimeSteps()][model.getLiborPeriodDiscretization().getNumberOfTimeSteps()];
		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFromGivenMatrix(new RandomVariableFactory(), model.getTimeDiscretization(), model.getLiborPeriodDiscretization(), volatility);
		
		//Create a correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(model.getTimeDiscretization(), model.getLiborPeriodDiscretization(), model.getNumberOfFactors(),0);
				
		//Combine volatility model and correlation model to a covariance model
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(model.getTimeDiscretization(),
						model.getLiborPeriodDiscretization(), volatilityModel, correlationModel);
		
		double displacementParameter = 0.5880313623110442;
		AbstractLIBORCovarianceModelParametric covarianceModelBlended = new BlendedLocalVolatilityModel(createRandomVariableFactoryAAD(),covarianceModel, displacementParameter, false);


		Map<String, Object> dataModified = new HashMap<>();
		dataModified.put("covarianceModel", covarianceModelBlended);
		return new LIBORModelMonteCarloSimulation((LIBORModelInterface)model.getModel().getCloneWithModifiedData(dataModified),process);
		
	}
	
	

}

