package initialmargin.isdasimm.test;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.stream.IntStream;

import net.finmath.exception.CalculationException;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.products.AbstractSIMMProduct;
import initialmargin.isdasimm.products.SIMMBermudanSwaption;
import initialmargin.isdasimm.products.SIMMSimpleSwap;
import initialmargin.isdasimm.products.SIMMSwaption;
import initialmargin.isdasimm.products.SIMMBermudanSwaption.ExerciseType;
import initialmargin.isdasimm.products.SIMMSwaption.DeliveryType;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.SensitivityMode;
import initialmargin.isdasimm.sensitivity.AbstractSIMMSensitivityCalculation.WeightMode;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.stochastic.RandomVariableInterface;


public class SensitivityApproximationTest {
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");

	// Model Paths 
	final static int numberOfPaths = 10;
	final static double notional = 100;

	public static enum TestProductType {Swaps, Swaptions, BermudanCallable, BermudanCancelable};

	public static void main(String[] args) throws CalculationException{

		/*
		 *  Create a Libor market Model
		 */

		AbstractRandomVariableFactory randomVariableFactory = SIMMTest.createRandomVariableFactoryAAD();

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


		LIBORModelMonteCarloSimulationInterface model = SIMMTest.createLIBORMarketModel(false,randomVariableFactory,numberOfPaths, 1 /*numberOfFactors*/, 
				discountCurve,
				forwardCurve);

		double[] exerciseDates = null;
		int[] numberOfPeriods = null;

		// Select TestProducts
		TestProductType testProductType = TestProductType.BermudanCancelable;

		switch(testProductType){
		case Swaps: 
			exerciseDates = new double[] {0.0};
			numberOfPeriods = new int[] {10,20,30,40};
			break;
		case Swaptions: 
			exerciseDates = new double[] {5.0, 10.0};
			numberOfPeriods = new int[]  {4, 8, 10, 12, 16};
		case BermudanCallable: 
		case BermudanCancelable:
			exerciseDates = new double[] {5.0, 8.0, 10.0};
			numberOfPeriods = new int[]  {20};
		}


		AbstractSIMMProduct[] products = createProducts(testProductType, exerciseDates, numberOfPeriods, forwardCurve, discountCurve);

		// Test Swaption Exact Vs Melting
		testSIMMProductApproximation(products, exerciseDates, numberOfPeriods, forwardCurve, discountCurve, model);

	}

	public static void testSIMMProductApproximation(AbstractSIMMProduct[] product, double[] exerciseDates, int[] numberOfPeriods, ForwardCurve forwardCurve, DiscountCurve discountCurve, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{

		double  timeStep = 0.1;
		boolean isUseAnalyticSwapSensis = false;
		boolean isUseTimeGridAdjustment = true;
		boolean isConsiderOISSensis     = true;

		System.out.println("Exercise in Y" + "\t" + "NumberSwapPeriods" + "\t" + "Time Exact" + "\t" + "Time Melting" + "\t" + "Time Interpolation" + "\t" + "Mean Deviation in % of IM Sum Melting" + "\t" + "Mean Deviation in % of IM Sum Interpolation" + "\t" + "RMSE Melting" + "\t" + "RMSE Interpolation");

		int productIndex = 0;
		for(int exerciseIndex = 0; exerciseIndex < exerciseDates.length; exerciseIndex++){
			for(int swapPeriodsIndex = 0; swapPeriodsIndex < numberOfPeriods.length; swapPeriodsIndex++){

				double finalIMTime=exerciseDates[exerciseIndex]+model.getLiborPeriodDiscretization().getTimeStep(0)*numberOfPeriods[swapPeriodsIndex];
				/*
				 * Calculate IM
				 */
				RandomVariableInterface[][] valuesSwaption = new RandomVariableInterface[3][(int)(finalIMTime/timeStep)+1];

				// 1) Exact
				long timeStart = System.currentTimeMillis();
				for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[0][i] = product[productIndex].getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Exact, WeightMode.Stochastic, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
				long timeEnd = System.currentTimeMillis();

				// 2) Melting
				long timeStartMelting = System.currentTimeMillis();
				for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[1][i] = product[productIndex].getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.LinearMelting, WeightMode.Stochastic, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
				long timeEndMelting = System.currentTimeMillis();

				// 3) Interpolation
				long timeStartInterpolation = System.currentTimeMillis();
				for(int i=0;i<finalIMTime/timeStep+1;i++) valuesSwaption[2][i] = product[productIndex].getInitialMargin(i*timeStep, model, "EUR", SensitivityMode.Interpolation, WeightMode.Stochastic, 1.0, isUseTimeGridAdjustment, isUseAnalyticSwapSensis, isConsiderOISSensis);
				long timeEndInterpolation = System.currentTimeMillis();

				double deviationSumMelting = 0;
				double deviationSumInterpolation = 0;
				double sumIM = 0; // The sum of IM values (always positive)
				double sumSquaredMelting = 0;
				double sumSquaredInterpolation = 0;

				// Print Result and calculate Deviations 
				for(int i=0;i<finalIMTime/timeStep+1;i++){
					System.out.println(valuesSwaption[0][i].getAverage() + "\t" + valuesSwaption[1][i].getAverage() + "\t" + valuesSwaption[2][i].getAverage());
					double errorMelting = (valuesSwaption[1][i].getAverage()-valuesSwaption[0][i].getAverage());
					deviationSumMelting += errorMelting;
					double errorInterpolation = (valuesSwaption[2][i].getAverage()-valuesSwaption[0][i].getAverage());
					deviationSumInterpolation += errorInterpolation;
					sumIM +=valuesSwaption[0][i].getAverage();
					sumSquaredMelting += errorMelting*errorMelting;
					sumSquaredInterpolation += errorInterpolation*errorInterpolation;
				}

				double averageDeviationInPercentOfSumMelting = deviationSumMelting/sumIM;
				double averageDeviationInPercentOfSumInterpolation = deviationSumInterpolation/sumIM;
				double RMSEMelting = Math.sqrt(sumSquaredMelting/valuesSwaption[0].length);
				double RMSEInterpolation = Math.sqrt(sumSquaredInterpolation/valuesSwaption[0].length);

				System.out.println(exerciseDates[exerciseIndex] + "\t" + numberOfPeriods[swapPeriodsIndex] + 
						"\t" + formatterTime.format((timeEnd-timeStart)/1000.0)+"s" + 
						"\t" + formatterTime.format((timeEndMelting-timeStartMelting)/1000.0)+"s"+ 
						"\t" + formatterTime.format((timeEndInterpolation-timeStartInterpolation)/1000.0)+"s" + 
						"\t" + formatterValue.format(averageDeviationInPercentOfSumMelting) + 
						"\t" + formatterValue.format(averageDeviationInPercentOfSumInterpolation) + 
						"\t" + RMSEMelting + "\t" + RMSEInterpolation);
				productIndex++;
			}
		}


	}

	public static AbstractSIMMProduct[] createProducts(TestProductType type, double[] exerciseDates,int[] periodNumber, ForwardCurve forwardCurve, DiscountCurve discountCurve) throws CalculationException{

		ArrayList<AbstractSIMMProduct> products = new ArrayList<>();
		double[]   fixingDates;
		double[]   paymentDates; 	
		double[]   swapRates;
		double[]   swapTenor;
		double[]   periodLength;
		double[]   periodNotionals;
		boolean[]   isPeriodStartDateExerciseDate;

		switch(type){

		case Swaps: 

			// 1) Swap Input
			double     startTime            = 0.0;	// Exercise date
			for(int i=0; i< periodNumber.length; i++){
				fixingDates     = new double[periodNumber[i]];
				paymentDates    = new double[periodNumber[i]];  	
				swapRates       = new double[periodNumber[i]];
				swapTenor       = new double[periodNumber[i]+1];

				// Fill data
				fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(n->startTime+n*0.5).toArray();
				paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(n->startTime+(n+1)*0.5).toArray();
				swapTenor = IntStream.range(0, periodNumber[i]+1).mapToDouble(n->startTime+n*0.5).toArray();
				Arrays.fill(swapRates, SIMMTest.getParSwaprate(forwardCurve, discountCurve, swapTenor)); 

				products.add(new SIMMSimpleSwap(fixingDates, paymentDates, swapRates, true /*isPayFix*/,notional, new String[]{"OIS", "Libor6m"}, "EUR"));
			}
			break;
		case Swaptions:

			for(int i=0; i<exerciseDates.length; i++){
				for(int j= 0; j<periodNumber.length; j++){

					fixingDates     = new double[periodNumber[j]];
					paymentDates    = new double[periodNumber[j]];  	
					swapRates       = new double[periodNumber[j]];
					swapTenor       = new double[periodNumber[j]+1];

					periodLength    = new double[paymentDates.length];
					periodNotionals = new double[periodLength.length];

					int index = i;
					// Set values
					fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(n->exerciseDates[index]+(n+1)*0.5).toArray();
					swapTenor = IntStream.range(0, periodNumber[j]+1).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					Arrays.fill(periodLength, 0.5);
					Arrays.fill(periodNotionals, notional);
					Arrays.fill(swapRates, SIMMTest.getParSwaprate(forwardCurve, discountCurve, swapTenor)); 

					products.add(new SIMMSwaption(exerciseDates[index], fixingDates, paymentDates, swapRates, notional, 
							DeliveryType.Physical, new String[]{"OIS","Libor6m"}, "EUR"));

				}
			}
			break;

		case BermudanCallable:

			for(int i=0; i<exerciseDates.length; i++){
				for(int j= 0; j<periodNumber.length; j++){

					fixingDates     = new double[periodNumber[j]];
					paymentDates    = new double[periodNumber[j]];  	
					swapRates       = new double[periodNumber[j]];
					swapTenor       = new double[periodNumber[j]+1];

					periodLength    = new double[paymentDates.length];
					periodNotionals = new double[periodLength.length];
					isPeriodStartDateExerciseDate = new boolean[periodLength.length];
					Arrays.fill(isPeriodStartDateExerciseDate, false);
					isPeriodStartDateExerciseDate[0]=true;
					isPeriodStartDateExerciseDate[4]=true;
					isPeriodStartDateExerciseDate[8]=true;
					isPeriodStartDateExerciseDate[12]=true;

					int index = i;
					// Set values
					fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(n->exerciseDates[index]+(n+1)*0.5).toArray();
					swapTenor = IntStream.range(0, periodNumber[j]+1).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					Arrays.fill(periodLength, 0.5);
					Arrays.fill(periodNotionals, notional);
					Arrays.fill(swapRates, SIMMTest.getParSwaprate(forwardCurve, discountCurve, swapTenor)); 
					products.add(new SIMMBermudanSwaption(fixingDates, periodLength, paymentDates, periodNotionals,
							swapRates, isPeriodStartDateExerciseDate, ExerciseType.Callable, new String[]{"OIS", "Libor6m"}, "EUR"));

				}
			}
			break;

		case BermudanCancelable:

			for(int i=0; i<exerciseDates.length; i++){
				for(int j= 0; j<periodNumber.length; j++){

					fixingDates     = new double[periodNumber[j]];
					paymentDates    = new double[periodNumber[j]];  	
					swapRates       = new double[periodNumber[j]];
					swapTenor       = new double[periodNumber[j]+1];

					periodLength    = new double[paymentDates.length];
					periodNotionals = new double[periodLength.length];
					isPeriodStartDateExerciseDate = new boolean[periodLength.length];
					Arrays.fill(isPeriodStartDateExerciseDate, false);
					isPeriodStartDateExerciseDate[0]=true;
					isPeriodStartDateExerciseDate[4]=true;
					isPeriodStartDateExerciseDate[8]=true;
					isPeriodStartDateExerciseDate[12]=true;

					int index = i;
					// Set values
					fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(n->exerciseDates[index]+(n+1)*0.5).toArray();
					swapTenor = IntStream.range(0, periodNumber[j]+1).mapToDouble(n->exerciseDates[index]+n*0.5).toArray();
					Arrays.fill(periodLength, 0.5);
					Arrays.fill(periodNotionals, notional);
					Arrays.fill(swapRates, SIMMTest.getParSwaprate(forwardCurve, discountCurve, swapTenor)); 
					products.add(new SIMMBermudanSwaption(fixingDates, periodLength, paymentDates, periodNotionals,
							swapRates, isPeriodStartDateExerciseDate, ExerciseType.Cancelable, new String[]{"OIS", "Libor6m"}, "EUR"));

				}
			}
			break;

		default:
			break;

		}
		return products.stream().toArray(AbstractSIMMProduct[]::new);

	}

}
