package initialmargin.isdasimm.test;

import java.util.Arrays;
import java.util.stream.IntStream;

import net.finmath.exception.CalculationException;
import initialmargin.isdasimm.aggregationscheme.CalculationSchemeInitialMarginISDA;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.products.SIMMSimpleSwap;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;

public class RiskWeightCalibrationTest {

	public final static int numberOfPaths = 100;
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
		
		
		/*
		 * Create calibration products
		 */
		SIMMSimpleSwap[] calibrationProducts = createCalibrationProducts(1000000 /*notional*/, new int[]{10, 20, 30} /*numberOfPeriods*/, 0.5 /*periodLength*/, forwardCurve, discountCurve);
        double[] targetIMValues = new double[]{25000,	38000,	55000};

        CalculationSchemeInitialMarginISDA isdaScheme = new CalculationSchemeInitialMarginISDA("EUR");
        System.out.println(calibrationProducts[0].getInitialMargin(0.0, model, "EUR").getAverage() + "\t" + calibrationProducts[1].getInitialMargin(0.0, model, "EUR").getAverage() + "\t" + calibrationProducts[2].getInitialMargin(0.0, model, "EUR").getAverage());
        double[] calibratedRiskWeights = isdaScheme.getRiskWeightsCalibrated(model, calibrationProducts, targetIMValues, null /*parameters*/);
	
        for(int i = 0; i<calibratedRiskWeights.length; i++) System.out.println(calibratedRiskWeights[i]);

		
	}

	public static SIMMSimpleSwap[] createCalibrationProducts(double notional, int[] numberOfPeriods, double periodLength, ForwardCurve forwardCurve, DiscountCurve discountCurve) throws CalculationException{

		SIMMSimpleSwap[] swaps = new SIMMSimpleSwap[numberOfPeriods.length];

		for(int swapIndex=0; swapIndex < numberOfPeriods.length; swapIndex++){
			double[]   fixingDatesSwap     = new double[numberOfPeriods[swapIndex]];
			double[]   paymentDatesSwap    = new double[numberOfPeriods[swapIndex]];  	
			double[]   swapRatesSwap       = new double[numberOfPeriods[swapIndex]];
			double[]   swapTenorSwap       = new double[numberOfPeriods[swapIndex]+1];

			fixingDatesSwap = IntStream.range(0, fixingDatesSwap.length).mapToDouble(i->i*0.5).toArray();
			paymentDatesSwap = IntStream.range(0, paymentDatesSwap.length).mapToDouble(i->(i+1)*periodLength).toArray();
			swapTenorSwap = IntStream.range(0, numberOfPeriods[swapIndex]+1).mapToDouble(i->i*periodLength).toArray();
			Arrays.fill(swapRatesSwap, LMMCalibrationTest.getParSwaprate(forwardCurve, discountCurve, swapTenorSwap));

			swaps[swapIndex] = new SIMMSimpleSwap(fixingDatesSwap, paymentDatesSwap, swapRatesSwap, true /*isPayFix*/,notional, new String[]{"OIS", "Libor6m"}, "EUR");
		}
		return swaps;

	}

}
