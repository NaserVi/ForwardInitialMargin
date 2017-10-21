package initialmargin.isdasimm;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import initialmargin.isdasimm.SIMMPortfolio.SensitivityMode;
import initialmargin.isdasimm.SIMMPortfolio.WeightToLiborAdjustmentMethod;
import initialmargin.regression.changedfinmath.products.Portfolio;
import initialmargin.simm.changedfinmath.products.Swaption;
import initialmargin.simm.SIMMAAD;
import initialmargin.simm.SIMMTestAAD;
import initialmargin.simm.SIMMAAD.CurrencyVolatility;
import initialmargin.simm.changedfinmath.LIBORMarketModel;
import initialmargin.simm.changedfinmath.LIBORMarketModelInterface;
import initialmargin.simm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.simm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.simm.changedfinmath.products.*;
import initialmargin.simm.changedfinmath.products.components.AbstractNotional;
import initialmargin.simm.changedfinmath.products.components.Notional;
import initialmargin.simm.changedfinmath.products.indices.AbstractIndex;
import initialmargin.simm.changedfinmath.products.indices.LIBORIndex;
import net.finmath.analytic.model.curves.DiscountCurve;
import net.finmath.analytic.model.curves.DiscountCurveInterface;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAAD;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;


public class SIMMTest {
   final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
   
   public static void main(String[] args) throws CalculationException{
	   
	   // Create a Libor market Model
	   AbstractRandomVariableFactory randomVariableFactory = createRandomVariableFactoryAAD();
   	   DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
   			                                                                            new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
   			                                                                            getRVAAD(new double[] {0.996 , 0.995, 0.994, 0.993, 0.98}) /*discountFactors*/);
   			                                                                            
   	   ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
					                                                              new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,					                                                            
					                                                              new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
					                                                              0.5/* tenor / period length */);
					
   	   LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(randomVariableFactory,1000/*numberOfPaths*/, 1 /*numberOfFactors*/, 
   				                                                            discountCurve,
   				                                                            forwardCurve,0.0 /* Correlation */);
   	   
  
   	   // (Bermudan) Swaption
  	   double     exerciseDate  = 5.0;	// Exercise date
  	   double[]   fixingDates   = {5.0, 5.5, 6.0, 6.5, 7.0,7.5,8.0,8.5};   // Vector of fixing dates (must be sorted)
  	   double[]   paymentDates  = {    5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0};	// Vector of payment dates (same length as fixing dates)
  	   double[]   periodLength  = new double[paymentDates.length]; Arrays.fill(periodLength, 0.5);
  	   double[]   periodNotionals = new double[periodLength.length]; Arrays.fill(periodNotionals, 100.0);
  	   boolean[]  isPeriodStartDateExerciseDate = new boolean[periodLength.length]; Arrays.fill(isPeriodStartDateExerciseDate, false);
  	   isPeriodStartDateExerciseDate[0]=true;
  	   
  	   double[]   swapRates     = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};// Vector of strikes
  		
  	   AbstractLIBORMonteCarloProduct swaption = new Swaption(exerciseDate,fixingDates,paymentDates,swapRates,100.0,"Physical");
       AbstractLIBORMonteCarloProduct bermudan = new BermudanSwaption(isPeriodStartDateExerciseDate,fixingDates,periodLength,paymentDates,periodNotionals, swapRates);
   	   
       // Create ClassifiedSIMMProduct Swap
	   AbstractLIBORMonteCarloProduct swap = SIMMTestAAD.createSwaps(new String[] {"5Y"})[0];
	   AbstractLIBORMonteCarloProduct swap2 = SIMMTestAAD.createSwaps(new String[] {"3Y"})[0];
	   
	   // Classify the products 
	   SIMMClassifiedProduct product1 = new SIMMClassifiedProduct(swaption,"RatesFX",new String[] {"InterestRate"}, new String[] {"OIS","Libor6m"},"EUR",null,false,false);
	   SIMMClassifiedProduct product2 = new SIMMClassifiedProduct(swap2,"RatesFX",new String[] {"InterestRate"}, new String[] {"OIS","Libor6m"},"EUR",null,false, false);
	   
	   SIMMPortfolio portfolioST = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
			                                         SensitivityMode.Stochastic,
			                                         WeightToLiborAdjustmentMethod.Constant);
	   SIMMPortfolio portfolioLM = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
                                                     SensitivityMode.LinearMelting,
                                                     WeightToLiborAdjustmentMethod.Constant);

	   double finalIMTime=9.0;
	   
//	   System.out.println("Survival Probabilities");
//	   for(int i=0;i<finalIMTime/0.125;i++) System.out.println(portfolioLM.getProducts()[0].getSurvivalProbability(i*0.125,model,true));
	   
	   // 1)Calculate forward sensis by AAD at each time point
	   double[] valuesST = new double[(int)(finalIMTime/0.125)];
	   
	   long timeSTStart = System.currentTimeMillis();
	     for(int i=0;i<finalIMTime/0.125;i++) valuesST[i] = portfolioST.getValue(i*0.125, model).getAverage();
	   long timeSTEnd = System.currentTimeMillis();
//	   System.out.println("Time with calculation of AAD sensis at each time point: " + formatterTime.format((timeSTEnd-timeSTStart)/1000.0)+"s");
//
//	   // 1) Melt sensis linearly 
//  	   double[] valuesLM = new double[(int)(finalIMTime/0.125)];
//  	   
//  	   long timeLMStart = System.currentTimeMillis();
//	     for(int i=0;i<finalIMTime/0.125;i++) valuesLM[i] = portfolioLM.getValue(i*0.125, model).getAverage();
//	   long timeLMEnd = System.currentTimeMillis();
//  	
//	   System.out.println("Time with Melting: " + formatterTime.format((timeLMEnd-timeLMStart)/1000.0)+"s");
//	          
	   System.out.println("IM Linear Melting" + "\t" + "IM Forward AAD Sensis");
       
	   for(int i=0;i<finalIMTime/0.125;i++){
    	   System.out.println(/*valuesLM[i] + "\t" + */valuesST[i]);
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
			appliedDiscountCurve = (DiscountCurveInterface) new DiscountCurveFromForwardCurve(forwardCurve);
		} else {
			appliedDiscountCurve = discountCurve;
		}
		/*
		 * Create a simulation time discretization
		 */
		double lastTime	= 10.0;
		double dt		= 0.125;

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
	
	public static RandomVariableInterface[] getRVAAD(double[] rates){
		RandomVariableInterface[] rv = new RandomVariableInterface[rates.length];
		for(int i=0;i<rv.length;i++) rv[i]=new RandomVariableDifferentiableAAD(rates[i]);
		return rv;
	}

	

               
         
	   
   }

