package initialmargin.isdasimm;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

import initialmargin.isdasimm.SIMMPortfolio.SensitivityMode;
import initialmargin.isdasimm.SIMMPortfolio.WeightToLiborAdjustmentMethod;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModel;
import initialmargin.isdasimm.changedfinmath.LIBORMarketModelInterface;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulation;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;
import initialmargin.isdasimm.changedfinmath.products.BermudanSwaption;
import initialmargin.isdasimm.changedfinmath.products.SimpleSwap;
import initialmargin.isdasimm.changedfinmath.products.Swap;
import initialmargin.isdasimm.changedfinmath.products.SwapLeg;
import initialmargin.isdasimm.changedfinmath.products.Swaption;
import initialmargin.isdasimm.changedfinmath.products.components.AbstractNotional;
import initialmargin.isdasimm.changedfinmath.products.components.Notional;
import initialmargin.isdasimm.changedfinmath.products.indices.AbstractIndex;
import initialmargin.isdasimm.changedfinmath.products.indices.LIBORIndex;
import initialmargin.isdasimm.old.SIMMTestAADold;
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
					
   	   LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(randomVariableFactory,500/*numberOfPaths*/, 1 /*numberOfFactors*/, 
   				                                                            discountCurve,
   				                                                            forwardCurve,0.0 /* Correlation */);
   	   
  
   	   
   	   // Create Products
  	   double     exerciseTime  = 6.0;	// Exercise date
  	   double     constantSwapRate = 0.01;
  	   int        numberOfPeriods = 6;
  	   double     notional        = 100;
  	   
  	   double[]   fixingDates     = new double[numberOfPeriods];
  	   double[]   paymentDates    = new double[numberOfPeriods];
  	   double[]   periodLength    = new double[paymentDates.length];
  	   double[]   periodNotionals = new double[periodLength.length];
  	   double[]   swapRates       = new double[numberOfPeriods];
  	   // for Bermudan
  	   boolean[]  isPeriodStartDateExerciseDate = new boolean[periodLength.length];
  	   // Set values
  	   fixingDates = IntStream.range(0, fixingDates.length).mapToDouble(i->exerciseTime+i*0.5).toArray();
  	   paymentDates = IntStream.range(0, paymentDates.length).mapToDouble(i->exerciseTime+(i+1)*0.5).toArray();
  	   Arrays.fill(periodLength, 0.5);
  	   Arrays.fill(periodNotionals, notional);
  	   Arrays.fill(swapRates, constantSwapRate); 
  	   Arrays.fill(isPeriodStartDateExerciseDate, false);
//  	   isPeriodStartDateExerciseDate[0]=true;
//  	   isPeriodStartDateExerciseDate[2]=true;
//  	   isPeriodStartDateExerciseDate[4]=true;
//  	   isPeriodStartDateExerciseDate[6]=true;
  	  
  	   // Create Products
  	   AbstractLIBORMonteCarloProduct simpleSwap = new SimpleSwap(fixingDates,paymentDates,swapRates,100);
  	   AbstractLIBORMonteCarloProduct swaption = new Swaption(exerciseTime,fixingDates,paymentDates,swapRates,100.0,"Physical");
       AbstractLIBORMonteCarloProduct bermudan = new BermudanSwaption(isPeriodStartDateExerciseDate,fixingDates,periodLength,paymentDates,periodNotionals, swapRates);
	   AbstractLIBORMonteCarloProduct swap = SIMMTestAADold.createSwaps(new String[]  {"5Y"})[0];
	   AbstractLIBORMonteCarloProduct swap2 = SIMMTestAADold.createSwaps(new String[] {"3Y"})[0];
	   
	   // Classify the products 
	   SIMMClassifiedProduct product1 = new SIMMClassifiedProduct(simpleSwap,"RatesFX",new String[] {"InterestRate"}, new String[] {"OIS","Libor6m"},"EUR","null",true,false);
	   SIMMClassifiedProduct product2 = new SIMMClassifiedProduct(swap2,"RatesFX",new String[] {"InterestRate"}, new String[] {"OIS","Libor6m"},"EUR",null,false, false);
	   
	   // Create SIMMPortfolios
	   double sensiResetStep = 2.0; // Time step at which Sensis are reset for melting: if > exercise Date we have no reset
	   SIMMPortfolio portfolioST = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
			                                         SensitivityMode.Stochastic,
			                                         WeightToLiborAdjustmentMethod.Constant, 0.0);
	   SIMMPortfolio portfolioIP = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
                                                     SensitivityMode.Interpolation,
                                                     WeightToLiborAdjustmentMethod.Constant, sensiResetStep);
	   SIMMPortfolio portfolioLB = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
			   										 SensitivityMode.LinearMelting,
			   										 WeightToLiborAdjustmentMethod.Constant, 50 /*sensiResetStep ( -> no reset)*/);
	   SIMMPortfolio portfolioLL = new SIMMPortfolio(new SIMMClassifiedProduct[] {product1},"EUR",
                                                     SensitivityMode.LinearMeltingOnLiborBuckets,
                                                     WeightToLiborAdjustmentMethod.Constant, sensiResetStep);

	   // Perform calculations
	   double finalIMTime=exerciseTime+0.5*numberOfPeriods;
	   double timeStep = 0.1;
	
	   //System.out.println("Survival Probabilities");
	   //for(int i=0;i<finalIMTime/0.125;i++) System.out.println(portfolioLM.getProducts()[0].getSurvivalProbability(i*0.125,model,true));
	   
	   // 1) Melt sensis linearly on buckets
  	   double[] valuesLB = new double[(int)(finalIMTime/timeStep)];
  	  
  	   long timeLBStart = System.currentTimeMillis();
	     for(int i=0;i<finalIMTime/timeStep;i++) valuesLB[i] = portfolioST.getValue(4.0, model).getAverage();
	   long timeLBEnd = System.currentTimeMillis();
  	
	   System.out.println("Time with Melting on Buckets: " + formatterTime.format((timeLBEnd-timeLBStart)/1000.0)+"s");
	      
	  // 2) Melt sensis linearly on LiborPeriodDiscretization
//     double[] valuesLL = new double[(int)(finalIMTime/timeStep)];
//  	  
//     long timeLLStart = System.currentTimeMillis();
//	     for(int i=0;i<finalIMTime/timeStep;i++) valuesLL[i] = portfolioLL.getValue(i*timeStep, model).getAverage();
//	   long timeLLEnd = System.currentTimeMillis();
  	
//	   System.out.println("Time with Melting on LiborPeriodDiscretization: " + formatterTime.format((timeLLEnd-timeLLStart)/1000.0)+"s");
	   
	   
	   // 3) Interpolate Sensitivities
	   double[] valuesIP = new double[(int)(finalIMTime/timeStep)];
	  
	   long timeIPStart = System.currentTimeMillis();
          for(int i=0;i<finalIMTime/timeStep;i++) valuesIP[i] = portfolioIP.getValue(i*timeStep, model).getAverage();
       long timeIPEnd = System.currentTimeMillis();
	
       System.out.println("Time with Interpolation: " + formatterTime.format((timeIPEnd-timeIPStart)/1000.0)+"s");
	   
       
	   // 4)Calculate forward sensis by AAD at each time point
	   double[] valuesST = new double[(int)(finalIMTime/timeStep)];
	   
	   long timeSTStart = System.currentTimeMillis();
	     for(int i=0;i<finalIMTime/timeStep;i++) valuesST[i] = portfolioST.getValue(i*timeStep, model).getAverage();
	   long timeSTEnd = System.currentTimeMillis();
	   
	   System.out.println("Time with calculation of AAD sensis at each time point: " + formatterTime.format((timeSTEnd-timeSTStart)/1000.0)+"s");


	   System.out.println("IM Linear Melting on Buckets" + "\t" + "IM Interpolation" + "\t" + "IM Forward AAD Sensis");
       
	   for(int i=0;i<finalIMTime/timeStep;i++){
    	   System.out.println(valuesLB[i] + "\t" + valuesIP[i] + "\t" + valuesST[i]);
       }
	   
	   
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

