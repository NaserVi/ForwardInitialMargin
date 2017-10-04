package initialmargin.simm;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;

import initialmargin.simm.SIMMAAD.CurrencyVolatility;
import initialmargin.simm.SIMMAAD.WeightToLiborAdjustmentMethod;
import initialmargin.simm.changedfinmath.LIBORMarketModel;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.automaticdifferentiation.backward.RandomVariableDifferentiableAADFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModelFromGivenMatrix;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Swap;
import net.finmath.montecarlo.interestrate.products.SwapLeg;
import net.finmath.montecarlo.interestrate.products.components.AbstractNotional;
import net.finmath.montecarlo.interestrate.products.components.Notional;
import net.finmath.montecarlo.interestrate.products.indices.AbstractIndex;
import net.finmath.montecarlo.interestrate.products.indices.LIBORIndex;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;

public class SIMMTestAAD {
	final static DecimalFormat formatterTime	= new DecimalFormat("0.000");
	final static DecimalFormat formatterIM  	= new DecimalFormat("0.00000000000");
	
	// LIBOR Market Model parameters
	private final static int numberOfPaths		= 500;
	private final static int numberOfFactors	= 1;
		
	public static void main(String[] args) throws SolverException, CloneNotSupportedException, CalculationException {
		
		 Map<String, Object> properties = new HashMap<String, Object>();
     	 properties.put("isGradientRetainsLeafNodesOnly", new Boolean(false));
     	 AbstractRandomVariableFactory randomVariableFactory = new RandomVariableDifferentiableAADFactory(new RandomVariableFactory(), properties);
     	 
     	 // Create a Libor market Model
     	 DiscountCurve discountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve",
     			                                                                            new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
     			                                                                            new double[] {0.996 , 0.995, 0.994, 0.993, 0.98} /*discountFactors*/);
     	 ForwardCurve  forwardCurve = ForwardCurve.createForwardCurveFromForwards("forwardCurve",
 					                                                              new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,					                                                            
 					                                                              new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
 					                                                              0.5/* tenor / period length */);
 					
     	 LIBORModelMonteCarloSimulationInterface model = createLIBORMarketModel(randomVariableFactory,numberOfPaths, numberOfFactors, 
     				                                                            discountCurve,
     				                                                            forwardCurve,0.0 /* Correlation */);
     	 
     	 // Create another Libor market model with steeper curves
     	 DiscountCurve discountCurveSteep = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve2",
     			                                                                                 new double[] {0.5 , 1.0, 2.0, 5.0, 30.0} /*times*/,
                                                                                                 new double[] {0.99 , 0.98, 0.97, 0.9, 0.7} /*discountFactors*/);
         ForwardCurve  forwardCurveSteep = ForwardCurve.createForwardCurveFromForwards("forwardCurve2",
                                                                                       new double[] {0.5 , 1.0, 2.0, 5.0, 30.0}	/* fixings of the forward */,
                                                                                       new double[] {0.02, 0.02, 0.02, 0.02, 0.02},
                                                                                       0.5/* tenor / period length */);

         LIBORModelMonteCarloSimulationInterface model2 = createLIBORMarketModel(randomVariableFactory,numberOfPaths, numberOfFactors, 
                                                                                 discountCurveSteep,
                                                                                 forwardCurveSteep, 0.0 /* Correlation */);
        
     	// IM Portfolio Products. First test: Simple IR Swap
 		AbstractLIBORMonteCarloProduct[] products = new Swap[1];
 		products = createSwaps(new String[] {"5Y"});
 		double lifeTime = 5.0;
 		
        // SIMM
		SIMMAAD SIMM_StochasticWeightAdj = new SIMMAAD(model, products, CurrencyVolatility.RegularCurrency, WeightToLiborAdjustmentMethod.Stochastic, true /*MulticurveFlag*/);
		SIMMAAD SIMM_ConstantWeightAdj   = new SIMMAAD(model, products, CurrencyVolatility.RegularCurrency, WeightToLiborAdjustmentMethod.Constant, true /*MulticurveFlag*/);
		
		SIMMAAD SIMM_StochasticWeightAdjSteep = new SIMMAAD(model2, products, CurrencyVolatility.RegularCurrency, WeightToLiborAdjustmentMethod.Stochastic, true /*MulticurveFlag*/);
		SIMMAAD SIMM_ConstantWeightAdjSteep   = new SIMMAAD(model2, products, CurrencyVolatility.RegularCurrency, WeightToLiborAdjustmentMethod.Constant, true /*MulticurveFlag*/);

		
		// --------------------------------------------------------------------------------------------------------------------------
		// 1) Compare swap rate forward sensis dV/dS using stochastic dL/dS and constant dL/dS
		// --------------------------------------------------------------------------------------------------------------------------
		
		double initialMarginTime = 1.5;
		
		// Get dVdL (same for both weight adjustment methods). dVdL is always stochastic.
		RandomVariableInterface[]   dVdL = SIMM_StochasticWeightAdj.getValueLiborSensitivities(initialMarginTime);
		
		// Get weakly stochastic dLdS(t) 
		RandomVariableInterface[][] dLdS_StochasticFlat  = SIMM_StochasticWeightAdj.getLiborSwapSensitivities(initialMarginTime);
		RandomVariableInterface[][] dLdS_StochasticSteep = SIMM_StochasticWeightAdjSteep.getLiborSwapSensitivities(initialMarginTime);
		
		// Get constant dLdS(t=0) 
		RandomVariableInterface[][] dLdS_ConstantFlat = SIMM_ConstantWeightAdj.getLiborSwapSensitivities(initialMarginTime);
		RandomVariableInterface[][] dLdS_ConstantSteep = SIMM_ConstantWeightAdjSteep.getLiborSwapSensitivities(initialMarginTime);
		
		System.out.println("Comparison of swap rate forward sensis at time " + formatterTime.format(initialMarginTime) + "for flat and steep discount curve");
		System.out.println("Maturity S     " + "\t" + "dVdS stochastic flat    " + "\t" +  "dVdS constant flat" + 
		                                       "\t" + "dVdS stochastic steep   " + "\t" +  "dVdS constant steep");
		// Calculate dVdS
		for(int swapIndex=0;swapIndex<(int)(lifeTime/0.5);swapIndex++){
			
		    RandomVariableInterface dVdS_StochasticFlat  = new RandomVariable(0.0);
			RandomVariableInterface dVdS_ConstantFlat    = new RandomVariable(0.0);
			RandomVariableInterface dVdS_StochasticSteep = new RandomVariable(0.0);
			RandomVariableInterface dVdS_ConstantSteep   = new RandomVariable(0.0);
			RandomVariableInterface factor;
			
			for(int liborIndex=0;liborIndex<dLdS_StochasticFlat.length-1;liborIndex++){
			    
				factor = dLdS_StochasticFlat[liborIndex][swapIndex]==null ? new RandomVariable(0.0) : dLdS_StochasticFlat[liborIndex][swapIndex];
			    dVdS_StochasticFlat=dVdS_StochasticFlat.addProduct(dVdL[liborIndex], factor);
			    
			    factor = dLdS_ConstantFlat[liborIndex][swapIndex]==null? new RandomVariable(0.0) : dLdS_ConstantFlat[liborIndex][swapIndex];
				dVdS_ConstantFlat = dVdS_ConstantFlat.addProduct(dVdL[liborIndex], factor);
				
				factor = dLdS_StochasticSteep[liborIndex][swapIndex]==null ? new RandomVariable(0.0) : dLdS_StochasticSteep[liborIndex][swapIndex];
			    dVdS_StochasticSteep=dVdS_StochasticSteep.addProduct(dVdL[liborIndex], factor);
			    
			    factor = dLdS_ConstantSteep[liborIndex][swapIndex]==null? new RandomVariable(0.0) : dLdS_ConstantSteep[liborIndex][swapIndex];
				dVdS_ConstantSteep = dVdS_ConstantSteep.addProduct(dVdL[liborIndex], factor);
			}
			System.out.println(formatterTime.format((swapIndex+1)*0.5) + "      \t" + 
			                   formatterIM.format(dVdS_StochasticFlat.getAverage())+ "   \t" +
			                   formatterIM.format(dVdS_ConstantFlat.getAverage()) + "    \t" +
			                   formatterIM.format(dVdS_StochasticSteep.getAverage())+ "  \t" +
			                   formatterIM.format(dVdS_ConstantSteep.getAverage()));
		
		}
		System.out.println("---------------------------------------------------------------------------------------------------------");
		
		// --------------------------------------------------------------------------------------------------------------------------
		// 2) Calculate forward Initial Margin over time 
		// --------------------------------------------------------------------------------------------------------------------------
		System.out.println("Forward Initial Margin over time");
		
		// Choose if discount curve should be ignored in SIMM (incorrect if ignored; just to see how long it takes with and without ignoring it)
		SIMM_StochasticWeightAdj.setIgnoreDiscountCurve(true); // Takes much longer if false. Different implementation / method required ?!
		SIMM_ConstantWeightAdj.setIgnoreDiscountCurve(true);
		
		double finalTime = 5.0; // The last time of the IM exposure to be calculated 
		double timeStep  = 0.1;
		
		// Perform Calculations
		// 1) Stochastic
		long timeStochasticStart = System.currentTimeMillis();
		RandomVariableInterface[] forwardIMStochastic = SIMM_StochasticWeightAdj.getForwardIM(finalTime, timeStep);
		long timeStochasticEnd   = System.currentTimeMillis();
		
		// 2) Constant
		long timeConstantStart = System.currentTimeMillis();
		RandomVariableInterface[] forwardIMConstant = SIMM_ConstantWeightAdj.getForwardIM(finalTime, timeStep);
		long timeConstantEnd   = System.currentTimeMillis();
		
		// Printing results
		System.out.println("Calculation Time Stochastic = " + formatterTime.format((timeStochasticEnd-timeStochasticStart)/1000.0) + "s");
		System.out.println("Calculation Time Constant   = " + formatterTime.format((timeConstantEnd-timeConstantStart)/1000.0) + "s");
		System.out.println("Time" + "\t" + "IM Stochastic Weights" + "\t" +   "IM Constant Weights");	
		for(int i = 0; i<(finalTime/timeStep);i++){
		   double forwardIMTime = i*timeStep;
		   System.out.println(formatterTime.format(forwardIMTime) + "\t" +
		                      formatterIM.format(forwardIMStochastic[i].getAverage())  + "      \t" +
		                      formatterIM.format(forwardIMConstant[i].getAverage()));
		}
		
		// Time Grid Adjustment by band matrix yields strange result. Check what happens without the Band Matrix Adjustment
		SIMM_ConstantWeightAdj.setUseTimeGridAdjustment(false);
		forwardIMConstant = SIMM_ConstantWeightAdj.getForwardIM(finalTime, timeStep);
		// Printing results
		System.out.println("Time" + "\t" + "IM WITHOUT Time Grid Adjustment by Band Matrix dL/dL (Constant Weight)");	
		for(int i = 0; i<(finalTime/timeStep);i++){
		   double forwardIMTime = i*timeStep;
		   System.out.println(formatterTime.format(forwardIMTime) + "\t" +
		                      /*formatterIM.format(*/forwardIMConstant[i].getAverage());
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
	
}

