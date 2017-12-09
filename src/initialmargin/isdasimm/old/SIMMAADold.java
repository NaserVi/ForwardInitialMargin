package initialmargin.isdasimm.old;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import net.finmath.analytic.model.curves.DiscountCurve;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import initialmargin.isdasimm.changedfinmath.*;
import initialmargin.isdasimm.changedfinmath.products.*;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;

public class SIMMAADold {
	
	private final LIBORModelMonteCarloSimulationInterface model;
	private final AbstractLIBORMonteCarloProduct[] portfolioProducts;
	private AbstractLIBORMonteCarloProduct currentProduct; // The product whose IM is being calculated 
	private Map<Long, RandomVariableInterface> gradientOfProduct; // Same for all evaluationTimes; Is reset for different products
	private final boolean isUseMultiCurve;
	private ConditionalExpectationEstimatorInterface conditionalExpectationOperator;
	private boolean isUseTimeGridAdjustment=true; // can be discarded later.. just to check how time grid adjustment affects IM
	private boolean isIgnoreDiscountCurve=false;  // can be discarded later.. just to check how discount curve influences IM
	private double[] discountCurvePillars = {0.5 , 1.0, 2.0, 5.0, 30.0};// shouldn't the discount curve know its pillars ?
	// SIMM parameters
	private final double deltaThreshold = 250000000;
	//risk factor in days
	private final int[] riskFactorsSIMM = new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950};
	private WeightToLiborAdjustmentMethod liborWeightMethod;
	
	public enum CurrencyVolatility{
		RegularCurrency,//Regular currency (USD, EUR, GBP, CHF, AUD, NZD, CAD, SEK, NOK, DKK, HKD, KRW, SGD, TWD)
		LowVolCurrency, //Low volatility currency (JPY)
		HighVolCurrency //High volatility currency (all other currencies)
	}
	
	public enum WeightToLiborAdjustmentMethod{
		Constant,  //Sets dL/dS(t=0) for all forward IM times, i.e. leave the weight adjustment dL/dS constant
		Stochastic //Calculate dL/dS(t) for all forward IM times, i.e. (weakly) stochastic weight adjustment 
	}
	
	private RandomVariableInterface[][] riskWeightToLiborAdjustments; // dL/dS(t=0). To be used if WeightToLiborAdjustmentMethod == Constant
	private double[] riskWeightsSIMM; // will be one of the following, depending on CurrencyVolatility
	private final double[] riskWeightsRegularCurrency = new double[]{0.0077, 0.0077, 0.0077, 0.0064, 0.0058, 0.0049, 0.0047, 0.0047,	0.0045,	0.0045,	0.0048,	0.0056};
	private final double[] riskWeightsLowVolCurrency  = new double[]{0.0010, 0.0010, 0.0010, 0.0010, 0.0013, 0.0016, 0.0018, 0.0020, 0.0025, 0.0022, 0.0022, 0.0023};
	private final double[] riskWeightsHighVolCurrency = new double[]{0.0089, 0.0089, 0.0089, 0.0094, 0.0104, 0.0099, 0.0096, 0.0099, 0.0087, 0.0097, 0.0097, 0.0098};
	
	private final double[][] correlationMatrixWithinSubCurve = new double[][]{
			{0    , 1    , 1    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
			{1    , 0    , 1    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
			{1    , 1    , 0    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
			{0.782, 0.782, 0.782, 0    , 0.84 , 0.739, 0.667, 0.569, 0.444, 0.375, 0.349, 0.296},
			{0.618, 0.618, 0.618, 0.84 , 0    , 0.917, 0.859, 0.757, 0.626, 0.555, 0.526, 0.471},
			{0.498, 0.498, 0.498, 0.739, 0.917, 0    , 0.976, 0.895, 0.749, 0.69 , 0.66 , 0.602},
			{0.438, 0.438, 0.438, 0.667, 0.859, 0.976, 0    , 0.958, 0.831, 0.779, 0.746, 0.69 },
			{0.361, 0.361, 0.361, 0.569, 0.757, 0.895, 0.958, 0    , 0.925, 0.893, 0.859, 0.812},
			{0.27 , 0.27 , 0.27 , 0.444, 0.626, 0.749, 0.831, 0.925, 0    , 0.98 , 0.961, 0.931},
			{0.196, 0.196, 0.196, 0.375, 0.555, 0.69 , 0.779, 0.893, 0.98 , 0    , 0.989, 0.97 },
			{0.174, 0.174, 0.174, 0.349, 0.526, 0.66 , 0.746, 0.859, 0.961, 0.989, 0    , 0.988},
			{0.129, 0.129, 0.129, 0.296, 0.471, 0.602, 0.69 , 0.812, 0.931, 0.97 , 0.988, 0    }
	};
	
	private final double[][] correlationMatrixBetweenSubCurves = new double[][]{

	{0.98,		0.98,		0.98,		0.76636,	0.60564,	0.48804,	0.42924,	0.35378,	0.2646,	0.19208,0.17052,0.12642},
	{0.98,		0.98,		0.98,		0.76636,	0.60564, 	0.48804,	0.42924,	0.35378,	0.2646,	0.19208,0.17052,0.12642},
	{0.98,		0.98,		0.98,		0.76636,	0.60564,	0.48804,	0.42924,	0.35378,	0.2646,	0.19208,0.17052,0.12642},
	{0.76636,	0.76636,	0.76636,	0.98,		0.8232,		0.72422,	0.65366,	0.55762,	0.43512,0.3675,	0.34202,0.29008},
	{0.60564,	0.60564, 	0.60564,	0.8232,		0.98,		0.89866,	0.84182, 	0.74186,	0.61348,0.5439,	0.51548,0.46158},
	{0.48804,	0.48804,	0.48804,	0.72422,	0.89866,	0.98,		0.95648,	0.8771,		0.73402,0.6762,	0.6468,	0.58996},
	{0.42924,	0.42924,	0.42924,	0.65366,	0.84182,	0.95648,	0.98,		0.93884,	0.81438,0.76342,0.73108,0.6762},
	{0.35378,	0.35378,	0.35378,	0.55762,	0.74186,	0.8771,		0.93884,	0.98,		0.9065,	0.87514,0.84182,0.79576},
	{0.2646,	0.2646,		0.2646,		0.43512,	0.61348,	0.73402,	0.81438,	0.9065,		0.98,	0.9604,	0.94178,0.91238},
	{0.19208,	0.19208,	0.19208,	0.3675,		0.5439,		0.6762,		0.76342,	0.87514,	0.9604,	0.98,	0.96922,0.9506},
	{0.17052,	0.17052,	0.17052,	0.34202,	0.51548,	0.6468,		0.73108,	0.84182,	0.94178,0.96922,0.98,	0.96824},
	{0.12642,	0.12642,	0.12642,	0.29008,	0.46158,	0.58996,	0.6762,		0.79576,	0.91238,0.9506,	0.96824,0.98}
	};

	
	/**Construct object of stochastic SIMM
	 * 
	 * @param model the interface implemented by Libor Market Model
	 * @param portfolioProducts array of products of the portfolio
	 * @param currencyVol type of currency
	 * @param liborWeightMethod The method for transforming weights (stochastic or constant)
	 * @param isUseMultiCurve True for multi curve model
	 * @throws CalculationException 
	 */
	public SIMMAADold(LIBORModelMonteCarloSimulationInterface model,
				   AbstractLIBORMonteCarloProduct[] portfolioProducts,
				   CurrencyVolatility currencyVol,
				   WeightToLiborAdjustmentMethod liborWeightMethod,
				   boolean isUseMultiCurve) throws CalculationException{ // should be already known by the model
		
		this.model= model;
		this.portfolioProducts = portfolioProducts;
		this.currentProduct = portfolioProducts[0];
		//double[] portfolioWeights = new double[portfolioProducts.length];
		//for(int i=0;i<portfolioWeights.length;i++) portfolioWeights[i]=1;
		//this.portfolio = new Portfolio(portfolioProducts, portfolioWeights);
		
		this.isUseMultiCurve = isUseMultiCurve;
		
		switch (currencyVol) {
		
		  case RegularCurrency:
			  this.riskWeightsSIMM = riskWeightsRegularCurrency;
			  break;
		
		  case HighVolCurrency:
			  this.riskWeightsSIMM = riskWeightsHighVolCurrency;
			  break;
		
		  case LowVolCurrency:
			  this.riskWeightsSIMM = riskWeightsLowVolCurrency;
			  break;
		}
		
		this.liborWeightMethod = liborWeightMethod;
		// Get Weight Adjustment for Libor Sensitivities if they remain constant 
		if(liborWeightMethod == WeightToLiborAdjustmentMethod.Constant){
			this.riskWeightToLiborAdjustments = getLiborSwapSensitivities(0.0 /*evaluationTime*/);
		} 
		
	}
	

	
	/**
	 * 
	 * @param finalTime The end of the initial margin interval
	 * @param timeStep  The step size for between initial margin times
	 * @return Returns a vector of forward initial margin
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getForwardIM(double finalTime, double timeStep) throws SolverException, CloneNotSupportedException, CalculationException{
		RandomVariableInterface[] initialMargin = new RandomVariableInterface[(int)(finalTime/timeStep)+1];
		Arrays.fill(initialMargin, new RandomVariable(0.0));
		// Loop over products
		for(int productIndex = 0; productIndex < getPortfolioProducts().length; productIndex++){
			setProduct(productIndex);
			setGradient(getPortfolioProduct(productIndex)); // gradient is not time dependent 
			// Loop over times
		    for(int timeIndex=0;timeIndex<initialMargin.length;timeIndex++){
		    	setConditionalExpectationOperator(timeIndex*timeStep);
			    initialMargin[timeIndex]=initialMargin[timeIndex].add(getInitialMargin(timeIndex*timeStep));
		    }
		}
		return initialMargin;
	}
	
	
//	private RandomVariableInterface getInitialMarginPKL(double evaluationTime)
//	{
//		SIMMSchemeMain.ParameterCollection dummyParameterSet = new SIMMSchemeMain.ParameterCollection();
//		SIMMSchemeMain isdaSimmScheme = new SIMMSchemeMain(this.portfolioProducts,this.model,dummyParameterSet,"EUR");
//		return isdaSimmScheme.getValue(evaluationTime);
//	}
	
	/**Calculate forward initial margin at specific time point
	 * 
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return Returns the forward initial margin
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException, CalculationException, CloneNotSupportedException
	 */
	//calculate initial margin at evaluation time
	private RandomVariableInterface getInitialMargin(double evaluationTime) throws SolverException, CalculationException, CloneNotSupportedException{
		RandomVariableInterface initialMargin = new RandomVariable(0.0);
	    // For the case of different sub yield curves for forward and discount rates (multi curve model)	
	    if(isUseMultiCurve){
	    	RandomVariableInterface[] forwardCurveSensitivity = getForwardCurveSensitivities(evaluationTime);
		    RandomVariableInterface[] discountCurveSensitivity = new RandomVariableInterface[forwardCurveSensitivity.length];
		    if(!isIgnoreDiscountCurve) discountCurveSensitivity = getDiscountCurveSensitivities(evaluationTime);
		    RandomVariableInterface concentrationThreshold = new RandomVariable(0.0); 
		
		    for(int riskFactor = 0; riskFactor < forwardCurveSensitivity.length; riskFactor++){
		    	if(isIgnoreDiscountCurve) discountCurveSensitivity[riskFactor] = new RandomVariable(0.0); 
			    concentrationThreshold = concentrationThreshold.add(forwardCurveSensitivity[riskFactor]).add(discountCurveSensitivity[riskFactor]); 
		    }
		    concentrationThreshold = concentrationThreshold.abs().div(deltaThreshold).sqrt();
		
		    for(int i = 0; i<forwardCurveSensitivity.length; i++){
			    discountCurveSensitivity[i] = discountCurveSensitivity[i].mult(riskWeightsSIMM[i]).mult(concentrationThreshold.floor(1.0));
			    forwardCurveSensitivity[i]	= forwardCurveSensitivity[i].mult(riskWeightsSIMM[i]);
		    }
		
		    for(int riskFactor1 = 0; riskFactor1 < forwardCurveSensitivity.length; riskFactor1++){
			    initialMargin = initialMargin.add(discountCurveSensitivity[riskFactor1].squared().add(forwardCurveSensitivity[riskFactor1].squared()));
			    for(int riskFactor2 = 0; riskFactor2 < discountCurveSensitivity.length; riskFactor2++){
				    RandomVariableInterface addend1 = discountCurveSensitivity[riskFactor1].mult(discountCurveSensitivity[riskFactor2]).mult(correlationMatrixWithinSubCurve[riskFactor1][riskFactor2]);
				    RandomVariableInterface addend2 = forwardCurveSensitivity[riskFactor1].mult(forwardCurveSensitivity[riskFactor2]).mult(correlationMatrixWithinSubCurve[riskFactor1][riskFactor2]);
				    RandomVariableInterface addend3 = discountCurveSensitivity[riskFactor1].mult(forwardCurveSensitivity[riskFactor2]).mult(correlationMatrixBetweenSubCurves[riskFactor1][riskFactor2]).mult(2.0);
				    initialMargin = initialMargin.add(addend1.add(addend2).add(addend3));
			    }
		    }
	    } else{ //For the case of a single curve model
		    RandomVariableInterface[] discountCurveSensitivity = getForwardCurveSensitivities(evaluationTime);//fixed sensitivity inputs: new double[]{0, 0, 0, 0, 0, 0, 0, -5000000, -10000000, -5000000, 0, 0}
		    RandomVariableInterface concentrationThreshold = new RandomVariable(0.0); 
		    for(int riskFactor = 0; riskFactor < discountCurveSensitivity.length; riskFactor++){
			    concentrationThreshold = concentrationThreshold.add(discountCurveSensitivity[riskFactor]);
		    }
		    concentrationThreshold = concentrationThreshold.abs().div(deltaThreshold).sqrt();
		    for(int i = 0; i<discountCurveSensitivity.length; i++){
			    discountCurveSensitivity[i] = concentrationThreshold.floor(1.0).mult(riskWeightsSIMM[i]).mult(discountCurveSensitivity[i]);
		    }
		    for(int riskFactor1 = 0; riskFactor1 < discountCurveSensitivity.length; riskFactor1++){	
			    initialMargin = initialMargin.add(discountCurveSensitivity[riskFactor1].squared());				   	
			    for(int riskFactor2 = 0; riskFactor2 < discountCurveSensitivity.length; riskFactor2++){	
				    initialMargin = initialMargin.add(discountCurveSensitivity[riskFactor1].mult(discountCurveSensitivity[riskFactor2]).mult(2* correlationMatrixWithinSubCurve[riskFactor1][riskFactor2]));
			    }
		    }
	    }
	return initialMargin.abs().sqrt();
	}
	
	/**Calculate forward curve sensitivities dV/dS with respect to all swap rates on all points of the LiborPeriodDiscretization.
	 * 
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return The sensitivities dV/dS for swap rates on forward curve
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	private RandomVariableInterface[] getForwardCurveSensitivities(double evaluationTime) throws SolverException, CloneNotSupportedException, CalculationException{
		// Calculate Sensitivity wrt Libors, already time grid adjusted
		RandomVariableInterface[] dVdL = getValueLiborSensitivities(evaluationTime);
		// the following line will be removed later. Just checking how timeGridAdjustment affects the result
		int timeGridIndicator = 0; if(!isUseTimeGridAdjustment && !onLiborPeriodDiscretization(evaluationTime)) timeGridIndicator = 1;
		
		RandomVariableInterface[] delta = new RandomVariableInterface[dVdL.length-timeGridIndicator];
		// Calculate time dependent risk weight adjustments for Libor rates if applicable
		RandomVariableInterface[][] dLdS;
		if(liborWeightMethod == WeightToLiborAdjustmentMethod.Stochastic){
			dLdS = getLiborSwapSensitivities(evaluationTime);
		} else {
			dLdS = riskWeightToLiborAdjustments;
		}
		// Calculate Sensitivities wrt Swaps
		for(int swapIndex = 0; swapIndex<dVdL.length-timeGridIndicator; swapIndex++){
			RandomVariableInterface dVdS  =new RandomVariable(0.0);
			RandomVariableInterface factor;
			for(int liborIndex=0;liborIndex<dVdL.length-timeGridIndicator;liborIndex++){
				      factor = dLdS[liborIndex][swapIndex]==null ?  new RandomVariable(0.0) : dLdS[liborIndex][swapIndex];
					  dVdS = dVdS.addProduct(dVdL[liborIndex+timeGridIndicator], factor);
		    }
			delta[swapIndex]=dVdS;
		}
		return getSensitivitiesOnBuckets(delta);
	}
	
	/**Performs rebucketing of sensitivities to the SIMM buckets. 
	 * 
	 * @param sensitivities The sensitivities wrt swap rates dV/dS
	 * @return The sensitivities on the SIMM buckets
	 */
	private RandomVariableInterface[] getSensitivitiesOnBuckets(RandomVariableInterface[] sensitivities){
		//rebucketing to SIMM structure(buckets: 2w, 1m, 3m, 6m, 1y, 2y, 3y, 5y, 10y, 15y, 20y, 30y)	
		int[] riskFactorDays = new int[sensitivities.length];
		RandomVariableInterface[] deltaSIMM = new RandomVariableInterface[riskFactorsSIMM.length];
		for(int i = 0;i<deltaSIMM.length;i++) deltaSIMM[i] = new RandomVariable(0.0);
		int counter = 0;
		for(int simmFactor =0; simmFactor<riskFactorsSIMM.length;simmFactor++){
			for(int i = counter; i<sensitivities.length; i++){
				
					// act/365 as default daycount convention

					riskFactorDays[i] = (int)Math.round(365 * model.getLiborPeriodDiscretization().getTime(i+1));	
					
					if(riskFactorDays[i] < riskFactorsSIMM[0]){
						deltaSIMM[0] = deltaSIMM[0].add(sensitivities[i]);
						counter++;
					}
					else{
						if(riskFactorDays[i] >= riskFactorsSIMM[riskFactorsSIMM.length-1]){
							deltaSIMM[deltaSIMM.length-1] = deltaSIMM[deltaSIMM.length-1].add(sensitivities[i]);
						}
					
						else{
							if(riskFactorDays[i] >= riskFactorsSIMM[simmFactor] && riskFactorDays[i] < riskFactorsSIMM[simmFactor+1]){
					
							deltaSIMM[simmFactor] = deltaSIMM[simmFactor].addProduct(sensitivities[i],((double)(riskFactorsSIMM[simmFactor+1] - riskFactorDays[i]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							deltaSIMM[simmFactor+1] = deltaSIMM[simmFactor+1].addProduct(sensitivities[i],((double)(riskFactorDays[i]-riskFactorsSIMM[simmFactor]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							counter++;
							}
							
							else{
							break;
							}
						}
					
					}
			}
			
		}
		
	return deltaSIMM;		
			
	}
	
	
	/**Calculates dV/dS 
	 * 
	 * @param evaluationTime The time at which the sensitivity is calculated
	 * @return The matrix dV/dS 
	 * @throws CalculationException
	 */
	public RandomVariableInterface[][] getLiborSwapSensitivities(double evaluationTime) throws CalculationException{
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		RandomVariableInterface[][] dLdS=null;
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		
		switch(getLiborWeightMethod()){
		  case Constant:{
	          evaluationTime = 0.0; 
			  int numberOfRemainingLibors = model.getNumberOfLibors();
			  dLdS = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
			
			  if(isUseMultiCurve){
			 	  // Calculate dLdS directly  
				  dLdS[0][0]=new RandomVariable(1.0);
				  double discountTime = evaluationTime+liborPeriodLength;
				  RandomVariableInterface sumDf = model.getNumeraire(discountTime).invert();
				  for(int liborIndex = 1; liborIndex<dLdS.length;liborIndex++){
					  discountTime +=model.getLiborPeriodDiscretization().getTimeStep(0);
					  RandomVariableInterface df = model.getNumeraire(discountTime).invert();
					  double denominator = df.getAverage();
				      dLdS[liborIndex][liborIndex-1]= new RandomVariable(-sumDf.getAverage()/denominator); // constant RV!
				      sumDf = sumDf.add(df);
				      dLdS[liborIndex][liborIndex] = new RandomVariable(sumDf.getAverage()/denominator);
				    } 
				} else { // Single Curve
					RandomVariableInterface[][] dLdP = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
					RandomVariableInterface numeraireAtEvalTime = model.getNumeraire(evaluationTime); //1
					double discountTime = evaluationTime + liborPeriodLength;
				    RandomVariableInterface numeraireAtDiscountTime = model.getNumeraire(discountTime);
				    RandomVariableInterface P1 = new RandomVariable(numeraireAtEvalTime.div(numeraireAtDiscountTime).getAverage());
					// Calculate dLdP
					dLdP[0][0]= P1.squared().pow(-1.0).mult(-1.0);
					RandomVariableInterface sumDf = P1;
					for(int liborIndex = 1;liborIndex<dLdP.length;liborIndex++){
						discountTime += liborPeriodLength;
						numeraireAtDiscountTime = model.getNumeraire(discountTime);
						RandomVariableInterface P2 = new RandomVariable(numeraireAtEvalTime.div(numeraireAtDiscountTime).getAverage());
					    dLdP[liborIndex][liborIndex-1] = P2.invert();
					    dLdP[liborIndex][liborIndex]   = P2.squared().invert().mult(P1).mult(-1.0);
					    sumDf=sumDf.add(P2);
					    P1 = P2;
					}
					// Get dLdS
					dLdS = multiply(dLdP,getBondSwapSensitivity(evaluationTime));
				} // end single curve
			} // end case 1 
		  break;
		  
		  case Stochastic:{
			    // Get index of first Libor starting >= evaluationTime
		        int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
				int numberOfRemainingLibors = model.getNumberOfLibors()-nextLiborIndex;
				dLdS = new RandomVariableInterface [numberOfRemainingLibors][numberOfRemainingLibors];
				
				if(isUseMultiCurve){
			 	    // Calculate dLdS directly  
				    dLdS[0][0]=new RandomVariable(1.0);
				    double discountTime = evaluationTime+liborPeriodLength;
				    RandomVariableInterface sumDf = model.getNumeraire(discountTime).invert();
				    for(int liborIndex = 1; liborIndex<dLdS.length;liborIndex++){
					    discountTime +=model.getLiborPeriodDiscretization().getTimeStep(0);
					    RandomVariableInterface df = model.getNumeraire(discountTime).invert();
					    RandomVariableInterface denominator = df.getConditionalExpectation(conditionalExpectationOperator);
				        dLdS[liborIndex][liborIndex-1]=sumDf.getConditionalExpectation(conditionalExpectationOperator).div(denominator).mult(-1.0);//dLdS[liborIndex][liborIndex-1]=-sumDf.getConditionalExpectation(conditionalExpectationOperator).getAverage()/denominator;
				        sumDf = sumDf.add(df);
				        dLdS[liborIndex][liborIndex] = sumDf.getConditionalExpectation(conditionalExpectationOperator).div(denominator);
				    }   
				} else{ // Single Curve
					RandomVariableInterface numeraireAtEvalTime = model.getNumeraire(evaluationTime);
					double discountTime = evaluationTime+liborPeriodLength;
					RandomVariableInterface numeraireAtDiscountTime = model.getNumeraire(discountTime);
					RandomVariableInterface P1 = numeraireAtEvalTime.div(numeraireAtDiscountTime).getConditionalExpectation(conditionalExpectationOperator);
					RandomVariableInterface[][] dLdP = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
					// Calculate dLdP
					dLdP[0][0]= P1.squared().pow(-1.0).mult(-1.0);
					RandomVariableInterface sumDf = P1;
					for(int liborIndex = 1;liborIndex<dLdP.length;liborIndex++){
						discountTime += liborPeriodLength;
						numeraireAtDiscountTime = model.getNumeraire(discountTime);
						RandomVariableInterface P2 = numeraireAtEvalTime.div(numeraireAtDiscountTime).getConditionalExpectation(conditionalExpectationOperator);
					    dLdP[liborIndex][liborIndex-1] = P2.invert();
					    dLdP[liborIndex][liborIndex]   = P2.squared().invert().mult(P1).mult(-1.0);
					    sumDf=sumDf.add(P2);
					    P1 = P2;
					}
					// Get dLdS
					dLdS = multiply(dLdP,getBondSwapSensitivity(evaluationTime));
				} // end single curve
		  } // end case 2
		break;
		default:
			break;
		} // end switch
		
		return dLdS;
	}
	
	
	
	/**Since dV/dL is wrt the incorrect Libor times this function provides a matrix dL/dL to be multiplied with dV/dL in order to 
	 * have the correct libor times starting at evaluationTime. 
	 * @param evaluationTime The time at which the adjustment should be calculated.
	 * @return Pseudo Inverse of derivative band matrix; Identity matrix in case of evaluationTime on LiborPeriodDiscretization; 
	 * @throws CalculationException
	 */
	private RandomVariableInterface[][] getLiborTimeGridAdjustment(double evaluationTime) throws CalculationException{
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime);
		
		// If evaluationTime lies on Libor Time Grid - return identity matrix
		if (onLiborPeriodDiscretization(evaluationTime)) {
			RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
			for(int i=0;i<dLdL.length;i++) dLdL[i][i]=new RandomVariable(1.0);
		    return dLdL;
		}
		
		// Calculate dLdL. It is a (n-1)x n Matrix!
		RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors+1];
		double swapTenorLength = model.getLiborPeriodDiscretization().getTimeStep(0); // Model must have same tenor as swap!
		double timeOfFirstLiborPriorToEval = getPreviousLiborTime(evaluationTime);
		int timeIndexAtEvaluationTime = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int timeIndexAtFirstLiborPriorToEval = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(timeOfFirstLiborPriorToEval);
		
		for(int liborIndex = 0; liborIndex <numberOfRemainingLibors; liborIndex++){
			double liborTime = evaluationTime+liborIndex*swapTenorLength; // t+j*\Delta T
		    int    previousLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(liborTime);
		    double previousLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex);
		    double firstNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+1);
		    double secondNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+2);
		    double factor1 = (secondNextLiborTime-(liborTime+swapTenorLength))/(secondNextLiborTime-firstNextLiborTime);
		    double factor2 = (liborTime-previousLiborTime)/(firstNextLiborTime-previousLiborTime);
		    int    timeIndex = liborIndex==0 ? timeIndexAtFirstLiborPriorToEval : timeIndexAtEvaluationTime;
		    // Get Libors. We take the average here to avoid calculation of pseudo inverse of a matrix of RV. Should be fixed.
		    RandomVariableInterface previousLibor = model.getLIBOR(timeIndex, previousLiborIndex);     
		    RandomVariableInterface nextLibor     = model.getLIBOR(timeIndex, previousLiborIndex + 1); 
		    RandomVariableInterface logInterpol = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).log().mult(-factor1);
		                            logInterpol = logInterpol.add(previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).log().mult(-factor2)).exp();
		    // Set derivatives
		    dLdL[liborIndex][liborIndex]   = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).mult(logInterpol).mult(1-factor2);// dLdL_i-1
		    dLdL[liborIndex][liborIndex+1] = previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).mult(logInterpol).mult(1-factor1);
		}
		
		// dLdL is (n-1) x n matrix. Get PseudoInverse for all paths and then put it back together as RV
		return getPseudoInverse(dLdL);
	}
	
	
	/**Calculates the row vector dV/dL
	 * 
	 * @param evaluationTime The time at which the forward sensistivity dVdL is calculated
	 * @return The forward sensisivity dVdL (as a row vector)
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getValueLiborSensitivities(double evaluationTime) throws CalculationException{
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		RandomVariableDifferentiableInterface numeraire = (RandomVariableDifferentiableInterface) model.getNumeraire(evaluationTime);
		
		// Calculate forward sensitivities
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime);
		int numberOfSensis = evaluationTime == getNextLiborTime(evaluationTime) ? numberOfRemainingLibors : numberOfRemainingLibors+1;
		RandomVariableInterface[] valueLiborSensitivities = new RandomVariableInterface[numberOfSensis];// exclude last libor
		int timeIndexAtEval = model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		
		// Set all entries of dVdL
		// Set dVdL for last libor which is already fixed (if applicable)
		int timeGridIndicator = 0;
		int lastLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		if(numberOfSensis!=numberOfRemainingLibors){
			timeGridIndicator = 1;
			double lastLiborTime = model.getLiborPeriodDiscretization().getTime(lastLiborIndex);
			RandomVariableInterface lastLibor = model.getLIBOR(model.getTimeDiscretization().getTimeIndex(lastLiborTime), lastLiborIndex);
			RandomVariableInterface dVdL = getProductValueDerivative(lastLibor);
			valueLiborSensitivities[0] = dVdL.mult(numeraire);
		}
		
		for(int liborIndex=lastLiborIndex+timeGridIndicator;liborIndex<model.getNumberOfLibors(); liborIndex++){
			RandomVariableInterface liborAtTimeIndex = model.getLIBOR(timeIndexAtEval, liborIndex);
		    RandomVariableInterface dVdL = getProductValueDerivative(liborAtTimeIndex);
		    valueLiborSensitivities[liborIndex-lastLiborIndex] = dVdL.mult(numeraire).getConditionalExpectation(conditionalExpectationOperator);
		}
		
		if(isUseTimeGridAdjustment){
		// Up to now dVdL is wrt the Libors on the LiborPeriodDiscretization. Adjust it such that we have dVdL wrt Libors starting at evaluationTime 
		RandomVariableInterface[][] dLdL = getLiborTimeGridAdjustment(evaluationTime);
		RandomVariableInterface[] dVdLAdjusted = multiply(valueLiborSensitivities,dLdL);
		return dVdLAdjusted; 
		} else return valueLiborSensitivities;
	}
	
	/**Calculates dV/dS where S are swap rates of the discount curve.
	 * 
	 * @param evaluationTime The time at which dVdS is calculated
	 * @return The row vector of sensitivities wrt swap rates from discount curve.
	 * @throws CalculationException 
	 */
	public RandomVariableInterface[] getDiscountCurveSensitivities(double evaluationTime) throws CalculationException{
		// We calculate dV/dP * dP/dS. dV/dP is at t=0 since the curve starts at 0 !?!
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		final double shift = 0.0001;
		
		// Remove first entry from pillars if it is at time 0.
		int index = discountCurvePillars[0]==0 ? 1 : 0;
		double[] pillars = new double[discountCurvePillars.length-index];
		for(int i=0;i<pillars.length;i++) pillars[i]=discountCurvePillars[i+index];
		
		RandomVariableInterface  value = currentProduct.getValue(evaluationTime, model); //.mult(model.getNumeraire(evaluationTime));
		Map<Long, RandomVariableInterface> gradientOfProduct = ((RandomVariableDifferentiableInterface) value).getGradient();
		int numberOfP = getNumberOfRemainingLibors(evaluationTime);
		
	
		int lastPillarIndex = evaluationTime>pillars[0] ? new TimeDiscretization(pillars).getTimeIndexNearestLessOrEqual(evaluationTime) : 0;
		
        double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		RandomVariableInterface[] dVdP = new RandomVariableInterface[pillars.length-lastPillarIndex];//numberOfP];
        DiscountCurve discountCurve = (DiscountCurve) model.getModel().getDiscountCurve();
        // Define new pillars. Another option is to calculate dV/dP wrt original Pillars, and do interpolation: dV/dP*dP/d\tilde{P}
//		TimeDiscretization curveTimes = new TimeDiscretization(evaluationTime, numberOfP, liborPeriodLength);
        TimeDiscretization curveTimes = new TimeDiscretization(pillars);
		RandomVariableInterface[] discountFactors = new RandomVariableInterface[pillars.length];//numberOfP+1];
		// get discount factors
		for(int i=0;i<pillars.length;i++) discountFactors[i]=discountCurve.getDiscountFactor(pillars[i]);//curveTimes.getTime(i));
		// dV(t)/dP(T_i;0)
		for(int i=lastPillarIndex;i<pillars.length;i++){
//			discountFactors[i]+=shift;
//			DiscountCurve newDiscountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve", pillars/*curveTimes.getAsDoubleArray()*/, discountFactors); 
//		    // get clone of LMM with shifted curve. 
//			Map<String,Object> dataModified = new HashMap<String,Object>();
//		    dataModified.put("discountCurve", newDiscountCurve);
//		    LIBORModelMonteCarloSimulationInterface newModel = (LIBORModelMonteCarloSimulation) model.getCloneWithModifiedData(dataModified);
//		    dVdP[i-lastPillarIndex] = currentProduct.getValue(evaluationTime, newModel).mult(newModel.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
//            dVdP[i-lastPillarIndex]=dVdP[i-lastPillarIndex].sub(value).div(shift).mult(discountCurve.getDiscountFactor(evaluationTime));
//		    discountFactors[i]-=shift;
			dVdP[i-lastPillarIndex] = getDerivative(gradientOfProduct, discountFactors[i]).getConditionalExpectation(conditionalExpectationOperator);
		}
		// Get dP(T_i;0)/dP(t+i\delta T;0): Linear interpolation on log value per time
		RandomVariableInterface[][] dPdP = new RandomVariableInterface[numberOfP][pillars.length-lastPillarIndex];
		for(int i=0;i<dPdP.length;i++){
			double discountTime = evaluationTime + (i+1) * liborPeriodLength;
			if(discountTime < pillars[0]) {
				double term = Math.pow(discountTime/pillars[0],2.0);
				dPdP[i][0]=discountFactors[0].invert().mult(term).mult(discountFactors[0].log().mult(term).exp());
				continue;
			}
			// Get upper and lower index
			int lowerIndex = curveTimes.getTimeIndexNearestLessOrEqual(discountTime); // as 0 is included in time discretization but nut in pillars
			lowerIndex = lowerIndex <0 ? 0 : lowerIndex;
			int upperIndex = lowerIndex+1;
			double delta = (discountTime-pillars[lowerIndex])/curveTimes.getTimeStep(lowerIndex);
			RandomVariableInterface summand1 = discountFactors[lowerIndex].log().mult((1-delta)/pillars[lowerIndex]);
			RandomVariableInterface summand2 = discountFactors[upperIndex].log().mult(delta/pillars[upperIndex]);
			RandomVariableInterface factor   = summand1.add(summand2).mult(discountTime).exp();
			//Math.exp(((1-delta)/pillars[lowerIndex]*Math.log(discountFactors[lowerIndex])+delta/pillars[upperIndex]*Math.log(discountFactors[upperIndex]))*discountTime);
			dPdP[i][lowerIndex-lastPillarIndex]=factor.div(discountFactors[lowerIndex]).mult((1-delta)/pillars[lowerIndex]*discountTime);
			dPdP[i][upperIndex-lastPillarIndex]=factor.div(discountFactors[upperIndex]).mult(delta/pillars[upperIndex]*discountTime);
		}
        dVdP = multiply(dVdP,getPseudoInverse(dPdP));
		RandomVariableInterface[][] dPdS = getBondSwapSensitivity(evaluationTime);
		RandomVariableInterface[] dVdS = multiply(dVdP,dPdS); // multiply(dVdN,dNdP)
		return getSensitivitiesOnBuckets(dVdS);
	}
	
	
	/**Calculates dPdS in a single curve context. Used for calculating sensis with respect to discount curve.
	 * 
	 * @param evaluationTime
	 * @return
	 * @throws CalculationException 
	 */
	private RandomVariableInterface[][] getBondSwapSensitivity(double evaluationTime) throws CalculationException{
		int numberOfBonds = getNumberOfRemainingLibors(evaluationTime);
		RandomVariableInterface sum= new RandomVariable(0.0);
		RandomVariableInterface[][] dSdP = new RandomVariableInterface[numberOfBonds][numberOfBonds];
		for(int bondIndex=0;bondIndex<dSdP[0].length;bondIndex++){
			RandomVariableInterface bond = model.getNumeraire(evaluationTime+(bondIndex+1)*0.5).invert().mult(model.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
		    sum = sum.add(bond);
		    for(int swapIndex=0;swapIndex<dSdP.length;swapIndex++){
		    	if(swapIndex<bondIndex) dSdP[swapIndex][bondIndex] = new RandomVariable(0.0);
		    	else if(swapIndex==bondIndex) dSdP[swapIndex][bondIndex] = sum.add(1.0).sub(bond).mult(-1.0).div(sum.squared());
		    	else dSdP[swapIndex][bondIndex] = bond.sub(1.0).div(sum.squared());    	
		    }
		} 
		return getPseudoInverse(dSdP); // PseudoInverse == Inverse for n x n matrix.
	}
	
	
	/**Calculate Pseudo Inverse of matrix of type RandomVariableInterface[][]
	 * 
	 * @param matrix The matrix for which the pseudo inverse is calculated
	 * @return The pseudo inverse of the matrix
	 */
    private RandomVariableInterface[][] getPseudoInverse(RandomVariableInterface[][] matrix){
    	double[][][] inv = new double[matrix[0].length][matrix.length][model.getNumberOfPaths()];
		double[][] matrixOnPath = new double[matrix.length][matrix[0].length];
		for(int pathIndex=0; pathIndex<model.getNumberOfPaths(); pathIndex++){
			// Get double[][] matrix on path
			for(int i=0;i<matrixOnPath.length;i++){
				for(int j=0;j<matrixOnPath[0].length;j++){
					matrixOnPath[i][j]=matrix[i][j]==null ? 0 : matrix[i][j].get(pathIndex);
				}
			}
		    // Get Pseudo Inverse 
		    RealMatrix pseudoInverse = new SingularValueDecomposition(MatrixUtils.createRealMatrix(matrixOnPath)).getSolver().getInverse();
		    for(int j=0;j<pseudoInverse.getColumnDimension();j++){
			    double[] columnValues = pseudoInverse.getColumn(j);
			    for(int i=0;i<pseudoInverse.getRowDimension();i++){
				    inv[i][j][pathIndex]= columnValues[i];
			    }
		    }
		}
		// Wrap to RandomVariableInterface[][]
		RandomVariableInterface[][] pseudoInverse = new RandomVariableInterface[matrix[0].length][matrix.length];
		for(int i=0;i<pseudoInverse.length; i++){
			for(int j=0;j<pseudoInverse[0].length; j++){
				pseudoInverse[i][j] = new RandomVariable(0.0 /*should be evaluationTime*/,inv[i][j]);
			}
		}
		return pseudoInverse;
    }
   

	/*
	 * Some useful functions 
	 */
	public static RandomVariableInterface[][] multiply(RandomVariableInterface[][] A,RandomVariableInterface[][] B){
		RandomVariableInterface[][] AB = new RandomVariableInterface[A.length][B.length];
		RandomVariableInterface ABproduct;
		for(int i=0;i<A.length;i++){
			for(int j=0; j<B.length; j++){
				AB[i][j] = new RandomVariable(0.0);
				for(int k=0;k<B.length;k++) {
					if(A[i][k]==null || B[k][j]==null) {ABproduct = new RandomVariable(0.0);}
					else {ABproduct = A[i][k].mult(B[k][j]);}
					AB[i][j]=AB[i][j].add(ABproduct);
				}
			}
		}
		return AB;
	}
	
	public static RandomVariableInterface[] multiply(RandomVariableInterface[] A,RandomVariableInterface[][] B){
		RandomVariableInterface[] AB = new RandomVariableInterface[B[0].length];
		RandomVariableInterface ABproduct;
		for(int i=0;i<B[0].length;i++){
				AB[i] = new RandomVariable(0.0);
				for(int k=0;k<A.length;k++) {
					if(A[k]==null || B[k][i]==null) {ABproduct = new RandomVariable(0.0);}
					else {ABproduct = A[k].mult(B[k][i]);}
					AB[i]=AB[i].add(ABproduct);
				}
		}
		return AB;
	}
	
	/**Calculates the derivative of the current portfolio product with respect to the specified parameter, dV/dX 
	 * 
	 * @param parameter The parameter with respect to which the derivative is calculated
	 * @return dV/dX 
	 * @throws CalculationException
	 */
	private RandomVariableInterface getProductValueDerivative(RandomVariableInterface parameter) throws CalculationException{
		if(gradientOfProduct==null) setGradient(this.currentProduct);
		RandomVariableInterface derivative = this.gradientOfProduct.get(((RandomVariableDifferentiableInterface)parameter).getID());
		return derivative==null ? new RandomVariable(0.0) : derivative;
	}
	
	private RandomVariableInterface getDerivative(Map<Long, RandomVariableInterface> gradient, RandomVariableInterface parameter){
		RandomVariableInterface derivative = gradient.get(((RandomVariableDifferentiableInterface)parameter).getID());
		return derivative==null ? new RandomVariable(0.0) : derivative;
	}
	
	private int getNumberOfRemainingLibors(double evaluationTime){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getNumberOfLibors()-nextLiborIndex;
	}
	
	private double getNextLiborTime(double evaluationTime){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex);
	}
	
	private double getPreviousLiborTime(double evaluationTime){
		if(evaluationTime==0) return 0.0;
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex-1);
	}
	
	/*
	// This function is used to calculate the input of conditional expectation (i.e. the F_t measurable RV).
	// Currently not used: We take cond. exp. wrt only short libor and long libor.
	private RandomVariableInterface[] getRemainingLibors(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		// Ask the model for its discretisation
	    int timeIndex	= model.getTimeDiscretization().getTimeIndexNearestLessOrEqual(evaluationTime);
		// Get all libors at timeIndex which are not yet fixed (others null)
		ArrayList<RandomVariableInterface> liborsAtTimeIndex = new ArrayList<RandomVariableInterface>();
		int firstLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		if(model.getLiborPeriodDiscretization().getTime(firstLiborIndex)>evaluationTime) liborsAtTimeIndex.add(model.getLIBOR(evaluationTime, evaluationTime, model.getLiborPeriodDiscretization().getTime(firstLiborIndex)));
		for(int i=firstLiborIndex;i<model.getNumberOfLibors();i++) {
			    liborsAtTimeIndex.add(model.getLIBOR(timeIndex,i));
		}
		return liborsAtTimeIndex.toArray(new RandomVariableInterface[liborsAtTimeIndex.size()]);
	}
	*/	
		
	private static ArrayList<RandomVariableInterface> getRegressionBasisFunctions(RandomVariableInterface[] libors, int order) {
		ArrayList<RandomVariableInterface> basisFunctions = new ArrayList<RandomVariableInterface>();
		// Create basis functions - here: 1, S, S^2, S^3, S^4
		for(int liborIndex=0; liborIndex<libors.length;liborIndex++){
		  for(int powerOfRegressionMonomial=0; powerOfRegressionMonomial<=order; powerOfRegressionMonomial++) {
			  basisFunctions.add(libors[liborIndex].pow(powerOfRegressionMonomial));
		  }
		}
		return basisFunctions;
	}
	
	/*
	 *  Getters and Setters
	 */
	public AbstractLIBORMonteCarloProduct[] getPortfolioProducts() {
		return portfolioProducts;
	}
	
	public AbstractLIBORMonteCarloProduct getPortfolioProduct(int index) {
		return portfolioProducts[index];
	}

	public String getDiscountCurveName() {
		return model.getModel().getDiscountCurve().getName();
	}


	public String getForwardCurveName() {
		return model.getModel().getForwardRateCurve().getName();
	}
	
	public WeightToLiborAdjustmentMethod getLiborWeightMethod(){
		return this.liborWeightMethod;
	}
	
	public LIBORModelMonteCarloSimulationInterface getModel(){
		return this.model;
	}
	
	public void setProduct(int productIndex){
		this.currentProduct = portfolioProducts[productIndex];
	}
	
	public AbstractLIBORMonteCarloProduct getProduct(){
		return this.currentProduct;
	}
	
	private boolean onLiborPeriodDiscretization(double evaluationTime){
		return (evaluationTime == getNextLiborTime(evaluationTime));
	}
	
	private void setConditionalExpectationOperator(double evaluationTime) throws CalculationException{
		// Create a conditional expectation estimator with some basis functions (predictor variables) for conditional expectation estimation.
        RandomVariableInterface[] regressor = new RandomVariableInterface[2];
        regressor[0]= model.getLIBOR(evaluationTime, evaluationTime,evaluationTime+model.getLiborPeriodDiscretization().getTimeStep(0));
		regressor[1]= model.getLIBOR(evaluationTime, evaluationTime, model.getLiborPeriodDiscretization().getTime(model.getNumberOfLibors()-1));
       	ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(regressor, 2);
//		Alternative definition of regressors
//		RandomVariableInterface[] libors = getRemainingLibors(evaluationTime, model);
//		ArrayList<RandomVariableInterface> basisFunctions = getRegressionBasisFunctions(libors, 1 /*polyNomialOrder*/);
//		ConditionalExpectationEstimatorInterface conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));
       	this.conditionalExpectationOperator = new MonteCarloConditionalExpectationRegression(basisFunctions.toArray(new RandomVariableInterface[0]));

	}
	
	private void setGradient(AbstractLIBORMonteCarloProduct product) throws CalculationException{
		RandomVariableDifferentiableInterface productValue = (RandomVariableDifferentiableInterface) product.getValue(model.getTime(0), model);
		Map<Long, RandomVariableInterface> gradientOfProduct = productValue.getGradient();
		this.gradientOfProduct = gradientOfProduct;
	}

	// Delete later..
	public void setUseTimeGridAdjustment(boolean method){ 
		this.isUseTimeGridAdjustment = method;
	}
	// Delete later..
	public void setIgnoreDiscountCurve(boolean method){ 
		this.isIgnoreDiscountCurve = method;
	}
	
	private RandomVariableInterface[][] getPseudoInverse(double[][] matrix){
		RealMatrix pseudoInverse = new SingularValueDecomposition(MatrixUtils.createRealMatrix(matrix)).getSolver().getInverse();
		RandomVariableInterface[][] inv = new RandomVariableInterface[matrix[0].length][matrix.length];
		for(int j=0;j<pseudoInverse.getColumnDimension();j++){
		    double[] columnValues = pseudoInverse.getColumn(j);
		    for(int i=0;i<pseudoInverse.getRowDimension();i++){
			    inv[i][j]= new RandomVariable(columnValues[i]);
		    }		    
	    }
		return inv;
	}

}
