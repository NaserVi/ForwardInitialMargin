package initialmargin.isdasimm;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.montecarlo.AbstractMonteCarloProduct;
import net.finmath.montecarlo.MonteCarloSimulationInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.automaticdifferentiation.RandomVariableDifferentiableInterface;
import net.finmath.montecarlo.conditionalexpectation.MonteCarloConditionalExpectationRegression;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.ConditionalExpectationEstimatorInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;


import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class SIMMSchemeMain extends AbstractLIBORMonteCarloProduct{

	
	// The class SIMM contains two other classes of which SIMM class variables are generated.
	public static class ParameterCollection{
        public ParameterCollection(){//constructor
        }
        
        public Double[][]               CrossRiskClassCorrelationMatrix;

        public Map<String,String>       MapFXCategory;

        final public String[]           CreditMaturityBuckets = {"1y","2y","3y","5y","10y"};
        final public String[]           IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};

        public Double             IRCorrelationCrossCurrency = .27;

        public Map<String,String>       IRCurrencyMap;

        public Map<String,Double[][] > MapRiskClassCorrelationIntraBucketMap;
        public Map<String,Double[][] > MapRiskClassCorrelationCrossBucketMap;
        public Map<String,Map<String,Map<String,Double[][]> > > MapRiskClassThresholdMap;
        public Map<String,Map<String,Map<String,Double[][]> > > MapRiskClassRiskweightMap;


        public void setIRCurrencyMap(Map<String, String> IRCurrencyMap) {
            this.IRCurrencyMap = IRCurrencyMap;
        }

        public ParameterCollection setCrossRiskClassCorrelationMatrix(Double[][] crossRiskClassCorrelationMatrix) {
            CrossRiskClassCorrelationMatrix = crossRiskClassCorrelationMatrix;
            return this;
        }

    }
	
	
	public class ClassifiedSIMMProduct {
	    
		private AbstractMonteCarloProduct product; // can be IR or other (equity, commodity)
		private Map<Long, RandomVariableInterface> gradientOfProduct; // Same for all evaluationTimes; Is reset for different products
		private double lastEvaluationTime;
		// Product classification within ISDA SIMM
		private String   productClass;
		private String[] riskClass; 
		private String[] curveIndexNames;   // e.g. OIS & Libor6m
		private String   currency;
		private boolean  hasOptionality; // e.g. Swap has no curvature risk
		private String   bucketKey;
		
		private HashMap<String/*RiskClass*/,HashMap<String/*curveIndexName*/,
		                HashMap<String/*maturityBucket*/,RandomVariableInterface>>> deltaSensitivities = null; // currently only for InterestRate riskClass
		private RandomVariableInterface vegaSensitivity=null;
		
		public ClassifiedSIMMProduct(AbstractMonteCarloProduct product,
				                 String   productClass,
				                 String[] riskClass, // One product may contribute to several risk Classes
				                 String[] curveIndexNames,
				                 String   currency,
				                 String   bucketKey,
				                 boolean  hasOptionality){
			this.product=product;
			this.productClass = productClass; 
			this.riskClass = riskClass;
			this.curveIndexNames = curveIndexNames;
			this.currency=currency;
			this.hasOptionality = hasOptionality;
			this.bucketKey = bucketKey;
		
		}
		
		public RandomVariableInterface getSensitivity(String productClass, 
				                                      String riskClass, 
				                                      String maturityBucket, // only for IR and Credit risk class, null otherwise
				                                      String curveIndexName, // null if riskClass not IR
				                                      String bucketKey,      // currency for IR otherwise bucket nr.
				                                      String riskType, double evaluationTime, MonteCarloSimulationInterface model) throws CalculationException, SolverException, CloneNotSupportedException{
			
			RandomVariableInterface result=null;
			if(!hasOptionality && riskType!="delta") return new RandomVariable(0.0);
			if(gradientOfProduct==null) setGradient(model); // needs to be set only once
			if(evaluationTime!=lastEvaluationTime) {
				deltaSensitivities.clear();
				vegaSensitivity = null;
			}
			if(productClass==this.productClass && 
			   Arrays.asList(this.riskClass).contains(riskClass)){
			   
			   switch(riskType){
			      case("delta"): 
				      switch(riskClass){
				          case("InterestRate"):
				 
				              if(Arrays.asList(this.curveIndexNames).contains(curveIndexName) & bucketKey==this.currency){
					             // There exists a sensitivity. Check if the sensitivities (on all maturityBuckets) have already been calculated for given riskClass and riskType)
					             if(!deltaSensitivities.containsKey(riskClass) || !deltaSensitivities.get(riskClass).containsKey(curveIndexName)){
						            // The sensitivities need to be calculated for the given riskClass and riskType
						            RandomVariableInterface[] maturityBucketSensis = doCalculateDeltaSensitivitiesIR(curveIndexName, this, evaluationTime); // currently only for riskClass IR
						  
						            HashMap<String,HashMap<String,RandomVariableInterface>> curveIndexNameSensiMap = new HashMap<String,HashMap<String,RandomVariableInterface>>();
						            HashMap<String,RandomVariableInterface> bucketSensitivities = new HashMap<String,RandomVariableInterface>();
						  
						            for(int i=0;i<parameterCollection.IRMaturityBuckets.length;i++) bucketSensitivities.put(parameterCollection.IRMaturityBuckets[i], maturityBucketSensis[i]);
						            curveIndexNameSensiMap.put(curveIndexName,bucketSensitivities);
						            deltaSensitivities.put(riskClass, curveIndexNameSensiMap);
					             }
					             result = deltaSensitivities.get(riskClass).get(curveIndexName).get(maturityBucket);
				               } else result = new RandomVariable(0.0); // There exists no delta Sensi for risk Class InterestRate
				           break;
			               case("CreditQ"):
			               case("CreditNonQ"):
			               case("FX"):
			               case("Commodity"):
			               case("Equity"): result = null;
				       }
			      
			      case("vega"): 
				  case("curvature"):
				       switch(riskClass){
				           case("InterestRate"): //if(vegaSensitivity!=null) vegaSensitivity = getVegaSensitivityIR(curveIndexNames, product, evaluationTime, model);
				           case("CreditQ"):
				           case("CreditNonQ"):
				           case("FX"):
				           case("Commodity"):
				           case("Equity"): result=null;
				       }
			   }
			}
		  return result;
		} // end getSensitivity()


		private void setGradient(MonteCarloSimulationInterface model) throws CalculationException{
			RandomVariableDifferentiableInterface productValue = (RandomVariableDifferentiableInterface) product.getValue(model.getTime(0), model);
			Map<Long, RandomVariableInterface> gradientOfProduct = productValue.getGradient();
			this.gradientOfProduct = gradientOfProduct;
		}
		
		public String getProductClass(){
			return productClass;
		}
		
		public String[] getRiskClasses(){
			return riskClass;
		}
		
		public String[] getCurveIndexNames(){
			return curveIndexNames;
		}

		public Map<Long, RandomVariableInterface> getGradient() {
			return gradientOfProduct;
		}
		public AbstractMonteCarloProduct getProduct(){
			return product;
		}
	} // end class ClassifiedSIMMProduct
	
	/*
	 *  Commencement of actual SIMM scheme 
	 */

    Map<AbstractMap.SimpleEntry<Double,SIMMSchemeSensitivitySet.Key>,RandomVariableInterface> CacheNetSensitivities;

    public Map<String,Double>   resultMap;

    private ClassifiedSIMMProduct[] portfolioProducts;
    
    private final LIBORModelMonteCarloSimulationInterface model;
    private RandomVariableInterface[][] riskWeightToLiborAdjustments;
    private boolean isUseTimeGridAdjustment;
    private boolean isIgnoreDiscountCurve=false;  // can be discarded later.. just to check how discount curve influences IM
	private double[] discountCurvePillars = {0.5 , 1.0, 2.0, 5.0, 30.0};// shouldn't the discount curve know its pillars ?
    private ConditionalExpectationEstimatorInterface conditionalExpectationOperator;

    private ParameterCollection parameterCollection;
    private String[] productClassKeys;
    private String[] riskClassKeys;
    private String[] IRCurveIndexNames;

    private String calculationCCY;
    
    public enum WeightToLiborAdjustmentMethod{
		Constant,  //Sets dL/dS(t=0) for all forward IM times, i.e. leave the weight adjustment dL/dS constant
		Stochastic //Calculate dL/dS(t) for all forward IM times, i.e. (weakly) stochastic weight adjustment 
	}
    private WeightToLiborAdjustmentMethod liborWeightMethod;

    // SIMM constructor
    public SIMMSchemeMain(ClassifiedSIMMProduct[] portfolioProducts, 
    					  LIBORModelMonteCarloSimulationInterface model, 
    					  ParameterCollection parameterCollection, 
    					  String calculationCCY,
    					  WeightToLiborAdjustmentMethod method) throws CalculationException{
        this.resultMap = new HashMap<>();
        this.calculationCCY = calculationCCY;
        this.model = model;
        CacheNetSensitivities = new HashMap<>();
        this.portfolioProducts = portfolioProducts;
        this.parameterCollection = parameterCollection;
        this.liborWeightMethod=method;
        
        // Screen portfolio products for relevant product classes, risk classes and curveIndexNames
        ArrayList<String> relevantProductClasses = new ArrayList<String>();
        ArrayList<String> relevantCurveIndices = new ArrayList<String>();
        ArrayList<String> relevantRiskClasses = new ArrayList<String>();
        
        for(int productIndex=0;productIndex<portfolioProducts.length;productIndex++){
        	String productClass = portfolioProducts[productIndex].getProductClass();
            if(!relevantProductClasses.contains(productClass)) relevantProductClasses.add(productClass);
            String[] curveIndex = portfolioProducts[productIndex].getCurveIndexNames();
            for(int i=0;i<curveIndex.length;i++){
            	if(!relevantCurveIndices.contains(curveIndex[i])) relevantCurveIndices.add(curveIndex[i]);
            }
            String[] riskClass = portfolioProducts[productIndex].getRiskClasses();
            for(int i=0;i<riskClass.length;i++){
            	if(!relevantRiskClasses.contains(riskClass[i])) relevantRiskClasses.add(curveIndex[i]);
            }
        }
        this.productClassKeys = (String[])relevantProductClasses.toArray();
        this.riskClassKeys   = (String[])relevantRiskClasses.toArray();
        this.IRCurveIndexNames = (String[])relevantCurveIndices.toArray();
        
        // Get Weight Adjustment for Libor Sensitivities if they remain constant 
        if(liborWeightMethod == WeightToLiborAdjustmentMethod.Constant){
     			this.riskWeightToLiborAdjustments = getLiborSwapSensitivities(0.0 /*evaluationTime*/);
        } 
    }



    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
        // Set liborWeightAdjustments dLdS
    	if(liborWeightMethod == WeightToLiborAdjustmentMethod.Stochastic) this.riskWeightToLiborAdjustments = getLiborSwapSensitivities(evaluationTime);
    	setConditionalExpectationOperator(evaluationTime);
    	RandomVariableInterface SIMMValue = null;
        for ( String productClass : productClassKeys) { // RatesFX, Credit etc.
            RandomVariableInterface SIMMProductValue = this.getSIMMProduct(productClass,evaluationTime);
            SIMMValue = SIMMValue == null ? SIMMProductValue : SIMMValue.add(SIMMProductValue);
        }
        return SIMMValue;
    }
    
    
    // for all times!
    public Map<String,Double> getResultMap() throws CalculationException{
        if (this.resultMap.entrySet().size()== 0)
            this.getValue(0.0, model);
        return this.resultMap;
    }

    // returns IM for productClass e.g. RatesFX
    public RandomVariableInterface      getSIMMProduct(String productClass, double atTime){

        Set<String> riskClassList = Stream.of(riskClassKeys).collect(Collectors.toSet());

        // contributions of risk classes to IM within that product class
        RandomVariableInterface[] contributions = new RandomVariableInterface[riskClassKeys.length];
        int i = 0;
        int pathDim = 1;//this.portfolioProducts[0].get//this.nettingset.getActiveProduct(this.nettingset.getActiveProductKeys()[0]).getSensitivitySet(atTime).getPathDimension();
        for ( String iRiskClass : riskClassKeys)
        {
            if ( riskClassList.contains(iRiskClass)) { // riskClassList == iRiskClass?
                RandomVariableInterface iIM = this.getIMForRiskClass(iRiskClass, productClass, atTime);
                contributions[i] = iIM;
            }
            else
                contributions[i] = new RandomVariable(atTime,pathDim,0.0);
            i++;
        }

        RandomVariableInterface simmProductClass = SIMMSchemeMain.getVarianceCovarianceAggregation(contributions,parameterCollection.CrossRiskClassCorrelationMatrix);
        resultMap.put(productClass,simmProductClass.getAverage());
        return simmProductClass;
    }


    public RandomVariableInterface getIMForRiskClass(String riskClassKey,String productClass, double atTime){
        RandomVariableInterface    deltaMargin = this.getDeltaMargin(riskClassKey,productClass,atTime);
        RandomVariableInterface    vegaMargin = this.getVegaMargin(riskClassKey,productClass,atTime);
        //RandomVariableInterface    curatureMargin = this.getDeltaMargin(riskClassKey,productClass, atTime);
        resultMap.put(productClass+"-"+riskClassKey+"-DeltaMargin",deltaMargin.getAverage());
        return deltaMargin.add(vegaMargin);//.add(vegaMargin).add(curatureMargin);
    }

    
    public RandomVariableInterface getDeltaMargin(String riskClassKey,String productClassKey, double atTime){
        RandomVariableInterface deltaMargin = null;
        //DeltaMarginSchemeNonIR DeltaScheme = new DeltaMarginSchemeNonIR(this,"Risk_IRCurve",productClassKey,
        //                                        this.riskClassRiskWeightMap.get(riskClassKey),this.riskClassCorrelationMap.get(riskClassKey),this.riskClassThresholdMap.get(riskClassKey));

        String riskType = "delta";
        if ( riskClassKey.equals("InterestRate"))
        {
            SIMMSchemeIRDelta DeltaScheme = new SIMMSchemeIRDelta(this,productClassKey);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
        /*else if ( riskClassKey.equals("CreditQ")){
            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditQ",productClassKey,riskType);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
        else if ( riskClassKey.equals("CreditNonQ")){
            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditNonQ",productClassKey,riskType);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
        else if ( riskClassKey.equals("Equity")){
            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Equity",productClassKey,riskType);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
        else if ( riskClassKey.equals("Commodity")){
            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Commodity",productClassKey,riskType);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
        else if ( riskClassKey.equals("FX")){
            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"FX",productClassKey,riskType);
            deltaMargin = DeltaScheme.getValue(atTime);
        }*/

        //MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,riskClassKey,productClassKey,riskType);
        return deltaMargin;
    }


    public RandomVariableInterface      getVegaMargin(String riskClassKey,String productClassKey,double atTime){
    	throw new RuntimeException();
        //MarginSchemeDeltaVega VegaScheme = new MarginSchemeDeltaVega(this,riskClassKey,productClassKey,"vega");
        //return VegaScheme.getValue(atTime);
    }

    public RandomVariableInterface      getCurvatureMargin(String riskClassKey,String productClassKey,double atTime){
        throw new RuntimeException();
    }



    public  static  RandomVariableInterface  getVarianceCovarianceAggregation(RandomVariableInterface[] contributions, Double[][] correlation){
        int i = 0;
        RandomVariableInterface value = null;
        for (RandomVariableInterface contribution1 : contributions) {
            int j=0;
            if ( contribution1!=null) {
                value = value == null ? contribution1.squared() : value.add(contribution1.squared());
                for (RandomVariableInterface contribution2 : contributions) {
                    if (contribution2 != null && i != j) {
                        RandomVariableInterface contribution = contribution1.mult(contribution2).mult(correlation[i][j]);
                        value = value == null ? contribution : value.add(contribution);
                    }
                    j++;
                }
            }
            i++;
        }
        value = value.sqrt();
        return value;
    }

    public  static  RandomVariableInterface  getVarianceCovarianceAggregation(RandomVariableInterface[] contributions, Double correlation){
        int i = 0;
        RandomVariableInterface value = null;
        for (RandomVariableInterface contribution1 : contributions) {
            int j=0;
            for (RandomVariableInterface contribution2 : contributions) {
                RandomVariableInterface contribution = null;
                if ( i!=j)
                    contribution= contribution1.mult(contribution2).mult(correlation);
                else
                    contribution= contribution1.mult(contribution2);
                value = value == null ? contribution : value.add(contribution);
                j++;
            }
            i++;
        }
        value = value.sqrt();
        return value;
    }

    // BUCKET IS CURRENCY FOR IR   risk Factor = index Name (e.g. Libor6m)
    public RandomVariableInterface   getNetSensitivity(String productClassKey, String riskClassKey, String maturityBucket,String riskFactor, String bucketKey,String riskType, double atTime) {

        //SIMMSchemeSensitivitySet.Key key = this.getKeySensitivitySet(iRateTenor,riskFactor,bucketKey,productClassKey,riskClassKey,riskType);
        
        RandomVariableInterface isdasimmsensiofAllProducts = Stream.of(this.portfolioProducts).map(
        		product->{try{
        					return product.getSensitivity(productClassKey, 
        							                      riskClassKey, 
                                                          maturityBucket, // only for IR and Credit risk class, null otherwise
                                                          riskFactor,     // CurveIndexName  
                                                          bucketKey,      // currency for IR otherwise null ?
                                                          riskType, atTime, model);
        					 }
        					catch(Exception e){ return null;
        						              }
        		          }).reduce(null, (e1, e2) -> {
                             if (e1 == null) return e2;
                             if (e2 == null) return e1;
                             return e1.add(e2);
                          });
        
        return isdasimmsensiofAllProducts;

    }
    
    
    /**Calculate the sensitivities dV/dS with respect to all swap rates for given product and curve. This applies to the risk class Interest Rates only.
	 * 
	 * @param curveIndexName The name of the curve to be considered (OIS, LiborXm)
	 * @param product The product whose sensitivity is to be considered.
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return The sensitivities dV/dS i.e. with respect to swap rates.
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	private RandomVariableInterface[] doCalculateDeltaSensitivitiesIR(String curveIndexName, 
                                                                      ClassifiedSIMMProduct product, 
                                                                      double evaluationTime) throws SolverException, CloneNotSupportedException, CalculationException{
		RandomVariableInterface[] delta;
		if(curveIndexName!="OIS"){ 
		   RandomVariableInterface[] dVdL = getValueLiborSensitivities(product, evaluationTime);
		   // the following line will be removed later. Just checking how timeGridAdjustment affects the result
		   int timeGridIndicator = 0; if(!isUseTimeGridAdjustment && !onLiborPeriodDiscretization(evaluationTime)) timeGridIndicator = 1;
		
		   delta = new RandomVariableInterface[dVdL.length-timeGridIndicator];
		   RandomVariableInterface[][] dLdS = this.riskWeightToLiborAdjustments;
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
		 } else { // CurveIndexName == OIS
			delta = getDiscountCurveSensitivities(product,evaluationTime);
		 }
		 return getSensitivitiesOnBuckets(delta, "InterestRate");
	}
	
	/**Performs rebucketing of sensitivities to the SIMM buckets by linear interpolation.
	 * 
	 * @param sensitivities The sensitivities wrt swap rates dV/dS
	 * @param riskClass The risk class
	 * @return The sensitivities on the SIMM maturity buckets
	 */
	private RandomVariableInterface[] getSensitivitiesOnBuckets(RandomVariableInterface[] sensitivities, String riskClass){
		//rebucketing to SIMM structure(buckets: 2w, 1m, 3m, 6m, 1y, 2y, 3y, 5y, 10y, 15y, 20y, 30y)	
		int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};
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
	
	
	/**Calculates dL/dS 
	 * 
	 * @param evaluationTime The time at which the sensitivity is calculated
	 * @return The matrix dL/dS 
	 * @throws CalculationException
	 */
	public RandomVariableInterface[][] getLiborSwapSensitivities(double evaluationTime) throws CalculationException{
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		RandomVariableInterface[][] dLdS=null;
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		
	    // Get index of first Libor starting >= evaluationTime
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int numberOfRemainingLibors = model.getNumberOfLibors()-nextLiborIndex;
		dLdS = new RandomVariableInterface [numberOfRemainingLibors][numberOfRemainingLibors];
					
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
	public RandomVariableInterface[] getValueLiborSensitivities(ClassifiedSIMMProduct product, double evaluationTime) throws CalculationException{
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
			RandomVariableInterface dVdL = getProductValueDerivative(product,lastLibor);
			valueLiborSensitivities[0] = dVdL.mult(numeraire);
		}
		
		for(int liborIndex=lastLiborIndex+timeGridIndicator;liborIndex<model.getNumberOfLibors(); liborIndex++){
			RandomVariableInterface liborAtTimeIndex = model.getLIBOR(timeIndexAtEval, liborIndex);
		    RandomVariableInterface dVdL = getProductValueDerivative(product,liborAtTimeIndex);
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
	public RandomVariableInterface[] getDiscountCurveSensitivities(ClassifiedSIMMProduct product, double evaluationTime) throws CalculationException{
		// We calculate dV/dP * dP/dS. dV/dP is at t=0 since the curve starts at 0 !?!
		if(conditionalExpectationOperator==null) setConditionalExpectationOperator(evaluationTime);
		final double shift = 0.0001;
		
		// Remove first entry from pillars if it is at time 0.
		int index = discountCurvePillars[0]==0 ? 1 : 0;
		double[] pillars = new double[discountCurvePillars.length-index];
		for(int i=0;i<pillars.length;i++) pillars[i]=discountCurvePillars[i+index];
		
		RandomVariableInterface  value = product.getProduct().getValue(evaluationTime, model).mult(model.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
		int numberOfP = getNumberOfRemainingLibors(evaluationTime);
		
	
		int lastPillarIndex = evaluationTime>pillars[0] ? new TimeDiscretization(pillars).getTimeIndexNearestLessOrEqual(evaluationTime) : 0;
		
        double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		RandomVariableInterface[] dVdP = new RandomVariableInterface[pillars.length-lastPillarIndex];//numberOfP];
        DiscountCurve discountCurve = (DiscountCurve) model.getModel().getDiscountCurve();
        // Define new pillars. Another option is to calculate dV/dP wrt original Pillars, and do interpolation: dV/dP*dP/d\tilde{P}
//		TimeDiscretization curveTimes = new TimeDiscretization(evaluationTime, numberOfP, liborPeriodLength);
        TimeDiscretization curveTimes = new TimeDiscretization(pillars);
		double[] discountFactors = new double[pillars.length];//numberOfP+1];
		// get discount factors
		for(int i=0;i<pillars.length;i++) discountFactors[i]=discountCurve.getDiscountFactor(pillars[i]);//curveTimes.getTime(i));
		// dV(t)/dP(T_i;0)
		for(int i=lastPillarIndex;i<pillars.length;i++){
			discountFactors[i]+=shift;
			DiscountCurve newDiscountCurve = DiscountCurve.createDiscountCurveFromDiscountFactors("discountCurve", pillars/*curveTimes.getAsDoubleArray()*/, discountFactors); 
		    // get clone of LMM with shifted curve. 
			Map<String,Object> dataModified = new HashMap<String,Object>();
		    dataModified.put("discountCurve", newDiscountCurve);
		    LIBORModelMonteCarloSimulationInterface newModel = (LIBORModelMonteCarloSimulation) model.getCloneWithModifiedData(dataModified);
		    dVdP[i-lastPillarIndex] = product.getProduct().getValue(evaluationTime, newModel).mult(newModel.getNumeraire(evaluationTime)).getConditionalExpectation(conditionalExpectationOperator);
            dVdP[i-lastPillarIndex]=dVdP[i-lastPillarIndex].sub(value).div(shift).mult(discountCurve.getDiscountFactor(evaluationTime));
		    discountFactors[i]-=shift;
		}
		// Get dP(T_i;0)/dP(t+i\delta T;0): Linear interpolation on log value per time
		double[][] dPdP = new double[numberOfP][pillars.length-lastPillarIndex];
		for(int i=0;i<dPdP.length;i++){
			double discountTime = evaluationTime + (i+1) * liborPeriodLength;
			if(discountTime < pillars[0]) {
				double term = Math.pow(discountTime/pillars[0],2.0);
				dPdP[i][0]=term/discountFactors[0]*Math.exp(term*Math.log(discountFactors[0]));
				continue;
			}
			// Get upper and lower index
			int lowerIndex = curveTimes.getTimeIndexNearestLessOrEqual(discountTime); // as 0 is included in time discretization but nut in pillars
			lowerIndex = lowerIndex <0 ? 0 : lowerIndex;
			int upperIndex = lowerIndex+1;
			double delta = (discountTime-pillars[lowerIndex])/curveTimes.getTimeStep(lowerIndex);
			double factor = Math.exp(((1-delta)/pillars[lowerIndex]*Math.log(discountFactors[lowerIndex])+delta/pillars[upperIndex]*Math.log(discountFactors[upperIndex]))*discountTime);
			dPdP[i][lowerIndex-lastPillarIndex]=(1-delta)/pillars[lowerIndex]*discountTime*factor/discountFactors[lowerIndex];
			dPdP[i][upperIndex-lastPillarIndex]=delta/pillars[upperIndex]*discountTime*factor/discountFactors[upperIndex];
		}
        dVdP = multiply(dVdP,getPseudoInverse(dPdP));
		RandomVariableInterface[][] dPdS = getBondSwapSensitivity(evaluationTime);
		RandomVariableInterface[] dVdS = multiply(dVdP,dPdS); // multiply(dVdN,dNdP)
		return getSensitivitiesOnBuckets(dVdS, "InterestRate");
	}
	
	
	/**Calculates dPdS in a single curve context. Used for calculating sensis with respect to discount curve.
	 * 
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @return The sensitivity of the discount curve (bonds) wrt to swap rates of the same curve.
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
	
	//----------------------------------------------------------------------------------------------------------------------------------
	// Some auxiliary functions
	//----------------------------------------------------------------------------------------------------------------------------------

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
	private RandomVariableInterface getProductValueDerivative(ClassifiedSIMMProduct product, RandomVariableInterface parameter) throws CalculationException{
		if(product.getGradient()==null) product.setGradient(model); // This should never be the case.
		RandomVariableInterface derivative = product.getGradient().get(((RandomVariableDifferentiableInterface)parameter).getID());
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
	public ClassifiedSIMMProduct[] getPortfolioProducts() {
		return portfolioProducts;
	}
	
	public AbstractMonteCarloProduct getPortfolioProduct(int index) {
		return portfolioProducts[index].getProduct();
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
    
    public String[] getIRCurveIndexNames(){
    	return this.IRCurveIndexNames;
    }
    
    public String    getCalculationCCY(){
        return this.calculationCCY;
    }

    
    // I think we do not need this ?!
    public int getPathDimension(){
        return 1;
    }

    public ParameterCollection getParameterCollection() {
        return parameterCollection;
    }

    public void setParameterCollection(ParameterCollection parameterCollection) {
        this.parameterCollection = parameterCollection;
    }

    
}