package initialmargin.isdasimm.aggregationscheme;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.ArrayUtils;

import net.finmath.exception.CalculationException;
import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.products.*;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.optimizer.OptimizerFactoryInterface;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.OptimizerInterface;
import net.finmath.optimizer.OptimizerInterface.ObjectiveFunction;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;

/**
 * @author Peter Kohl-Landgraf
 */

public class CalculationSchemeInitialMarginISDA {

	public class ParameterCollection{
	     public ParameterCollection(){//hard values inserted by Mario Viehmann
	     	
	     	// Set correlationMatrixWithinSubCurve
	     	Double[][] corrMatrix = new Double[correlationMatrixWithinSubCurve.length*5+2][correlationMatrixWithinSubCurve.length*5+2];
	     	//for(int i=0;i<curve.length;i++)curve[i]=ArrayUtils.toObject(this.correlationMatrixWithinSubCurve[i]);
	     	for(int i=0; i<corrMatrix.length; i++){
	     		int rowIndex = i/correlationMatrixWithinSubCurve.length;
	     		int row = i % correlationMatrixWithinSubCurve.length;
	     		for(int j=0; j<corrMatrix.length; j++){
	     			if(i==corrMatrix.length-1 || i==corrMatrix.length-2 || j==corrMatrix.length-1 || j==corrMatrix.length-2) {
	     				corrMatrix[i][j]=new Double(0.0);
	     				continue;
	     			}
	     		    
	     		    int colIndex = j/correlationMatrixWithinSubCurve.length;
	     		    double curveCorrelation = rowIndex == colIndex ? 1.0 : 0.982;
	     		    
	     		    int col = j % correlationMatrixWithinSubCurve.length;
	     		    corrMatrix[i][j]= correlationMatrixWithinSubCurve[row][col]*curveCorrelation;
	     		}
	     	}
	     	
	     	this.MapRiskClassCorrelationIntraBucketMap.put("InterestRate", corrMatrix);
	     	
	     	// Set IR Currency Map
	     	this.IRCurrencyMap = new HashMap<String,String>();
	     	this.IRCurrencyMap.put("EUR","Regular_Volatility_Currencies");
	     	// Set RiskClassThresholdMap
	     	Map<String, Double[][]> innerMap= new HashMap<String, Double[][]>(); Double[][] value = new Double[1][1];
	     	value[0][0] = new Double(250000000);
	     	innerMap.put("Regular_Volatility_Currencies", value); // or currency?
	     	Map<String,Map<String, Double[][]>> secondMap=new HashMap<String,Map<String, Double[][]>>();
	     	secondMap.put("InterestRate", innerMap);
	     	this.MapRiskClassThresholdMap.put("delta", secondMap);
	     	
	     	// Set RiskClassRiskWeightMap
	     	innerMap= new HashMap<String, Double[][]>(); Double[][] value1 = new Double[1][riskWeightsRegularCurrency.length];
	     	value1[0] = ArrayUtils.toObject(riskWeightsRegularCurrency);
	     	innerMap.put("Regular_Volatility_Currencies", value1); // or currency?
	     	value1 = new Double[1][riskWeightsRegularCurrency.length];
	     	value1[0] = ArrayUtils.toObject(riskWeightsLowVolCurrency);
	     	innerMap.put("Low_Volatility_Currencies", value1); // or currency?
	     	value1 = new Double[1][riskWeightsRegularCurrency.length];
	     	value1[0] = ArrayUtils.toObject(riskWeightsHighVolCurrency);
	     	innerMap.put("High_Volatility_Currencies", value1); // or currency?
	     	value1 = new Double[1][1];
	     	value1[0] = ArrayUtils.toObject(new double[] {0.0032});
	     	innerMap.put("inflation", value1); // or currency?
	     	value1 = new Double[1][1];
	     	value1[0] = ArrayUtils.toObject(new double[] {0.0018});
	     	innerMap.put("ccybasis", value1); // or currency?
	     	secondMap=new HashMap<String,Map<String, Double[][]>>();
	     	secondMap.put("InterestRate", innerMap);
	     	this.MapRiskClassRiskweightMap.put("delta", secondMap);
	     	
	     }
	     
	    private double[] riskWeightsRegularCurrency = new double[]{0.0077, 0.0077, 0.0077, 0.0064, 0.0058, 0.0049, 0.0047, 0.0047,	0.0045,	0.0045,	0.0048,	0.0056};
	     public void setRiskWeightsRegular(double[] weights){
	    	 riskWeightsRegularCurrency = weights;
	    	 Double[][] value = new Double[1][riskWeightsRegularCurrency.length];
	    	 value[0] = ArrayUtils.toObject(riskWeightsRegularCurrency);
	    	 MapRiskClassRiskweightMap.get("delta").get("InterestRate").put("Regular_Volatility_Currencies", value);
	     }
	    private final double[] riskWeightsLowVolCurrency  = new double[]{0.0010, 0.0010, 0.0010, 0.0010, 0.0013, 0.0016, 0.0018, 0.0020, 0.0025, 0.0022, 0.0022, 0.0023};
	 	private final double[] riskWeightsHighVolCurrency = new double[]{0.0089, 0.0089, 0.0089, 0.0094, 0.0104, 0.0099, 0.0096, 0.0099, 0.0087, 0.0097, 0.0097, 0.0098};

	     
	     private final Double[][] correlationMatrixWithinSubCurve = new Double[][]{
				{1.0    , 1.0    , 1.0    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
				{1.0    , 1.0    , 1.0    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
				{1.0    , 1.0    , 1.0    , 0.782, 0.618, 0.498, 0.438, 0.361, 0.27 , 0.196, 0.174, 0.129},
				{0.782, 0.782, 0.782, 1.0    , 0.84 , 0.739, 0.667, 0.569, 0.444, 0.375, 0.349, 0.296},
				{0.618, 0.618, 0.618, 0.84 , 1.0    , 0.917, 0.859, 0.757, 0.626, 0.555, 0.526, 0.471},
				{0.498, 0.498, 0.498, 0.739, 0.917, 1.0    , 0.976, 0.895, 0.749, 0.69 , 0.66 , 0.602},
				{0.438, 0.438, 0.438, 0.667, 0.859, 0.976, 1.0    , 0.958, 0.831, 0.779, 0.746, 0.69 },
				{0.361, 0.361, 0.361, 0.569, 0.757, 0.895, 0.958, 1.0    , 0.925, 0.893, 0.859, 0.812},
				{0.27 , 0.27 , 0.27 , 0.444, 0.626, 0.749, 0.831, 0.925, 1.0    , 0.98 , 0.961, 0.931},
				{0.196, 0.196, 0.196, 0.375, 0.555, 0.69 , 0.779, 0.893, 0.98 , 1.0    , 0.989, 0.97 },
				{0.174, 0.174, 0.174, 0.349, 0.526, 0.66 , 0.746, 0.859, 0.961, 0.989, 1.0    , 0.988},
				{0.129, 0.129, 0.129, 0.296, 0.471, 0.602, 0.69 , 0.812, 0.931, 0.97 , 0.988, 1.0    }
		    };
		
	     
	     public Double[][]               CrossRiskClassCorrelationMatrix;
	     
	     public Map<String,Double>       MapHistoricalVolaRatio;

	     public Map<String,String>       MapFXCategory;

	     final public String[]           ProductClassKeys = {"RatesFX","Credit","Equity","Commodity"};
	     final public String[]           RiskClassKeys = {"InterestRate","CreditQ","CreditNonQ","Equity","Commodity","FX"};
	     final public String[]           CreditMaturityBuckets = {"1y","2y","3y","5y","10y"};
	     final public String[]           IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};
	     final public String[]           IRCurveIndexNames = {"OIS","Libor1m","Libor3m","Libor6m","Libor12m"};

	     public Double             IRCorrelationCrossCurrency = .27;

	     public Map<String,String>       IRCurrencyMap;

	     public Map<String,Double[][]> MapRiskClassCorrelationIntraBucketMap = new HashMap<String,Double[][]>(); // cross tenor for DeltaIM
	     
	     public Map<String,Double[][] > MapRiskClassCorrelationCrossBucketMap = new HashMap<String,Double[][] >();
	     public Map<String,Map<String,Map<String,Double[][]> > > MapRiskClassThresholdMap = new HashMap<String,Map<String,Map<String,Double[][]> > >(); // risk class, currency, threshold
	     public Map<String,Map<String,Map<String,Double[][]> > > MapRiskClassRiskweightMap = new HashMap<String,Map<String,Map<String,Double[][]> > >();


	     public void setIRCurrencyMap(Map<String, String> IRCurrencyMap) { // currency and volatility class
	         this.IRCurrencyMap = IRCurrencyMap;
	     }

	     public ParameterCollection setCrossRiskClassCorrelationMatrix(Double[][] crossRiskClassCorrelationMatrix) {
	         CrossRiskClassCorrelationMatrix = crossRiskClassCorrelationMatrix;
	         return this;
	     }

	 }
	
	
	/*
	 *  Commencement of actual SIMM scheme 
	 */

    public Map<String,Double>   resultMap;

    private AbstractSIMMProduct[] products;
    
    private ParameterCollection parameterCollection;
    private String[] productClassKeys;
    private String[] riskClassKeys;
    private String[] IRCurveIndexNames;
    private String calculationCCY;

   
    // SIMM constructor
    public CalculationSchemeInitialMarginISDA(SIMMPortfolio portfolio, 
    					  //ParameterCollection parameterCollection, /* Uncomment this line if parameter collection constructor does not contain hard values */
    					  String calculationCCY) throws CalculationException{
        this.resultMap = new HashMap<>();
        this.calculationCCY = calculationCCY;
        this.products = portfolio.getProducts();
        this.parameterCollection = new ParameterCollection();
        
        // Inserted by Mario Viehmann: Screen portfolio products for relevant product classes, risk classes and curveIndexNames
        ArrayList<String> relevantProductClasses = new ArrayList<String>();
        ArrayList<String> relevantCurveIndices = new ArrayList<String>();
        ArrayList<String> relevantRiskClasses = new ArrayList<String>();
       
        for(int productIndex=0;productIndex<products.length;productIndex++){
        	String productClass = products[productIndex].getProductClass();
            if(!relevantProductClasses.contains(productClass)) relevantProductClasses.add(productClass);
            String[] curveIndex = products[productIndex].getCurveIndexNames();
            for(int i=0;i<curveIndex.length;i++){
            	if(!relevantCurveIndices.contains(curveIndex[i])) relevantCurveIndices.add(curveIndex[i]);
            }
            String[] riskClass = products[productIndex].getRiskClasses();
            for(int i=0;i<riskClass.length;i++){
            	if(!relevantRiskClasses.contains(riskClass[i])) relevantRiskClasses.add(riskClass[i]);
            }
        }
        this.productClassKeys = relevantProductClasses.toArray(new String[relevantProductClasses.size()]);
        this.riskClassKeys   = relevantRiskClasses.toArray(new String[relevantRiskClasses.size()]);
        this.IRCurveIndexNames = relevantCurveIndices.toArray(new String[relevantCurveIndices.size()]);
          
    }
    
    // SIMM constructor
    public CalculationSchemeInitialMarginISDA(AbstractSIMMProduct product, String calculationCCY) throws CalculationException{
        this.resultMap = new HashMap<>();
        this.calculationCCY = calculationCCY;
        this.parameterCollection = new ParameterCollection();
        this.products = new AbstractSIMMProduct[]{product};
        this.productClassKeys = new String[]{product.getProductClass()};
        this.riskClassKeys = product.getRiskClasses();
        this.IRCurveIndexNames = product.getCurveIndexNames();
    }

    public RandomVariableInterface getValue(double evaluationTime) throws CalculationException{
    	RandomVariableInterface SIMMValue = null;
        for ( String productClass : productClassKeys) { // RatesFX, Credit etc.
            RandomVariableInterface SIMMProductValue = this.getSIMMProduct(productClass,evaluationTime);
            SIMMValue = SIMMValue == null ? SIMMProductValue : SIMMValue.add(SIMMProductValue);
        }
        return SIMMValue;
    }
    
    
    // SIMM constructor 
    public CalculationSchemeInitialMarginISDA(String calculationCCY) throws CalculationException{
        this.resultMap = new HashMap<>();
        this.calculationCCY = calculationCCY;
        this.parameterCollection = new ParameterCollection();
    }
    
    
    public RandomVariableInterface getValue(AbstractSIMMProduct product, double evaluationTime) throws CalculationException{
    	RandomVariableInterface SIMMValue = null;
    	this.products = new AbstractSIMMProduct[]{product};
        this.productClassKeys = new String[]{product.getProductClass()};
        this.riskClassKeys = product.getRiskClasses();
        this.IRCurveIndexNames = product.getCurveIndexNames();
        for ( String productClass : productClassKeys) { // RatesFX, Credit etc.
            RandomVariableInterface SIMMProductValue = this.getSIMMProduct(productClass,evaluationTime);
            SIMMValue = SIMMValue == null ? SIMMProductValue : SIMMValue.add(SIMMProductValue);
        }
        return SIMMValue;
    }
    
    public void setRiskWeightsRegular(double[] weights){
    	this.parameterCollection.setRiskWeightsRegular(weights);
    }
    
    
    // for all times!
    public Map<String,Double> getResultMap() throws CalculationException{
        if (this.resultMap.entrySet().size()== 0)
            this.getValue(0.0);
        return this.resultMap;
    }

    // returns IM for productClass e.g. RatesFX
    public RandomVariableInterface      getSIMMProduct(String productClass, double atTime){

        Set<String> riskClassList = Stream.of(riskClassKeys).collect(Collectors.toSet());

        // contributions of risk classes to IM within that product class
        RandomVariableInterface[] contributions = new RandomVariableInterface[riskClassKeys.length];
        int i = 0;
        //this.portfolioProducts[0].get//this.nettingset.getActiveProduct(this.nettingset.getActiveProductKeys()[0]).getSensitivitySet(atTime).getPathDimension();
        for ( String iRiskClass : riskClassKeys)
        {
            if ( riskClassList.contains(iRiskClass)) { // riskClassList == iRiskClass?
                RandomVariableInterface iIM = this.getIMForRiskClass(iRiskClass, productClass, atTime);
                contributions[i] = iIM;
            }
            else
                contributions[i] = new RandomVariable(atTime,0.0);
            i++;
        }

        RandomVariableInterface simmProductClass = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(contributions,parameterCollection.CrossRiskClassCorrelationMatrix);
        resultMap.put(productClass,simmProductClass.getAverage());
        return simmProductClass;
    }


//    public RandomVariableInterface getIMForRiskClass(String riskClassKey,String productClass, double atTime){
//        RandomVariableInterface    deltaMargin = this.getDeltaMargin(riskClassKey,productClass,atTime);
//        //RandomVariableInterface    vegaMargin = this.getVegaMargin(riskClassKey,productClass,atTime);
//        //RandomVariableInterface    curatureMargin = this.getDeltaMargin(riskClassKey,productClass, atTime);
//        resultMap.put(productClass+"-"+riskClassKey+"-DeltaMargin",deltaMargin.getAverage());
//        return deltaMargin; //.add(vegaMargin);//.add(vegaMargin).add(curatureMargin);
//    }

    
//    public RandomVariableInterface getDeltaMargin(String riskClassKey,String productClassKey, double atTime){
//        RandomVariableInterface deltaMargin = null;
//        //DeltaMarginSchemeNonIR DeltaScheme = new DeltaMarginSchemeNonIR(this,"Risk_IRCurve",productClassKey,
//        //                                        this.riskClassRiskWeightMap.get(riskClassKey),this.riskClassCorrelationMap.get(riskClassKey),this.riskClassThresholdMap.get(riskClassKey));
//
//        String riskType = "delta";
//        if ( riskClassKey.equals("InterestRate"))
//        {
//            SIMMSchemeIRDelta DeltaScheme = new SIMMSchemeIRDelta(this,productClassKey);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        /*else if ( riskClassKey.equals("CreditQ")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditQ",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("CreditNonQ")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditNonQ",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("Equity")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Equity",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("Commodity")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Commodity",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("FX")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"FX",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }*/
//
//        //MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,riskClassKey,productClassKey,riskType);
//        return deltaMargin;
//    }
    
    
    public RandomVariableInterface      getIMForRiskClass(String riskClassKey,String productClass, double atTime){
        RandomVariableInterface    deltaMargin = this.getDeltaMargin(riskClassKey,productClass,atTime);
        RandomVariableInterface    vegaMargin = this.getVegaMargin(riskClassKey,productClass,atTime);
        //RandomVariableInterface    curatureMargin = this.getDeltaMargin(riskClassKey,productClass, atTime);
        resultMap.put(productClass+"-"+riskClassKey+"-DeltaMargin",deltaMargin.getAverage());
        resultMap.put(productClass+"-"+riskClassKey+"-VegaMargin",vegaMargin.getAverage());

        return deltaMargin.add(vegaMargin);//.add(vegaMargin).add(curatureMargin);

    }

    public RandomVariableInterface      getDeltaMargin(String riskClassKey,String productClassKey, double atTime){
        RandomVariableInterface deltaMargin = null;
        //DeltaMarginSchemeNonIR DeltaScheme = new DeltaMarginSchemeNonIR(this,"Risk_IRCurve",productClassKey,
        //                                        this.riskClassRiskWeightMap.get(riskClassKey),this.riskClassCorrelationMap.get(riskClassKey),this.riskClassThresholdMap.get(riskClassKey));

        String riskType = "delta";
        if ( riskClassKey.equals("InterestRate"))
        {
            MarginSchemeIRDelta DeltaScheme = new MarginSchemeIRDelta(this,productClassKey);
            deltaMargin = DeltaScheme.getValue(atTime);
        }
//        else if ( riskClassKey.equals("CreditQ")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditQ",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("CreditNonQ")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"CreditNonQ",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("Equity")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Equity",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("Commodity")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"Commodity",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }
//        else if ( riskClassKey.equals("FX")){
//            MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,"FX",productClassKey,riskType);
//            deltaMargin = DeltaScheme.getValue(atTime);
//        }

        //MarginSchemeDeltaVega DeltaScheme = new MarginSchemeDeltaVega(this,riskClassKey,productClassKey,riskType);
        return deltaMargin;
    }


    public RandomVariableInterface      getVegaMargin(String riskClassKey,String productClassKey,double atTime){
        // MarginSchemeDeltaVega VegaScheme = new MarginSchemeDeltaVega(this,riskClassKey,productClassKey,"vega");
        return new RandomVariable(0.0);//VegaScheme.getValue(atTime);
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
        
        RandomVariableInterface isdasimmsensiofAllProducts = Stream.of(products).map(
        		product->{try{
        					return product.getSensitivity(productClassKey, 
        							                      riskClassKey, 
                                                          maturityBucket, // only for IR and Credit risk class, null otherwise
                                                          riskFactor,     // CurveIndexName, can be inflation?! 
                                                          bucketKey,      // currency for IR otherwise null ?
                                                          riskType, atTime);
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
    
    // Inserted by Mario Viehmann
    public String[] getInterestRateDeltaBucketKeys(){
           ArrayList<String> relevantBuckets = new ArrayList<String>();
           for(int productIndex=0;productIndex<products.length;productIndex++){
        	  String bucket = products[productIndex].getCurrency();
              if(!relevantBuckets.contains(bucket)) relevantBuckets.add(bucket);
           }
           return relevantBuckets.toArray(new String[relevantBuckets.size()]);
    }
    
    public double[] getRiskWeightsCalibrated(final LIBORModelMonteCarloSimulationInterface model, final SIMMSimpleSwap[] calibrationProducts, final double[] calibrationTargetValues, Map<String,Object> calibrationParameters) throws CalculationException {

		if(calibrationParameters == null) calibrationParameters = new HashMap<String,Object>();
		Integer maxIterationsParameter	= (Integer)calibrationParameters.get("maxIterations");
		Double	parameterStepParameter	= (Double)calibrationParameters.get("parameterStep");
		Double	accuracyParameter		= (Double)calibrationParameters.get("accuracy");
		
		double[] initialParameters = Arrays.stream(this.parameterCollection.riskWeightsRegularCurrency).map(n->Math.log(n)).toArray();
		double[] lowerBound = new double[initialParameters.length];
		double[] upperBound = new double[initialParameters.length];
		double[] parameterStep = new double[initialParameters.length];
		double[] zero = new double[calibrationTargetValues.length];
		Arrays.fill(lowerBound, 0.0);
		Arrays.fill(upperBound, 1.0);
		Arrays.fill(parameterStep, parameterStepParameter != null ? parameterStepParameter.doubleValue() : 1E-5);
		Arrays.fill(zero, 0);

		int numberOfThreads = 2;
		OptimizerFactoryInterface optimizerFactoryParameter = (OptimizerFactoryInterface)calibrationParameters.get("optimizerFactory");

		int maxIterations	= maxIterationsParameter != null ? maxIterationsParameter.intValue() : 2000;
		double accuracy		= accuracyParameter != null ? accuracyParameter.doubleValue() : 1E-5;
		OptimizerFactoryInterface optimizerFactory = optimizerFactoryParameter != null ? optimizerFactoryParameter : new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		//int numberOfThreadsForProductValuation = 2 * Math.max(2, Runtime.getRuntime().availableProcessors());
		final ExecutorService executor = null;//Executors.newFixedThreadPool(numberOfThreadsForProductValuation);

		ObjectiveFunction calibrationError = new ObjectiveFunction() {			
			// Calculate ISDA SIMM IM 
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {

				parameters = Arrays.stream(parameters).map(n->Math.exp(n)).toArray();
				CalculationSchemeInitialMarginISDA.this.setRiskWeightsRegular(parameters);

				ArrayList<Future<Double>> valueFutures = new ArrayList<Future<Double>>(calibrationProducts.length);
				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					final int workerCalibrationProductIndex = calibrationProductIndex;
					Callable<Double> worker = new  Callable<Double>() {
						public Double call() throws SolverException {
							try {
								return calibrationProducts[workerCalibrationProductIndex].getInitialMargin(0.0 /*evaluationTime*/, model, CalculationSchemeInitialMarginISDA.this).getAverage();
							} catch (CalculationException e) {
								// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
								return calibrationTargetValues[workerCalibrationProductIndex];
							} catch (Exception e) {
								// We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
								return calibrationTargetValues[workerCalibrationProductIndex];
							}
						}
					};
					if(executor != null) {
						Future<Double> valueFuture = executor.submit(worker);
						valueFutures.add(calibrationProductIndex, valueFuture);
					}
					else {
						FutureTask<Double> valueFutureTask = new FutureTask<Double>(worker);
						valueFutureTask.run();
						valueFutures.add(calibrationProductIndex, valueFutureTask);
					}
				}
				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					try {
						double value = valueFutures.get(calibrationProductIndex).get();
						values[calibrationProductIndex] = value;
					}
					catch (InterruptedException e) {
						throw new SolverException(e);
					} catch (ExecutionException e) {
						throw new SolverException(e);
					}
				}
			}
		};

		OptimizerInterface optimizer = optimizerFactory.getOptimizer(calibrationError, initialParameters, lowerBound, upperBound, parameterStep, calibrationTargetValues);
		try {
			optimizer.run();
		}
		catch(SolverException e) {
			throw new CalculationException(e);
		}
		finally {
			if(executor != null) {
				executor.shutdown();
			}
		}

		// Get covariance model corresponding to the best parameter set.
		return Arrays.stream(optimizer.getBestFitParameters()).map(n->Math.exp(n)).toArray();
   	
	}
    
    
//    public Map<String,String[]>     getMapRiskClassRiskFactors(String riskTypeString, String bucketKey,double atTime){
////        String[] tradeList = nettingset.getActiveProductKeys();
////        Map<String,SensitivitySet> tradeSensitivityMap = Stream.of(tradeList).collect(Collectors.toMap(t->t,t->nettingset.getActiveProduct(t).getSensitivitySet(atTime)));
////
////        List<String> riskClassKeys = getRiskClassKeys(tradeSensitivityMap);
////        
//        List<String> riskClassKeys = new ArrayList<>(Arrays.asList(this.riskClassKeys));
//        if (!riskTypeString.equals("vega") ) {
//            Map<String, String[]> mapRiskClassRiskFactorKeys = new HashMap<>();
////        if ( riskTypeString.equals("delta")) {
//            riskClassKeys.stream().forEach(riskClass -> {
//                List<String> riskFactors = tradeSensitivityMap.keySet().stream()
//                        .flatMap(
//                                key -> tradeSensitivityMap.get(key).getKeySet().stream().filter(k -> k!=null &&  k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
//                        ).distinct().collect(Collectors.toList());
//                mapRiskClassRiskFactorKeys.put(riskClass, riskFactors.toArray(new String[riskFactors.size()]));
//            });
//            return mapRiskClassRiskFactorKeys;
////        }
//        }
//        else{
//            Map<String,String[]>     mapRiskClassRiskFactorKeys = new HashMap<>();
//            riskClassKeys.stream().forEach(riskClass -> {
//                if ( riskClass.equals("InterestRate") ) {
//                    List<String> riskFactors = tradeSensitivityMap.keySet().stream()
//                            .flatMap(
//                                    key -> tradeSensitivityMap.get(key).getKeySet().stream().filter(k -> k!=null && k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
//                            ).distinct().collect(Collectors.toList());
//                    mapRiskClassRiskFactorKeys.put(riskClass, riskFactors.toArray(new String[riskFactors.size()]));
//                }
//                else{
//                    List<String> riskFactors = tradeSensitivityMap.keySet().stream()
//                            .flatMap(
//                                    key -> tradeSensitivityMap.get(key).getKeySet().stream().filter(k ->  k!=null && k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
//                            ).distinct().collect(Collectors.toList());
//                    mapRiskClassRiskFactorKeys.put(riskClass, riskFactors.toArray(new String[riskFactors.size()]));
//                }
//            });
//            return mapRiskClassRiskFactorKeys;
//        }
//
//    }
}