package initialmargin.isdasimm;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;


import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SIMMSchemeMain  {

	public static class ParameterCollection{
        public ParameterCollection(){

        }
        public Double[][]               CrossRiskClassCorrelationMatrix;

        public Map<String,String>       MapFXCategory;

        final public String[]           ProductClassKeys = {"RatesFX","Credit","Equity","Commodity"};
        final public String[]           RiskClassKeys = {"InterestRate","CreditQ","CreditNonQ","Equity","Commodity","FX"};
        final public String[]           CreditMaturityBuckets = {"1y","2y","3y","5y","10y"};
        final public String[]           IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};
        final public String[]           IRCurveIndexNames = {"OIS","Libor1m","Libor3m","Libor6m","Libor12m"};

        public Double             IRCorrelationCrossCurrency;// = .27;

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

    Map<AbstractMap.SimpleEntry<Double,SIMMSchemeSensitivitySet.Key>,RandomVariableInterface> CacheNetSensitivities;

    public Map<String,Double>   resultMap;

    AbstractLIBORMonteCarloProduct[] portfolioProducts;
    private final LIBORModelMonteCarloSimulationInterface model;

    private ParameterCollection parameterCollection;

    private String calculationCCY;

    public SIMMSchemeMain(AbstractLIBORMonteCarloProduct[] portfolioProducts, 
    						LIBORModelMonteCarloSimulationInterface model, 
    						ParameterCollection parameterCollection, 
    						String calculationCCY){
        this.resultMap = new HashMap<>();
        this.calculationCCY = calculationCCY;
        this.portfolioProducts=portfolioProducts;
        this.model = model;
        CacheNetSensitivities = new HashMap<>();
       

        this.parameterCollection = parameterCollection;
     }



    

    public String    getCalculationCCY(){
        return this.calculationCCY;
    }

    public int getPathDimension(){
        return 1;
    }

    public ParameterCollection getParameterCollection() {
        return parameterCollection;
    }

    public void setParameterCollection(ParameterCollection parameterCollection) {
        this.parameterCollection = parameterCollection;
    }



    public RandomVariableInterface      getValue(double atTime){

        RandomVariableInterface SIMMValue = null;
        for ( String productClass : parameterCollection.ProductClassKeys) {
            RandomVariableInterface SIMMProductValue = this.getSIMMProduct(productClass,atTime);
            SIMMValue = SIMMValue == null ? SIMMProductValue : SIMMValue.add(SIMMProductValue);
        }
        return SIMMValue;
    }

    public Map<String,Double>      getResultMap(){
        if (this.resultMap.entrySet().size()== 0)
            this.getValue(0.0);
        return this.resultMap;
    }

    public RandomVariableInterface      getSIMMProduct(String productClass, double atTime){

        Set<String> riskClassList = Stream.of(this.parameterCollection.RiskClassKeys).collect(Collectors.toSet());

        RandomVariableInterface[] contributions = new RandomVariableInterface[parameterCollection.RiskClassKeys.length];
        int i = 0;
        int pathDim = 1;//this.portfolioProducts[0].get//this.nettingset.getActiveProduct(this.nettingset.getActiveProductKeys()[0]).getSensitivitySet(atTime).getPathDimension();
        for ( String iRiskClass : parameterCollection.RiskClassKeys)
        {
            if ( riskClassList.contains(iRiskClass)) {
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


    public RandomVariableInterface      getIMForRiskClass(String riskClassKey,String productClass, double atTime){
        RandomVariableInterface    deltaMargin = this.getDeltaMargin(riskClassKey,productClass,atTime);
        RandomVariableInterface    vegaMargin = this.getVegaMargin(riskClassKey,productClass,atTime);
        //RandomVariableInterface    curatureMargin = this.getDeltaMargin(riskClassKey,productClass, atTime);
        resultMap.put(productClass+"-"+riskClassKey+"-DeltaMargin",deltaMargin.getAverage());

        return deltaMargin.add(vegaMargin);//.add(vegaMargin).add(curatureMargin);

    }

    public RandomVariableInterface      getDeltaMargin(String riskClassKey,String productClassKey, double atTime){
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


    public     SIMMSchemeSensitivitySet.Key  getKeySensitivitySet(int iTenor,String riskFactorKey, String bucketKey, String productClass, String riskClass,String riskType)
    {
        if ( riskClass.equals("InterestRate") ) {

            String maturityBucket = parameterCollection.IRMaturityBuckets[iTenor];
            //String curveIndex = parameterCollection.IRCurveIndexNames[iRiskFactor];
            if (riskFactorKey.equals("inflation") || riskFactorKey.equals("ccybasis") || riskType.equals("vega"))
                return new SIMMSchemeSensitivitySet.Key("None",riskFactorKey,bucketKey,riskClass,riskType,productClass);
            else
                return new SIMMSchemeSensitivitySet.Key(maturityBucket,riskFactorKey,bucketKey,riskClass,riskType,productClass);
        }
        else if (riskClass.contains("Credit")){
            String maturityBucket = parameterCollection.CreditMaturityBuckets[iTenor];
            //String riskFactorKey = this.parameterCollection.MapRiskClassRiskFactorKeys.get(riskClass)[iRiskFactor];
            return new SIMMSchemeSensitivitySet.Key(maturityBucket,riskFactorKey,bucketKey,riskClass,riskType,productClass);
        }
        else if (riskClass.contains("FX")){ // BucketID is default "0"
            String maturityBucket = parameterCollection.CreditMaturityBuckets[iTenor];
            //String riskFactorKey = this.parameterCollection.MapRiskClassRiskFactorKeys.get(riskClass)[iRiskFactor];
            return new SIMMSchemeSensitivitySet.Key("None",riskFactorKey,"0",riskClass,riskType,productClass);
        }
        else {
            //String riskFactorKey = this.parameterCollection.MapRiskClassRiskFactorKeys.get(riskClass)[iRiskFactor];
            return new SIMMSchemeSensitivitySet.Key("None",riskFactorKey,bucketKey,riskClass,riskType,productClass);
        }

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

    // BUCKET IS CURRENCY FOR IR
    public RandomVariableInterface   getNetSensitivity(String productClassKey, String riskClassKey, int iRateTenor,String riskFactor, String bucketKey,String riskType, double atTime) {

        SIMMSchemeSensitivitySet.Key key = this.getKeySensitivitySet(iRateTenor,riskFactor,bucketKey,productClassKey,riskClassKey,riskType);
        
        RandomVariableInterface isdasimmsensiofAllProducts = Stream.of(this.portfolioProducts).map(
        		product->{try{
        					return product.getValue(atTime,model);
        					}
        					catch(Exception e){
        						return null;
        						}
        		} ).reduce(null, (e1, e2) -> {
            if (e1 == null) return e2;
            if (e2 == null) return e1;
            return e1.add(e2);
        });
        
        return isdasimmsensiofAllProducts;

      

    }

    
    /*
    public static List<String>   getRiskClassKeys(Map<String, SIMMSchemeSensitivitySet> tradeSensitivityMap)
    {
        //Map RIskClass -> RiskFactors
        List<java.lang.String> riskClasses = tradeSensitivityMap.keySet().stream()
                .flatMap(
                        key -> tradeSensitivityMap.get(key).getKeySet().stream().map(k->k.getRiskClass())
                ).distinct().collect(Collectors.toList());
        return riskClasses;
    }

    public Map<String,String[]>     getMapRiskClassBucketKeys(String riskTypeString){
        String[] tradeList = nettingset.getActiveProductKeys();
        Map<String,SIMMSchemeSensitivitySet> tradeSensitivityMap = Stream.of(tradeList).collect(Collectors.toMap(t->t,t->nettingset.getActiveProduct(t).getSensitivitySet(0.0)));

        List<String> riskClassKeys = getRiskClassKeys(tradeSensitivityMap);
        Map<String,String[]>     mapRiskClassBucketKeys = new HashMap<>();
        riskClassKeys.stream().forEach(riskClass ->{
            List<String> riskFactors = tradeSensitivityMap.keySet().stream()
                    .flatMap(
                            key -> tradeSensitivityMap.get(key).getKeySet().stream().filter(k->k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString)).map(k->k.getBucketKey())
                    ).distinct().collect(Collectors.toList());
            mapRiskClassBucketKeys.put(riskClass,riskFactors.toArray(new String[riskFactors.size()]));
        });
        return mapRiskClassBucketKeys;

    }

    public Map<String,String[]>     getMapRiskClassRiskFactors(String riskTypeString, String bucketKey){
        String[] tradeList = nettingset.getActiveProductKeys();
        Map<String,SIMMSchemeSensitivitySet> tradeSensitivityMap = Stream.of(tradeList).collect(Collectors.toMap(t->t,t->nettingset.getActiveProduct(t).getSensitivitySet(0.0)));

        List<String> riskClassKeys = getRiskClassKeys(tradeSensitivityMap);
        Map<String,String[]>     mapRiskClassRiskFactorKeys = new HashMap<>();
        riskClassKeys.stream().forEach(riskClass ->{
            List<String> riskFactors = tradeSensitivityMap.keySet().stream()
                    .flatMap(
                            key -> tradeSensitivityMap.get(key).getKeySet().stream().filter(k->k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k->k.getRiskFactorKey())
                    ).distinct().collect(Collectors.toList());
            mapRiskClassRiskFactorKeys.put(riskClass,riskFactors.toArray(new String[riskFactors.size()]));
        });
        return mapRiskClassRiskFactorKeys;

    }
*/




}




