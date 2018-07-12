package xva.initialmargin;

import java.util.Map;
import java.util.stream.Stream;

public class SIMMParameter {
    public SIMMParameter(){

    }

    final static public String   inflationKey = "inflation";
    final static public String   ccyBasisKey = "ccybasis";

    public enum RiskType {
        Delta,
        Vega,
        Curvature
    }

    public enum ProductClass {
        RatesFX,
        Credit,
        Equity,
        Commodity
    }

    public enum RiskClass {
        InterestRate,
        CreditQ,
        CreditNonQ,
        Equity,
        Commodity,
        FX
    }

    public enum RatesCurveNames{
        OIS,
        Libor1m,
        Libor3m,
        Libor6m,
        Libor12m,
    }


    public static String[]  getProductClassKeys(){
        return Stream.of(ProductClass.values()).map(ProductClass::name).toArray(String[]::new);
    }

    public static String[]  getRiskClassKeys(){
        return Stream.of(RiskClass.values()).map(RiskClass::name).toArray(String[]::new);
    }

    public static String[]  getRateCurveKeys(){
        return Stream.of(RatesCurveNames.values()).map(RatesCurveNames::name).toArray(String[]::new);
    }


    public Double[][]               CrossRiskClassCorrelationMatrix;

    public Map<String,String> MapFXCategory;

    public Map<String,Double>       MapHistoricalVolaRatio;

//    final public String[]           ProductClassKeys = {"RatesFX","Credit","Equity","Commodity"};
//    final public String[]           RiskClassKeys = {"InterestRate","CreditQ","CreditNonQ","Equity","Commodity","FX"};
    final public String[]           CreditMaturityBuckets = {"1y","2y","3y","5y","10y"};
    final public String[]           IRMaturityBuckets = {"2w","1m","3m","6m","1y","2y","3y","5y","10y","15y","20y","30y"};
//    final public String[]           IRCurveIndexNames = {"OIS","Libor1m","Libor3m","Libor6m","Libor12m"};

    //       final public Double             IRCorrelationCrossCurveIndex = 0.982;
    public Double             IRCorrelationCrossCurrency;// = .27;

    public Map<String,String>       IRCurrencyMap;


    public Map<RiskClass,Double[][] > MapRiskClassCorrelationIntraBucketMap;
    public Map<RiskClass,Double[][] > MapRiskClassCorrelationCrossBucketMap;
    public Map<RiskClass,Map<String,Map<String,Double[][]> > > MapRiskClassThresholdMap;
    public Map<RiskClass,Map<String,Map<String,Double[][]> > > MapRiskClassRiskweightMap;


    public void setIRCurrencyMap(Map<String, String> IRCurrencyMap) {
        this.IRCurrencyMap = IRCurrencyMap;
    }


    public SIMMParameter setCrossRiskClassCorrelationMatrix(Double[][] crossRiskClassCorrelationMatrix) {
        CrossRiskClassCorrelationMatrix = crossRiskClassCorrelationMatrix;
        return this;
    }

}
