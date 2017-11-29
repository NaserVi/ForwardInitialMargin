package initialmargin.isdasimm.aggregationscheme;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


/**
 * Created by xn04210 on 08/08/2017.
 */

/*
IR:
 Currency --> Bucket
 CurveTenor --> RiskFactor


Todo:
    * Residual Bucket Implementierung
    * Correcte RIskWeight Auswahl wenn None/Residual
    * FX and Equity Thresholds
 */

public class MarginSchemeIRDelta {
    CalculationSchemeInitialMarginISDA calculationSchemeInitialMarginISDA;
    String      productClassKey;
    String      riskClassKey;
    String[]    bucketKeys;
    final String riskTypeKey = "delta";

    public MarginSchemeIRDelta(CalculationSchemeInitialMarginISDA calculationSchemeInitialMarginISDA,
                               String productClassKey){
        this.calculationSchemeInitialMarginISDA = calculationSchemeInitialMarginISDA;
        this.riskClassKey = "InterestRate";
        this.productClassKey = productClassKey;
        // Only modification:
        this.bucketKeys = calculationSchemeInitialMarginISDA.getInterestRateDeltaBucketKeys();
    }

    public RandomVariableInterface getValue(double atTime){

        if (this.bucketKeys.length==0)
            return new RandomVariable(atTime,this.calculationSchemeInitialMarginISDA.getPathDimension(),0.0);

        RandomVariableInterface[] S1Contributions = new RandomVariableInterface[this.bucketKeys.length];
        RandomVariableInterface[] KContributions = new RandomVariableInterface[this.bucketKeys.length];
        int i=0;

        RandomVariableInterface[] concentrationFactors = new RandomVariableInterface[this.bucketKeys.length];
        for (String bucketKey : this.bucketKeys)
        {
            RandomVariableInterface[][] netSensitivities = this.getNetSensitivities(bucketKey,atTime);
            concentrationFactors[i] = getConcentrationRiskFactor(bucketKey,netSensitivities,atTime);
            RandomVariableInterface K1 = this.getAggregatedSensitivityForBucket(bucketKey,netSensitivities,concentrationFactors[i],atTime);
            RandomVariableInterface S1 = this.getFactorS(bucketKey,K1,netSensitivities,concentrationFactors[i],atTime);
            S1Contributions[i] = S1;
            KContributions[i] = K1;
            i++;
        }

        RandomVariableInterface deltaMargin = null;
        RandomVariableInterface VarCovar = null;
        Double[][] correlationMatrix = null;

        double singleCorrelation = calculationSchemeInitialMarginISDA.getParameterCollection().IRCorrelationCrossCurrency;
        correlationMatrix = new Double[this.bucketKeys.length][this.bucketKeys.length];
        for (i = 0; i< bucketKeys.length;i++)
        for (int j = 0; j< bucketKeys.length;j++)
            if ( i!=j)
                correlationMatrix[i][j] = getParameterG(concentrationFactors[i],concentrationFactors[j]).getAverage()*singleCorrelation;
        VarCovar = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(S1Contributions, correlationMatrix);


        /*Adjustment on Diagonal*/
        VarCovar = VarCovar.squared();
        RandomVariableInterface SSumSQ = null;
        RandomVariableInterface KSumSQ = null;
        for ( int k = 0;k<S1Contributions.length;k++){
            SSumSQ = SSumSQ == null ? SSumSQ = S1Contributions[k].squared() : SSumSQ.add(S1Contributions[k].squared());
            KSumSQ = KSumSQ == null ? KSumSQ = KContributions[k].squared() : KSumSQ.add(KContributions[k].squared());
        }
        VarCovar = VarCovar.sub(SSumSQ).add(KSumSQ);
        deltaMargin = VarCovar.sqrt();


        return deltaMargin;
    }


    private RandomVariableInterface[][] getNetSensitivities(String bucketKey, double atTime){
        int nTenors = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length;
        int nCurves = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames.length; //calculationSchemeInitialMarginISDA.getIRCurveIndexNames().length;
        RandomVariableInterface[][] netSensitivities = new RandomVariableInterface[nCurves][nTenors];
        
        for (int iCurve = 0; iCurve <nCurves; iCurve++)
        {
            String curveKey = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames[iCurve];
            
                for (int iTenor = 0; iTenor < nTenors; iTenor++) {
                    String maturityBucketKey = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[iTenor];               
                    netSensitivities[iCurve][iTenor] = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey, this.riskClassKey, maturityBucketKey, curveKey, bucketKey, "delta", atTime);
                }
            
        }
        return netSensitivities;
    }

    private RandomVariableInterface getAggregatedSensitivityForBucket(String bucketKey, RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double atTime){
        RandomVariableInterface aggregatedSensi = null;

        int nTenors = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length;
        int nCurves = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames.length;

        int dimensionTotal=nTenors*nCurves+2;
        RandomVariableInterface[] contributions = new RandomVariableInterface[dimensionTotal];

        for (int iCurve = 0; iCurve <nCurves; iCurve++)
            for (int iTenor = 0; iTenor <nTenors; iTenor++)
            {
                String curveKey = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames[iCurve];
                RandomVariableInterface iBucketSensi = this.getWeightedNetSensitivity(iTenor, iCurve,curveKey, bucketKey,netSensitivities,concentrationRiskFactor, atTime);
                contributions[iCurve*nTenors+iTenor] = iBucketSensi;
            }

        RandomVariableInterface inflationSensi = this.getWeightedNetSensitivity(0,0,"inflation",bucketKey,netSensitivities,concentrationRiskFactor,atTime);
        RandomVariableInterface ccyBasisSensi = this.getWeightedNetSensitivity(0,0,"ccybasis",bucketKey,netSensitivities,concentrationRiskFactor,atTime);
        contributions[dimensionTotal-2] = inflationSensi;
        contributions[dimensionTotal-1] = ccyBasisSensi;

        Double[][] crossTenorCorrelation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(riskClassKey);

        aggregatedSensi = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(contributions, crossTenorCorrelation);


        return aggregatedSensi;
    }



    private RandomVariableInterface   getWeightedNetSensitivity(int iRateTenor,int iIndex, String indexName,  String bucketKey, RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double atTime)
    {
        double riskWeight = 0;

        if (!indexName.equals("inflation") && !indexName.equals("ccybasis"))
        {
            Optional<Map.Entry<String, String>> optional = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurrencyMap.entrySet().stream().filter(entry -> entry.getKey().contains(bucketKey)).findAny();
            String currencyMapKey;
            if (!optional.isPresent())
                currencyMapKey = "High_Volatility_Currencies";
            else
                currencyMapKey = optional.get().getValue();
            currencyMapKey = currencyMapKey.replace("_Traded", "").replace("_Well", "").replace("_Less", "");
            Double[] riskWeights = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(riskTypeKey).get("InterestRate").get(currencyMapKey)[0];
            riskWeight = riskWeights[iRateTenor];
            RandomVariableInterface netSensi =  netSensitivities[iIndex][iRateTenor];
            if (netSensi!=null) {
                return netSensi.mult(riskWeight).mult(concentrationRiskFactor);
            }
            else
                return new RandomVariable(atTime,this.calculationSchemeInitialMarginISDA.getPathDimension(),0.0);
        }
        else { /* Inflation or CCYBasis*/
            riskWeight = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(riskTypeKey).get("InterestRate").get(indexName)[0][0];
            String maturityBucket = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[iRateTenor];
            RandomVariableInterface netSensi =  calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,maturityBucket, indexName, bucketKey,this.riskTypeKey, atTime);
            if (netSensi!=null)
            {
                netSensi = netSensi.mult(riskWeight);
                if ( !indexName.equals("ccybasis")) {
                    netSensi = netSensi.mult(concentrationRiskFactor);
                }
                return  netSensi;
            }
            else
                return null;
        }
    }

    public RandomVariableInterface     getParameterG(RandomVariableInterface CR1, RandomVariableInterface CR2){
        RandomVariableInterface min = CR1.barrier(CR1.sub(CR2),CR2,CR1);
        RandomVariableInterface max = CR1.barrier(CR1.sub(CR2),CR1,CR2);
        return min.div(max);

    }



    public RandomVariableInterface getFactorS(String bucketKey, RandomVariableInterface K,RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double atTime){
        RandomVariableInterface sum = this.getWeightedSensitivitySum(bucketKey,netSensitivities,concentrationRiskFactor, atTime);
        RandomVariableInterface S1 = K.barrier(sum.sub(K),K,sum);
        RandomVariableInterface KNegative = K.mult(-1);
        S1 = S1.barrier(S1.sub(KNegative),S1,KNegative);
        return S1;
    }

    private RandomVariableInterface   getWeightedSensitivitySum(String bucketKey,RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double atTime){
        RandomVariableInterface aggregatedSensi = new RandomVariable(atTime,0.0);

        for (int iIndex = 0; iIndex < calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames.length; iIndex++) {
            for (int iTenor = 0; iTenor < calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length; iTenor++) {
                String key = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames[iIndex];
                RandomVariableInterface summand = getWeightedNetSensitivity(iTenor,iIndex, key, bucketKey,netSensitivities,concentrationRiskFactor, atTime);
                if (summand !=null)
                    aggregatedSensi = aggregatedSensi == null ? aggregatedSensi = summand : aggregatedSensi.add(summand);
            }
        }
        RandomVariableInterface inflationSensi = this.getWeightedNetSensitivity(0,0,"inflation",bucketKey,netSensitivities,concentrationRiskFactor,atTime);
        RandomVariableInterface ccyBasisSensi = this.getWeightedNetSensitivity(0,0,"ccybasis",bucketKey,netSensitivities,concentrationRiskFactor,atTime);
        return aggregatedSensi.add(inflationSensi.add(ccyBasisSensi));
    }


    public RandomVariableInterface getConcentrationRiskFactor(String bucketKey, RandomVariableInterface[][] netSensitivities, double atTime){
        RandomVariableInterface sensitivitySum = new RandomVariable(atTime,0.0);
        for (int iIndex = 0; iIndex < calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames.length; iIndex++) {
            for (int iTenor = 0; iTenor < calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length; iTenor++) {
                String key = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurveIndexNames[iIndex];
                RandomVariableInterface summand = netSensitivities[iIndex][iTenor];//calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,iTenor,key,bucketKey,"delta",atTime);//"ccybasis",bucketKey,"delta",atTime);//getWeightedNetSensitivity(iTenor, key, bucketKey, atTime);
                if (summand !=null)
                    sensitivitySum = sensitivitySum == null ? sensitivitySum = summand : sensitivitySum.add(summand);
            }
        }
        RandomVariableInterface inflationSensi = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,"","inflation",bucketKey,"delta",atTime);
        if ( sensitivitySum !=null && inflationSensi !=null)
            sensitivitySum = sensitivitySum.add(inflationSensi); // Inflation Sensi are included in Sum, CCYBasis not

        Optional<Map.Entry<String,String> > optional = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurrencyMap.entrySet().stream().filter(entry->entry.getKey().contains(bucketKey)).findAny();
        String currencyMapKey;
        if (!optional.isPresent())
            currencyMapKey="High_Volatility_Currencies";
        else
            currencyMapKey = optional.get().getValue();

        double concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(this.riskTypeKey).get(riskClassKey).get(currencyMapKey)[0][0];
        RandomVariableInterface CR = (sensitivitySum.abs().div(concentrationThreshold)).sqrt();
        CR = CR.barrier(CR.sub(1.0), CR, 1.0);
        return CR;
    }



}
