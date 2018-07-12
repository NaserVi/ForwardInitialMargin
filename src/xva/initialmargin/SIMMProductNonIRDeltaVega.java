package xva.initialmargin;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;
import xva.initialmargin.simmsensitivityproviders.SIMMSensitivityProviderInterface;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SIMMProductNonIRDeltaVega {
    String                  calculationCCY;
    String                  productClassKey;
    SIMMParameter.RiskClass riskClassKey;
    String[]                activeBucketKeys;
    SIMMParameter.RiskType  riskTypeKey;
    final SIMMHelper        helper;
    final SIMMParameter     parameterSet;
    private SIMMSensitivityProviderInterface simmSensitivitivityProvider;

    public SIMMProductNonIRDeltaVega(SIMMSensitivityProviderInterface simmSensitivitivityProvider,
                                     String riskClassKey,
                                     String productClassKey,
                                     String riskTypeKey, SIMMParameter parameterSet, String calculationCCY, double atTime){
        this.calculationCCY = calculationCCY;
        this.helper = new SIMMHelper(simmSensitivitivityProvider.getSIMMTradeSpecs());
        this.parameterSet = parameterSet;
        this.simmSensitivitivityProvider = simmSensitivitivityProvider;
        this.riskClassKey = SIMMParameter.RiskClass.valueOf(riskClassKey);
        this.productClassKey = productClassKey;
        this.riskTypeKey = SIMMParameter.RiskType.valueOf(riskTypeKey);
        this.activeBucketKeys = helper.getRiskClassBucketKeyMap(riskTypeKey,atTime).get(riskClassKey).stream().filter(e->!e.equals("Residual")).toArray(String[]::new);
    }

    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){


        RandomVariableInterface deltaMargin = new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0);
        Double[][] correlationMatrix = parameterSet.MapRiskClassCorrelationCrossBucketMap.get(this.riskClassKey);

        int length = correlationMatrix.length == 1 ? this.activeBucketKeys.length : correlationMatrix.length;

        if(this.activeBucketKeys.length>0) {
            RandomVariableInterface[] KContributions = new RandomVariableInterface[length];
            RandomVariableInterface[] S1Contributions = new RandomVariableInterface[length];
            for (int iBucket = 0; iBucket < activeBucketKeys.length; iBucket++)
            {
                String bucketKey =null;
                int bucketIndex = 0;
                if (riskClassKey.equals(SIMMParameter.RiskClass.Equity) || riskClassKey.equals(SIMMParameter.RiskClass.Commodity)) {
                    bucketKey = new Integer(activeBucketKeys[iBucket]).toString();
                    bucketIndex = Integer.parseInt(activeBucketKeys[iBucket]);
                }
                else {
                    bucketKey = this.activeBucketKeys[iBucket];
                    bucketIndex = iBucket;
                }

                /*Check whether we have risk factors in that bucket*/
                Set<String> activeRiskFactorKeys = helper.getRiskClassRiskFactorMap(this.riskTypeKey.name(),bucketKey,evaluationTime).get(riskClassKey);
                if (activeRiskFactorKeys != null && activeRiskFactorKeys.size()>0) {
                    Map<String,RandomVariableInterface>  netSensitivityMap = this.getRiskFactorNetSensitivityMap(bucketKey, activeRiskFactorKeys, evaluationTime,model);
                    RandomVariableInterface K1 = getAggregatedSensitivityForBucket(bucketKey, netSensitivityMap, evaluationTime);
                    RandomVariableInterface sumWeigthedNetSensi = this.getRiskFactorWeightedNetSensitivityMap(bucketKey, netSensitivityMap, evaluationTime).values().stream().reduce((rv1,rv2)->rv1.add(rv2)).orElseGet(()->new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0));//this.getWeightedSensitivitySum(bucketKey, weightedNetSensitivities, evaluationTime);
                    RandomVariableInterface S1 = K1.barrier(sumWeigthedNetSensi.sub(K1), K1, sumWeigthedNetSensi);
                    RandomVariableInterface KNegative = K1.mult(-1);
                    S1 = S1.barrier(S1.sub(KNegative), S1, KNegative);
                    S1Contributions[bucketIndex] = S1;
                    KContributions[bucketIndex] = K1;
                }
            }

            RandomVariableInterface VarCovar = helper.getVarianceCovarianceAggregation(S1Contributions, correlationMatrix);

            if (VarCovar == null)
                return deltaMargin;
            else {
            /*Adjustment on Diagonal*/
                VarCovar = VarCovar.squared();
                RandomVariableInterface SSumSQ = null;
                RandomVariableInterface KSumSQ = null;
                for (int k = 0; k < S1Contributions.length; k++) {
                    if (S1Contributions[k] != null) {
                        SSumSQ = SSumSQ == null ? SSumSQ = S1Contributions[k].squared() : SSumSQ.add(S1Contributions[k].squared());
                        KSumSQ = KSumSQ == null ? KSumSQ = KContributions[k].squared() : KSumSQ.add(KContributions[k].squared());
                    }
                }
                VarCovar = VarCovar.sub(SSumSQ).add(KSumSQ);
                deltaMargin = VarCovar.sqrt();
            }
        }

        /* RESIDUAL TERM*/
        if (!this.riskClassKey.equals(SIMMParameter.RiskClass.FX)) {
            String bucketKey = "Residual";
            Set<String> activeRiskFactorKeys = this.helper.getRiskClassRiskFactorMap(this.riskTypeKey.name(),bucketKey,evaluationTime).get(riskClassKey);
            if (activeRiskFactorKeys != null && activeRiskFactorKeys.size()>0) {
                Map<String,RandomVariableInterface>  netSensitivityMap = this.getRiskFactorNetSensitivityMap(bucketKey, activeRiskFactorKeys, evaluationTime,model);
                Map<String,RandomVariableInterface> weightedNetSensitivityMap = this.getRiskFactorWeightedNetSensitivityMap(bucketKey, netSensitivityMap, evaluationTime);
                RandomVariableInterface KResidual = getAggregatedSensitivityForBucket(bucketKey, weightedNetSensitivityMap, evaluationTime);
                deltaMargin = deltaMargin.add(KResidual);
            }
        }

        return deltaMargin;
    }


    private Map<String,RandomVariableInterface>   getRiskFactorNetSensitivityMap(String bucketKey,Set<String> activeRiskFactorKeys, double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
        return activeRiskFactorKeys.stream().collect(Collectors.toMap(activeRiskFactorKey->activeRiskFactorKey,activeRiskFactor->this.getNetSensitivity(activeRiskFactor,bucketKey,evaluationTime,model) ));
    }

    private Map<String,RandomVariableInterface>  getConcentrationFactorMap(String bucketKey,Map<String,RandomVariableInterface> riskFactorNetSensitivityMap, double atTime){

        return
                riskFactorNetSensitivityMap.entrySet().stream().collect(Collectors.toMap(entry->entry.getKey(),entry->{
                        return  this.getConcentrationRiskFactor(entry.getValue(),entry.getKey(),bucketKey,atTime); }));
    }


    private Map<String,RandomVariableInterface>  getRiskFactorWeightedNetSensitivityMap(String bucketKey,Map<String,RandomVariableInterface> riskFactorNetSensitivityMap, double evaluationTime){

        Map<String,RandomVariableInterface> concentrationFactors = this.getConcentrationFactorMap(bucketKey,riskFactorNetSensitivityMap,evaluationTime);


        Map<String,RandomVariableInterface> weightedSensiMap =riskFactorNetSensitivityMap.entrySet().stream().collect(Collectors.toMap(entry->entry.getKey(),entry->{
            RandomVariableInterface netSensi  = entry.getValue();
            RandomVariableInterface concentrationRiskFactor = concentrationFactors.get(entry.getKey());
            return this.getWeightedNetSensitivity(netSensi,concentrationRiskFactor,entry.getKey(),bucketKey, evaluationTime);
        }));

        return weightedSensiMap;
    }


    private RandomVariableInterface getAggregatedSensitivityForBucket(String bucketKey, Map<String,RandomVariableInterface>  netSensitivityMap, double evaluationTime){
        RandomVariableInterface aggregatedSensi = null;
        Double[][] correlationMatrix = new Double[ netSensitivityMap.size()][netSensitivityMap.size()];
        Map<String,RandomVariableInterface> weightedNetSensitivityMap = this.getRiskFactorWeightedNetSensitivityMap(bucketKey, netSensitivityMap, evaluationTime);
        String[] activeRiskFactorKeys = netSensitivityMap.keySet().stream().toArray(String[]::new);
        if (riskClassKey.equals(SIMMParameter.RiskClass.InterestRate) && riskTypeKey.equals(SIMMParameter.RiskType.Vega)) {
            correlationMatrix = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get("InterestRate_Tenor");
            RandomVariableInterface[] contributionsReDim = new RandomVariableInterface[correlationMatrix.length];
            int nTenors = this.parameterSet.IRMaturityBuckets.length;
            for (int iIndex = 0; iIndex < activeRiskFactorKeys.length; iIndex++)
                for (int i = 0; i < nTenors; i++) {
                    if (activeRiskFactorKeys[iIndex].equals(this.parameterSet.IRMaturityBuckets[i]))
                        contributionsReDim[i] = weightedNetSensitivityMap.get(activeRiskFactorKeys[iIndex]);
                    if (activeRiskFactorKeys[iIndex].equals(SIMMParameter.inflationKey))
                        contributionsReDim[contributionsReDim.length-2]=weightedNetSensitivityMap.get(activeRiskFactorKeys[iIndex]);
                }
            aggregatedSensi = helper.getVarianceCovarianceAggregation(contributionsReDim,correlationMatrix);
        }
        else {
            Double correlation = 0.0;
            if (riskClassKey.equals(SIMMParameter.RiskClass.FX))
                correlation = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][0];
            else if (riskClassKey.equals(SIMMParameter.RiskClass.CreditQ) || riskClassKey.equals(SIMMParameter.RiskClass.CreditNonQ))
                correlation = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][1];
            else {
                int bucketNr = 0;
                try {
                    bucketNr = (int) Double.parseDouble(bucketKey);
                    correlation = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][bucketNr];
                } catch (Exception e) {
                    bucketNr = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0].length - 1;
                    correlation = this.parameterSet.MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][bucketNr];
                }
            }
            Map<String,RandomVariableInterface> concentrationFactors = this.getConcentrationFactorMap(bucketKey,netSensitivityMap,evaluationTime);
            RandomVariableInterface[] weightedNetSensitivitesArray = new RandomVariableInterface[activeRiskFactorKeys.length];
            for (int i = 0; i < activeRiskFactorKeys.length; i++) {
                String iRiskFactorKey = activeRiskFactorKeys[i];
                weightedNetSensitivitesArray[i] = netSensitivityMap.get(iRiskFactorKey);
                for (int j = 0; j < activeRiskFactorKeys.length; j++)
                    if (i != j) {
                        String jRiskFactorKey = activeRiskFactorKeys[j];
                        correlationMatrix[i][j] = getParameterF(concentrationFactors.get(iRiskFactorKey), concentrationFactors.get(jRiskFactorKey)).getAverage() * correlation;
                    }

            }
            aggregatedSensi = helper.getVarianceCovarianceAggregation(weightedNetSensitivitesArray,correlationMatrix);

        }

        return aggregatedSensi;
    }


    private RandomVariableInterface   getWeightedNetSensitivity(RandomVariableInterface netSensi, RandomVariableInterface concentrationRiskFactor, String riskFactorKey ,String bucketKey, double atTime)
    {
        int bucketIndex = 0;
        try{
            bucketIndex=(int) Double.parseDouble(bucketKey);

        }
        catch(NumberFormatException e){
            bucketIndex= this.parameterSet.MapRiskClassRiskweightMap.get(this.riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0].length-1;
        }
        bucketIndex = Math.max(0,bucketIndex);

        Double[][] riskWeights = this.parameterSet.MapRiskClassRiskweightMap.get(this.riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue();
        double riskWeight = riskWeights[0][bucketIndex];

        double riskWeightAdjustment = this.getRiskWeightAdjustment(bucketIndex);

        return netSensi != null ? netSensi.mult(riskWeight).mult(riskWeightAdjustment).mult(concentrationRiskFactor) : null;
    }




    public RandomVariableInterface     getParameterF(RandomVariableInterface CR1, RandomVariableInterface CR2){
        RandomVariableInterface min = CR1.barrier(CR1.sub(CR2),CR2,CR1);
        RandomVariableInterface max = CR1.barrier(CR1.sub(CR2),CR1,CR2);
        return min.div(max);

    }


    public RandomVariableInterface getConcentrationRiskFactor(RandomVariableInterface netSensi, String riskFactorKey, String bucketKey, double atTime){

        double concentrationThreshold = 1.0E12;
        int bucketIndex =0;
        if (riskClassKey.equals(SIMMParameter.RiskClass.FX))
        {
            Map<String,String> FXMap = this.parameterSet.MapFXCategory;
            String category = null;
            String defaultCategory = "Category3";
            if ( riskFactorKey.length()==3 ) {
                category = FXMap.containsKey(riskFactorKey) ? FXMap.get(riskFactorKey) : defaultCategory;
            }
            else /* Usually Vega Case */
            {
                String str1=riskFactorKey.substring(0,3);
                String str2=riskFactorKey.substring(3,6);
                String category1 = FXMap.containsKey(str1) ? FXMap.get(str1) : defaultCategory;
                String category2 = FXMap.containsKey(str2) ? FXMap.get(str2) : defaultCategory;
                category = category1+"-"+category2;
            }
            concentrationThreshold = this.parameterSet.MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).get(category)[0][0];

        }
        else{
            try {
                bucketIndex = (int) Double.parseDouble(bucketKey) ;
                concentrationThreshold = this.parameterSet.MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0][bucketIndex];
            }
            catch(Exception e){
                if( bucketKey.equals("Residual")) { //!NumberUtils.isNumber(bucketKey))/*Usually RESIDUAL*/ {
                    bucketIndex = this.parameterSet.MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0].length - 1;
                    concentrationThreshold = this.parameterSet.MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0][bucketIndex];
                }
                else {
                    String key = this.parameterSet.IRCurrencyMap.get(bucketKey);
                    if ( key==null)
                        key = "High_Volatility_Currencies";
                    try {
                        concentrationThreshold = this.parameterSet.MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).get(key)[0][0];
                    }
                    catch(Exception e1){
                        concentrationThreshold = 1.0E12;
                    }
                }
            }

        }

        double riskWeightAdjustment = this.getRiskWeightAdjustment(bucketIndex);
        netSensi = netSensi.mult(riskWeightAdjustment);
        RandomVariableInterface CR = (netSensi.abs().div(concentrationThreshold)).sqrt();
        CR = CR.barrier(CR.sub(1.0), CR, 1.0);
        return CR;
    }


    public double       getRiskWeightAdjustment(int bucketIndex){
        if ( this.riskTypeKey.equals(SIMMParameter.RiskType.Vega) && (this.riskClassKey.equals(SIMMParameter.RiskClass.FX) || this.riskClassKey.equals(SIMMParameter.RiskClass.Equity) || this.riskClassKey.equals(SIMMParameter.RiskClass.Commodity))) {
            Double[][] deltaRiskWeights = this.parameterSet.MapRiskClassRiskweightMap.get(SIMMParameter.RiskType.Delta).get(riskClassKey).entrySet().iterator().next().getValue();
            double deltaRiskWeight = deltaRiskWeights[0][bucketIndex];
            double riskWeight = deltaRiskWeight;
            if ( this.parameterSet.MapHistoricalVolaRatio.containsKey(this.riskClassKey))
            {
                double historicalVolaRatio = this.parameterSet.MapHistoricalVolaRatio.get(this.riskClassKey);
                riskWeight = riskWeight * historicalVolaRatio;
            }
            return riskWeight;
        }
        else
            return 1.0;

    }

    public  RandomVariableInterface     getNetSensitivity(String riskFactorKey,String bucketKey,double evaluationTime, LIBORModelMonteCarloSimulationInterface model){

        if ( riskClassKey.equals(SIMMParameter.RiskClass.FX)) /* Sensitivities against Calculation CCY should be zero*/
            if(this.calculationCCY.equals(riskFactorKey) && riskFactorKey.length()==3){
                return new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0);
            }
        if ( riskClassKey.equals(SIMMParameter.RiskClass.CreditQ) || riskClassKey.equals(SIMMParameter.RiskClass.CreditNonQ)) {
            return Arrays.stream(this.parameterSet.CreditMaturityBuckets).map(maturityBucket->this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey, this.riskClassKey.name(),this.riskTypeKey.name(),bucketKey, maturityBucket, riskFactorKey,   evaluationTime,model)).reduce((r1, r2)->r1.add(r2)).orElseGet(()->new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0));

        }
        else {
            if (this.riskTypeKey.equals(SIMMParameter.RiskType.Delta))
                return this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey, this.riskClassKey.name(), this.riskTypeKey.name(), bucketKey,"",riskFactorKey,  evaluationTime, model);
            else if (this.riskTypeKey.equals(SIMMParameter.RiskType.Vega) && riskFactorKey.equals(SIMMParameter.inflationKey))
                return this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey, this.riskClassKey.name(), this.riskTypeKey.name(), bucketKey,"",riskFactorKey, evaluationTime, model);
            else {

                return Arrays.stream(this.parameterSet.IRMaturityBuckets).map(maturityBucket->this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey, this.riskClassKey.name(),this.riskTypeKey.name(),bucketKey, maturityBucket, riskFactorKey,   evaluationTime,model)).reduce((r1, r2)->r1.add(r2)).orElseGet(()->new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0));
            }
        }
    }

}
