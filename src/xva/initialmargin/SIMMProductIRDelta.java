package xva.initialmargin;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;
import xva.initialmargin.simmsensitivityproviders.SIMMSensitivityProviderInterface;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class SIMMProductIRDelta extends AbstractLIBORMonteCarloProduct {
    final SIMMParameter.RiskClass riskClassKey = SIMMParameter.RiskClass.InterestRate;
    final String            riskTypeKey = "Delta";
    final String            productClassKey;
    final Set<String>       currencyKeys;
    final SIMMParameter     parameterSet;
    final SIMMHelper        helper;
    private SIMMSensitivityProviderInterface simmSensitivitivityProvider;

    public SIMMProductIRDelta(SIMMSensitivityProviderInterface simmSensitivitivityProvider, String productClassKey, SIMMParameter parameterSet, double atTime){

        this.helper = new SIMMHelper(simmSensitivitivityProvider.getSIMMTradeSpecs());
        this.productClassKey = productClassKey;
        this.currencyKeys = this.helper.getRiskClassBucketKeyMap(this.riskTypeKey,atTime).get(this.riskClassKey);
        this.parameterSet = parameterSet;
        this.simmSensitivitivityProvider = simmSensitivitivityProvider;


    }

    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){

        if (this.currencyKeys.size()==0)
            return new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0);

        RandomVariableInterface[] S1Contributions = new RandomVariableInterface[this.currencyKeys.size()];
        RandomVariableInterface[] KContributions = new RandomVariableInterface[this.currencyKeys.size()];
        int i=0;

        RandomVariableInterface[] concentrationFactors = new RandomVariableInterface[this.currencyKeys.size()];
        for (String bucketKey : this.currencyKeys)
        {
            RandomVariableInterface[][] netSensitivities = this.getNetSensitivityMatrix(bucketKey,evaluationTime,model);
            concentrationFactors[i] = getConcentrationRiskFactor(bucketKey,netSensitivities,evaluationTime,model);
            RandomVariableInterface K1 = this.getAggregatedSensitivityForBucket(bucketKey,netSensitivities,concentrationFactors[i],evaluationTime,model);
            RandomVariableInterface S1 = this.getFactorS(bucketKey,K1,netSensitivities,concentrationFactors[i],evaluationTime,model);
            S1Contributions[i] = S1;
            KContributions[i] = K1;
            i++;
        }

        double singleCorrelation = parameterSet.IRCorrelationCrossCurrency;
        Double[][] correlationMatrix = new Double[this.currencyKeys.size()][this.currencyKeys.size()];
        for (i = 0; i< currencyKeys.size(); i++)
            for (int j = 0; j< currencyKeys.size(); j++)
                if ( i!=j)
                    correlationMatrix[i][j] = getParameterG(concentrationFactors[i],concentrationFactors[j]).getAverage()*singleCorrelation;
        RandomVariableInterface VarCovar = helper.getVarianceCovarianceAggregation(S1Contributions, correlationMatrix);


        /*Adjustment on Diagonal*/
        VarCovar = VarCovar.squared();
        RandomVariableInterface SSumSQ = null;
        RandomVariableInterface KSumSQ = null;
        for ( int k = 0;k<S1Contributions.length;k++){
            SSumSQ = SSumSQ == null ? SSumSQ = S1Contributions[k].squared() : SSumSQ.add(S1Contributions[k].squared());
            KSumSQ = KSumSQ == null ? KSumSQ = KContributions[k].squared() : KSumSQ.add(KContributions[k].squared());
        }
        VarCovar = VarCovar.sub(SSumSQ).add(KSumSQ);
        RandomVariableInterface deltaMargin = VarCovar.sqrt();

        return deltaMargin;
    }


    private RandomVariableInterface[][] getNetSensitivityMatrix(String bucketKey, double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
        int nTenors = parameterSet.IRMaturityBuckets.length;
        String[] curveKeys = parameterSet.getRateCurveKeys();
        RandomVariableInterface[][] netSensitivities = new RandomVariableInterface[curveKeys.length][nTenors];
        Set<String> activeCurveKeys = helper.getRiskClassRiskFactorMap(this.riskTypeKey,bucketKey,evaluationTime).get(this.riskClassKey);
        for (int iCurve = 0; iCurve <curveKeys.length; iCurve++)
        {
            String curveKey = curveKeys[iCurve];
            if ( activeCurveKeys.contains(curveKey)) {
                for (int iTenor = 0; iTenor < nTenors; iTenor++) {
                    String maturityBucketKey = parameterSet.IRMaturityBuckets[iTenor];
                    Set<SIMMTradeSpecification> selectedTrades = helper.getTadeSelection(productClassKey,riskClassKey.name(),evaluationTime);

                    netSensitivities[iCurve][iTenor] = this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey,this.riskClassKey.name(),this.riskTypeKey,bucketKey,maturityBucketKey,curveKey,evaluationTime,model);//thiscalculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey, this.riskClassKey, iTenor, curveKey, bucketKey, "delta", atTime);
                }
            }
        }

        return netSensitivities;
    }

    private RandomVariableInterface getAggregatedSensitivityForBucket(String bucketKey, RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double evaluationTime, LIBORModelMonteCarloSimulationInterface model){

        int nTenors = parameterSet.IRMaturityBuckets.length;
        String[] curveKeys = parameterSet.getRateCurveKeys();

        int dimensionTotal=nTenors*curveKeys.length+2;
        RandomVariableInterface[] contributions = new RandomVariableInterface[dimensionTotal];

        for (int iCurve = 0; iCurve <curveKeys.length; iCurve++)
            for (int iTenor = 0; iTenor <nTenors; iTenor++)
            {
                String curveKey = curveKeys[iCurve];
                RandomVariableInterface iBucketSensi = this.getWeightedNetSensitivity(iTenor, iCurve,curveKey, bucketKey,netSensitivities,concentrationRiskFactor, evaluationTime,model);
                contributions[iCurve*nTenors+iTenor] = iBucketSensi;
            }

        RandomVariableInterface inflationSensi = this.getWeightedNetSensitivity(0,0,SIMMParameter.inflationKey,bucketKey,netSensitivities,concentrationRiskFactor,evaluationTime,model);
        RandomVariableInterface ccyBasisSensi = this.getWeightedNetSensitivity(0,0,SIMMParameter.ccyBasisKey,bucketKey,netSensitivities,concentrationRiskFactor,evaluationTime,model);
        contributions[dimensionTotal-2] = inflationSensi;
        contributions[dimensionTotal-1] = ccyBasisSensi;

        Double[][] crossTenorCorrelation = parameterSet.MapRiskClassCorrelationIntraBucketMap.get(riskClassKey);

        RandomVariableInterface aggregatedSensi = helper.getVarianceCovarianceAggregation(contributions, crossTenorCorrelation);

        return aggregatedSensi;
    }



    private RandomVariableInterface   getWeightedNetSensitivity(int iRateTenor,int iIndex, String indexName,  String bucketKey, RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor,  double evaluationTime,LIBORModelMonteCarloSimulationInterface model)
    {
        double riskWeight = 0;

        if (!indexName.equals(SIMMParameter.inflationKey) && !indexName.equals(SIMMParameter.ccyBasisKey))
        {
            Optional<Map.Entry<String, String>> optional = parameterSet.IRCurrencyMap.entrySet().stream().filter(entry -> entry.getKey().contains(bucketKey)).findAny();
            String currencyMapKey;
            currencyMapKey = !optional.isPresent() ? "High_Volatility_Currencies" : optional.get().getValue();
            currencyMapKey = currencyMapKey.replace("_Traded", "").replace("_Well", "").replace("_Less", "");
            Double[] riskWeights = parameterSet.MapRiskClassRiskweightMap.get(riskTypeKey).get(SIMMParameter.RiskClass.InterestRate).get(currencyMapKey)[0];
            riskWeight = riskWeights[iRateTenor];
            RandomVariableInterface netSensi =  netSensitivities[iIndex][iRateTenor];
            return netSensi != null ?
                    netSensi.mult(riskWeight).mult(concentrationRiskFactor) :
                    new RandomVariable(evaluationTime, model.getNumberOfPaths(), 0.0);
        }
        else { /* Inflation or CCYBasis*/
            riskWeight = parameterSet.MapRiskClassRiskweightMap.get(riskTypeKey).get(SIMMParameter.RiskClass.InterestRate).get(indexName)[0][0];

            String maturityBucketKey = this.parameterSet.IRMaturityBuckets[iRateTenor];
            RandomVariableInterface netSensi =  this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey,this.riskClassKey.name(),this.riskTypeKey,bucketKey,maturityBucketKey, indexName, evaluationTime,model);
            if (netSensi!=null)
            {
                netSensi = netSensi.mult(riskWeight);
                if ( !indexName.equals(SIMMParameter.ccyBasisKey)) {
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



    public RandomVariableInterface getFactorS(String bucketKey, RandomVariableInterface K,RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor, double evaluationTime,LIBORModelMonteCarloSimulationInterface model){
        RandomVariableInterface sum = this.getWeightedSensitivitySum(bucketKey,netSensitivities,concentrationRiskFactor, evaluationTime,model);
        RandomVariableInterface S1 = K.barrier(sum.sub(K),K,sum);
        RandomVariableInterface KNegative = K.mult(-1);
        S1 = S1.barrier(S1.sub(KNegative),S1,KNegative);
        return S1;
    }

    private RandomVariableInterface   getWeightedSensitivitySum(String bucketKey,RandomVariableInterface[][] netSensitivities,RandomVariableInterface concentrationRiskFactor,double evaluationTime,LIBORModelMonteCarloSimulationInterface model){
        RandomVariableInterface aggregatedSensi = null;
        String[] curveKeys = parameterSet.getRateCurveKeys();

        for (int iIndex = 0; iIndex < curveKeys.length; iIndex++) {
            for (int iTenor = 0; iTenor < parameterSet.IRMaturityBuckets.length; iTenor++) {
                String key = curveKeys[iIndex];
                RandomVariableInterface summand = getWeightedNetSensitivity(iTenor,iIndex, key, bucketKey,netSensitivities,concentrationRiskFactor, evaluationTime,model);
                if (summand !=null)
                    aggregatedSensi = aggregatedSensi != null ? aggregatedSensi.add(summand) : new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0);
            }
        }
        RandomVariableInterface inflationSensi = this.getWeightedNetSensitivity(0,0,SIMMParameter.inflationKey,bucketKey,netSensitivities,concentrationRiskFactor,evaluationTime,model);
        RandomVariableInterface ccyBasisSensi = this.getWeightedNetSensitivity(0,0,SIMMParameter.ccyBasisKey,bucketKey,netSensitivities,concentrationRiskFactor,evaluationTime,model);
        return aggregatedSensi.add(inflationSensi.add(ccyBasisSensi));
    }


    public RandomVariableInterface getConcentrationRiskFactor(String bucketKey, RandomVariableInterface[][] netSensitivities,double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
        RandomVariableInterface sensitivitySum = null;
        String[] curveKeys = parameterSet.getRateCurveKeys();
        for (int iIndex = 0; iIndex < curveKeys.length; iIndex++) {
            for (int iTenor = 0; iTenor < parameterSet.IRMaturityBuckets.length; iTenor++) {
                String key = curveKeys[iIndex];
                RandomVariableInterface summand = netSensitivities[iIndex][iTenor];
                if (summand !=null)
                    sensitivitySum = sensitivitySum != null ? sensitivitySum.add(summand) : new RandomVariable(evaluationTime,model.getNumberOfPaths(),0.0);
            }
        }
        RandomVariableInterface inflationSensi = this.simmSensitivitivityProvider.getSIMMSensitivity(this.productClassKey,this.riskClassKey.name(),this.riskTypeKey,bucketKey,"","inflation",evaluationTime,model);
        if ( sensitivitySum !=null && inflationSensi !=null)
            sensitivitySum = sensitivitySum.add(inflationSensi); // Inflation Sensi are included in Sum, CCYBasis not

        Optional<Map.Entry<String,String> > optional = parameterSet.IRCurrencyMap.entrySet().stream().filter(entry->entry.getKey().contains(bucketKey)).findAny();
        String currencyMapKey;
        currencyMapKey = !optional.isPresent() ? "High_Volatility_Currencies" : optional.get().getValue();

        double concentrationThreshold = parameterSet.MapRiskClassThresholdMap.get(this.riskTypeKey).get(riskClassKey).get(currencyMapKey)[0][0];
        RandomVariableInterface CR = (sensitivitySum.abs().div(concentrationThreshold)).sqrt();
        CR = CR.barrier(CR.sub(1.0), CR, 1.0);
        return CR;
    }

}
