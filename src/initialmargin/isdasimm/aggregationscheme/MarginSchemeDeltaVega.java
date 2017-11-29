package initialmargin.isdasimm.aggregationscheme;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.Map;
import java.util.stream.Stream;


public class MarginSchemeDeltaVega {


    CalculationSchemeInitialMarginISDA calculationSchemeInitialMarginISDA;
    String      productClassKey;
    String      riskClassKey;
    String[]    activeBucketKeys;
    String      riskTypeKey;

    public MarginSchemeDeltaVega(CalculationSchemeInitialMarginISDA calculationSchemeInitialMarginISDA,
                                 String riskClassKey,
                                 String productClassKey, String riskTypeKey){
        this.calculationSchemeInitialMarginISDA = calculationSchemeInitialMarginISDA;
        this.riskClassKey = riskClassKey;
        this.productClassKey = productClassKey;
        this.riskTypeKey = riskTypeKey;
        //Only Modification
        this.activeBucketKeys = calculationSchemeInitialMarginISDA.getInterestRateDeltaBucketKeys();
    }

    public RandomVariableInterface getValue(double atTime){

        RandomVariableInterface deltaMargin = new RandomVariable(atTime,this.calculationSchemeInitialMarginISDA.getPathDimension(),0.0);

        Double[][] correlationMatrix = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationCrossBucketMap.get(this.riskClassKey);

        int length = correlationMatrix.length == 1 ? this.activeBucketKeys.length : correlationMatrix.length;
        RandomVariableInterface[] S1Contributions = new RandomVariableInterface[length];
        RandomVariableInterface[] KContributions = new RandomVariableInterface[length];

        if(this.activeBucketKeys.length>0) {
            for (int iBucket = 0; iBucket < activeBucketKeys.length;iBucket++)
            {
                String bucketKey =null;
                int bucketIndex = 0;
                if (riskClassKey.equals("Equity") || riskClassKey.equals("Commodity")) {
                    bucketKey = new Integer(activeBucketKeys[iBucket]).toString();
                    bucketIndex = Integer.parseInt(activeBucketKeys[iBucket]);
                }
                else {
                    bucketKey = this.activeBucketKeys[iBucket];
                    bucketIndex = iBucket;
                }

                /*Check whether we have risk factors in that bucket*/
                String[] activeRiskFactorKeys = calculationSchemeInitialMarginISDA.getMapRiskClassRiskFactors(this.riskTypeKey,bucketKey,atTime).get(riskClassKey);
                if (activeRiskFactorKeys != null && activeRiskFactorKeys.length>0) {
                    RandomVariableInterface[] netSensitivities = this.getNetSensitivities(bucketKey, activeRiskFactorKeys, atTime);
                    RandomVariableInterface[] concentrationRiskFactors = this.getConcentrationRiskFactors(bucketKey, activeRiskFactorKeys, netSensitivities, atTime);
                    RandomVariableInterface[] weightedNetSensitivities = this.getWeightedNetSensitivities(bucketKey, activeRiskFactorKeys, netSensitivities, concentrationRiskFactors, atTime);
                    RandomVariableInterface K1 = getAggregatedSensitivityForBucket(bucketKey, activeRiskFactorKeys, weightedNetSensitivities, concentrationRiskFactors, atTime);
                    RandomVariableInterface sum = this.getWeightedSensitivitySum(bucketKey, weightedNetSensitivities, atTime);
                    RandomVariableInterface S1 = K1.barrier(sum.sub(K1), K1, sum);
                    RandomVariableInterface KNegative = K1.mult(-1);
                    S1 = S1.barrier(S1.sub(KNegative), S1, KNegative);
                    S1Contributions[bucketIndex] = S1;
                    KContributions[bucketIndex] = K1;
                }
            }

            RandomVariableInterface VarCovar = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(S1Contributions, correlationMatrix);

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
        if (!this.riskClassKey.equals("FX")) {
            String bucketKey = "Residual";
            String[] activeRiskFactorKeys = calculationSchemeInitialMarginISDA.getMapRiskClassRiskFactors(this.riskTypeKey,bucketKey,atTime).get(riskClassKey);
            if (activeRiskFactorKeys != null && activeRiskFactorKeys.length>0) {
                RandomVariableInterface[] netSensitivities = this.getNetSensitivities(bucketKey, activeRiskFactorKeys, atTime);
                RandomVariableInterface[] concentrationRiskFactors = this.getConcentrationRiskFactors(bucketKey, activeRiskFactorKeys, netSensitivities, atTime);
                RandomVariableInterface[] weightedNetSensitivities = this.getWeightedNetSensitivities(bucketKey, activeRiskFactorKeys, netSensitivities, concentrationRiskFactors, atTime);
                RandomVariableInterface KResidual = this.getAggregatedSensitivityForBucket(bucketKey, activeRiskFactorKeys, weightedNetSensitivities, concentrationRiskFactors, atTime);
                deltaMargin = deltaMargin == null ? deltaMargin = KResidual : deltaMargin.add(KResidual);
            }
        }

        return deltaMargin;
    }


    private RandomVariableInterface[]   getNetSensitivities(String bucketKey,String[] activeRiskFactorKeys, double atTime){
        RandomVariableInterface[] values = new RandomVariableInterface[activeRiskFactorKeys.length];
        for (int iIndex = 0; iIndex < activeRiskFactorKeys.length; iIndex++) {
            String riskFactorKey = activeRiskFactorKeys[iIndex];
            values[iIndex]  = this.getNetSensitivity(riskFactorKey,bucketKey,atTime);

        }
        return values;

    }

    private RandomVariableInterface[]   getConcentrationRiskFactors(String bucketKey,String[] activeRiskFactorKeys, RandomVariableInterface[] netSensitivities, double atTime){
        RandomVariableInterface[] values = new RandomVariableInterface[activeRiskFactorKeys.length];
        for (int iIndex = 0; iIndex < activeRiskFactorKeys.length; iIndex++) {
            String riskFactorKey = activeRiskFactorKeys[iIndex];
            values[iIndex]  = this.getConcentrationRiskFactor(netSensitivities[iIndex],riskFactorKey,bucketKey,atTime);

        }
        return values;
    }


    private RandomVariableInterface[]   getWeightedNetSensitivities(String bucketKey,String[] activeRiskFactorKeys, RandomVariableInterface[] netSensitivities, RandomVariableInterface[] concentrationRiskFactors, double atTime){
        RandomVariableInterface[] contributions = new RandomVariableInterface[activeRiskFactorKeys.length];
        for (int iIndex = 0; iIndex < activeRiskFactorKeys.length; iIndex++) {
            String riskFactorKey = activeRiskFactorKeys[iIndex];
            RandomVariableInterface netSensi  =netSensitivities[iIndex];
            RandomVariableInterface concentrationRiskFactor = concentrationRiskFactors[iIndex];
            contributions[iIndex] = this.getWeightedNetSensitivity(netSensi,concentrationRiskFactor,riskFactorKey,bucketKey, atTime);
        }
        return contributions;

    }


    private RandomVariableInterface getAggregatedSensitivityForBucket(String bucketKey, String[] activeRiskFactorKeys, RandomVariableInterface[] weightedNetSensitivites,RandomVariableInterface[] concentrationRiskFactors, double atTime){
        RandomVariableInterface aggregatedSensi = null;
        Double correlation = 0.0;
        Double[][] correlationMatrix = new Double[ activeRiskFactorKeys.length][activeRiskFactorKeys.length];
        if (riskClassKey.equals("InterestRate") && riskTypeKey.equals("vega")) {
            correlationMatrix = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get("InterestRate_Tenor");
            RandomVariableInterface[] contributionsReDim = new RandomVariableInterface[correlationMatrix.length];
            int nTenors = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length;
            for (int iIndex = 0; iIndex < activeRiskFactorKeys.length; iIndex++)
                for (int i = 0; i < nTenors; i++) {
                    if (activeRiskFactorKeys[iIndex].equals(calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[i]))
                        contributionsReDim[i] = weightedNetSensitivites[iIndex];
                    if (activeRiskFactorKeys[iIndex].equals("inflation"))
                        contributionsReDim[contributionsReDim.length-2]=weightedNetSensitivites[iIndex];
                }
            aggregatedSensi = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(contributionsReDim,correlationMatrix);
        }
        else {
            if (riskClassKey.contains("FX"))
                correlation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][0];
            else if (riskClassKey.contains("Credit"))
                correlation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][1];
            else {
                int bucketNr = 0;
                try {
                    bucketNr = (int) Double.parseDouble(bucketKey);
                    correlation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][bucketNr];
                } catch (Exception e) {
                    bucketNr = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0].length - 1;
                    correlation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(this.riskClassKey)[0][bucketNr];
                }
            }
            for (int i = 0; i < activeRiskFactorKeys.length; i++)
                for (int j = 0; j < activeRiskFactorKeys.length; j++)
                    if ( i!=j) {
                        correlationMatrix[i][j] = getParameterF(concentrationRiskFactors[i], concentrationRiskFactors[j]).getAverage() * correlation;
                    }
            aggregatedSensi = CalculationSchemeInitialMarginISDA.getVarianceCovarianceAggregation(weightedNetSensitivites,correlationMatrix);

        }

        return aggregatedSensi;
    }


    private RandomVariableInterface   getWeightedNetSensitivity(RandomVariableInterface netSensi, RandomVariableInterface concentrationRiskFactor, String riskFactorKey ,String bucketKey, double atTime)
    {
        double riskWeight = 0;
        int bucketIndex = 0;
        try{
            bucketIndex=(int) Double.parseDouble(bucketKey);

        }
        catch(NumberFormatException e){
            bucketIndex= calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(this.riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0].length-1;
        }
        bucketIndex = Math.max(0,bucketIndex);

        Double[][] riskWeights = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(this.riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue();
        riskWeight = riskWeights[0][bucketIndex];

        double riskWeightAdjustment = this.getRiskWeightAdjustment(bucketIndex);

        if (netSensi!=null)
            return netSensi.mult(riskWeight).mult(riskWeightAdjustment).mult(concentrationRiskFactor);
        else
            return null;
    }




    public RandomVariableInterface     getParameterF(RandomVariableInterface CR1, RandomVariableInterface CR2){
        RandomVariableInterface min = CR1.barrier(CR1.sub(CR2),CR2,CR1);
        RandomVariableInterface max = CR1.barrier(CR1.sub(CR2),CR1,CR2);
        return min.div(max);

    }



    private RandomVariableInterface   getWeightedSensitivitySum(String bucketKey, RandomVariableInterface[] weightedSensitivities, double atTime){
        RandomVariableInterface aggregatedSensi = null;

        String[] activeRiskFactorKeys = calculationSchemeInitialMarginISDA.getMapRiskClassRiskFactors(this.riskTypeKey,bucketKey,atTime).get(riskClassKey);
        if (activeRiskFactorKeys == null || activeRiskFactorKeys.length==0)
            return new RandomVariable(atTime,this.calculationSchemeInitialMarginISDA.getPathDimension(),0.0);
        for (int iRiskFactor = 0; iRiskFactor < activeRiskFactorKeys.length; iRiskFactor++) {
            String key = activeRiskFactorKeys[iRiskFactor];
            RandomVariableInterface summand = weightedSensitivities[iRiskFactor];//OLD_getWeightedNetSensitivity(key, bucketKey, atTime);
            aggregatedSensi = aggregatedSensi == null ? aggregatedSensi = summand : aggregatedSensi.add(summand);
        }

        return aggregatedSensi;


    }


    public RandomVariableInterface getConcentrationRiskFactor(RandomVariableInterface netSensi, String riskFactorKey, String bucketKey, double atTime){

        double concentrationThreshold = 1.0E12;
        int bucketIndex =0;
        if (riskClassKey.equals("FX"))
        {
            Map<String,String> FXMap = calculationSchemeInitialMarginISDA.getParameterCollection().MapFXCategory;
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
            concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).get(category)[0][0];

        }
        else{
            try {
                bucketIndex = (int) Double.parseDouble(bucketKey) ;
                concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0][bucketIndex];
            }
            catch(Exception e){
                if( bucketKey.equals("Residual")) { //!NumberUtils.isNumber(bucketKey))/*Usually RESIDUAL*/ {
                    bucketIndex = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0].length - 1;
                    concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).entrySet().iterator().next().getValue()[0][bucketIndex];
                }
                else {
                    String key = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurrencyMap.get(bucketKey);
                    if ( key==null)
                        key = "High_Volatility_Currencies";
                    try {
                        concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(riskTypeKey).get(riskClassKey).get(key)[0][0];
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
        if ( this.riskTypeKey.equals("vega") && (this.riskClassKey.equals("FX") || this.riskClassKey.equals("Equity") || this.riskClassKey.equals("Commodity"))) {
            Double[][] deltaRiskWeights = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get("delta").get(riskClassKey).entrySet().iterator().next().getValue();
            double deltaRiskWeight = deltaRiskWeights[0][bucketIndex];
            double riskWeight = deltaRiskWeight;
            if ( this.calculationSchemeInitialMarginISDA.getParameterCollection().MapHistoricalVolaRatio.containsKey(this.riskClassKey))
            {
                double historicalVolaRatio = this.calculationSchemeInitialMarginISDA.getParameterCollection().MapHistoricalVolaRatio.get(this.riskClassKey);
                riskWeight = riskWeight * historicalVolaRatio;
            }
            return riskWeight;
        }
        else
            return 1.0;

    }

    public  RandomVariableInterface     getNetSensitivity(String riskFactorKey,String bucketKey,double atTime){
        RandomVariableInterface netSensi = new RandomVariable(0.0,calculationSchemeInitialMarginISDA.getPathDimension(),0.0);

        if ( riskClassKey.equals("FX")) /* Sensitivities against Calculation CCY should be zero*/
            if(calculationSchemeInitialMarginISDA.getCalculationCCY().equals(riskFactorKey) && riskFactorKey.length()==3){
                return new RandomVariable(atTime,calculationSchemeInitialMarginISDA.getPathDimension(),0.0);
            }

        if ( riskClassKey.contains("Credit")) {
            for (String maturityBucket : this.calculationSchemeInitialMarginISDA.getParameterCollection().CreditMaturityBuckets) {
                RandomVariableInterface contribution = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey, this.riskClassKey, maturityBucket, riskFactorKey, bucketKey, this.riskTypeKey, atTime);
                netSensi = netSensi == null ? netSensi = contribution : netSensi.add(contribution);
            }
        }
        else
            if ( this.riskTypeKey.equals("delta"))
                netSensi = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,"", riskFactorKey, bucketKey,this.riskTypeKey, atTime);
            else if ( this.riskTypeKey.equals("vega") && riskFactorKey.equals("inflation"))
                netSensi = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,"", riskFactorKey, bucketKey,this.riskTypeKey, atTime);
            else{
                for (String maturityBucket : this.calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets)
                {
                    RandomVariableInterface contrib = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey, this.riskClassKey, maturityBucket, riskFactorKey, bucketKey, this.riskTypeKey, atTime);
                    netSensi = netSensi==null ? contrib : netSensi.add(contrib);
                }
            }

        return netSensi;
    }

}
