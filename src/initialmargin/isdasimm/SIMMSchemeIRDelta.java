package initialmargin.isdasimm;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.stream.Stream;


public class SIMMSchemeIRDelta {
	
    SIMMSchemeMain calculationSchemeInitialMarginISDA;
    String      productClassKey;
    String      riskClassKey;
    String[]    bucketKeys;
    final String riskTypeKey = "delta";
    private Map<String/*bucketKey*/, RandomVariableInterface> concentrationRiskFactorMap = new HashMap<String/*bucketKey*/, RandomVariableInterface>(); 

    public SIMMSchemeIRDelta(SIMMSchemeMain calculationSchemeInitialMarginISDA,
                               String productClassKey){
        this.calculationSchemeInitialMarginISDA = calculationSchemeInitialMarginISDA;
        this.riskClassKey = "InterestRate";
        this.productClassKey = productClassKey;
        this.bucketKeys = calculationSchemeInitialMarginISDA.getInterestRateDeltaBucketKeys();
    }

    public RandomVariableInterface getValue(double atTime){

        if (this.bucketKeys.length==0)
            return new RandomVariable(atTime,this.calculationSchemeInitialMarginISDA.getPathDimension(),0.0);

        RandomVariableInterface[] S1Contributions = new RandomVariableInterface[this.bucketKeys.length];
        RandomVariableInterface[] KContributions = new RandomVariableInterface[this.bucketKeys.length];
        int i=0;
        for (String bucketKey : this.bucketKeys)
        {
            RandomVariableInterface K1 = this.getAggregatedSensitivityForBucket(bucketKey,atTime);
            RandomVariableInterface S1 = this.getFactorS(bucketKey,K1,atTime);
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
                correlationMatrix[i][j] = getParameterG(bucketKeys[i],bucketKeys[j],atTime).getAverage()*singleCorrelation;
        VarCovar = SIMMSchemeMain.getVarianceCovarianceAggregation(S1Contributions, correlationMatrix);


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

        if ( Double.isNaN(deltaMargin.get(0)))
            System.out.print("");
        return deltaMargin;
    }



    private RandomVariableInterface getAggregatedSensitivityForBucket(String bucketKey, double atTime){
        RandomVariableInterface aggregatedSensi = null;


        //if ( riskClassKey.equals("InterestRate")) {
        /**
         * Changed to first across curves with sub curve correlation
         * Then across tenors
         * therefore riskfactor is tenor,
         */

        int nTenors = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length;
        int nCurves = calculationSchemeInitialMarginISDA.getIRCurveIndexNames().length;

        // Inserted 
        RandomVariableInterface[][] subCurveContributions = new RandomVariableInterface[nCurves+2][nTenors];
        
        int dimensionTotal=nTenors*nCurves+2;
        RandomVariableInterface[] contributions = new RandomVariableInterface[dimensionTotal];

        for (int iCurve = 0; iCurve <nCurves; iCurve++)
            for (int iTenor = 0; iTenor <nTenors; iTenor++)
            {
                String curveKey = calculationSchemeInitialMarginISDA.getIRCurveIndexNames()[iCurve];
                RandomVariableInterface iBucketSensi = this.getWeightedNetSensitivity(iTenor, curveKey, bucketKey, atTime);
                contributions[iCurve*nTenors+iTenor] = iBucketSensi;
                // Inserted 
                subCurveContributions[iCurve][iTenor]=iBucketSensi;
            }

        // use this if the inflation sensitivity is calculated by doCalculateDeltaSensitivityIR
        RandomVariableInterface inflationSensi = new RandomVariable(0.0);//this.getWeightedNetSensitivity(0,"inflation",bucketKey,atTime);
        RandomVariableInterface ccyBasisSensi =  new RandomVariable(0.0);//this.getWeightedNetSensitivity(0,"ccybasis",bucketKey,atTime);
        contributions[dimensionTotal-2] = inflationSensi;
        contributions[dimensionTotal-1] = ccyBasisSensi;
        
        // Inserted
        subCurveContributions[nCurves][0]=inflationSensi;
        subCurveContributions[nCurves+1][0]=ccyBasisSensi;

        Double[][] crossTenorCorrelation = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassCorrelationIntraBucketMap.get(riskClassKey);

        aggregatedSensi = getVarianceCovarianceAggregation(subCurveContributions, crossTenorCorrelation);//SIMMSchemeMain.getVarianceCovarianceAggregation(contributions, crossTenorCorrelation);
        //}    

        if ( Double.isNaN(aggregatedSensi.get(0)))
            System.out.print("");
        return aggregatedSensi;
    }
    
    
 // new variance aggregation
    private RandomVariableInterface  getVarianceCovarianceAggregation(RandomVariableInterface[][] subCurveContribution, Double[][] correlation){
        double subCurveCorrelation = 0.98;
        int nTenors = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length;
        int nCurves = calculationSchemeInitialMarginISDA.getIRCurveIndexNames().length;
        RandomVariableInterface value = new RandomVariable(0.0);
        
        for(int curveIndex1=0;curveIndex1<nCurves;curveIndex1++){
        		for(int tenorIndex1=0;tenorIndex1<nTenors;tenorIndex1++){
        			if(subCurveContribution[curveIndex1][tenorIndex1]!=null){
             		   value=value.add(subCurveContribution[curveIndex1][tenorIndex1].squared());
             		}
        			
        			for(int curveIndex2=0;curveIndex2<nCurves;curveIndex2++){
        			for(int tenorIndex2=0;tenorIndex2<nTenors;tenorIndex2++){ 
        				if(subCurveContribution[curveIndex1][tenorIndex1]!=null && subCurveContribution[curveIndex1][tenorIndex2]!=null){
        			
        				if(curveIndex1==curveIndex2 && (tenorIndex1!=tenorIndex2)){
        					value=value.add(subCurveContribution[curveIndex1][tenorIndex1].mult(subCurveContribution[curveIndex1][tenorIndex2]).mult(correlation[tenorIndex1][tenorIndex2]));
        				} else if(curveIndex1!=curveIndex2){
        					value=value.add(subCurveContribution[curveIndex1][tenorIndex1].mult(subCurveContribution[curveIndex2][tenorIndex2]).mult(correlation[tenorIndex1][tenorIndex2]*subCurveCorrelation));
        				}
        			    }
        		    
        		}
        	}
        }
        }
        // inflation and ccy risk factors
        value=value.add(subCurveContribution[nCurves-2][0].squared()).add(subCurveContribution[nCurves-1][0].squared());
        //...
        return value.sqrt();
    
    }



    private RandomVariableInterface   getWeightedNetSensitivity(int iRateTenor,String indexName,String bucketKey, double atTime)
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
            currencyMapKey = currencyMapKey.replace("_Traded", "").replace("_Well", "").replace("_Less", ""); // no replacement takes place ?!
            Double[] riskWeights = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(riskTypeKey).get("InterestRate").get(currencyMapKey)[0];
            riskWeight = riskWeights[iRateTenor];
        }
        else {
            riskWeight = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassRiskweightMap.get(riskTypeKey).get("InterestRate").get(indexName)[0][0];
        }

        String maturityBucket = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[iRateTenor];
        RandomVariableInterface netSensi =  calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,maturityBucket, indexName, bucketKey,this.riskTypeKey, atTime);
        if (netSensi!=null) {
            if ( indexName.equals("ccybasis"))
                return netSensi.mult(riskWeight);
            else {
                RandomVariableInterface concentrationRiskFactor = this.getConcentrationRiskFactor(bucketKey,atTime);
                //if(netSensi.abs().getAverage()>0)System.out.println("RW: " + riskWeight);
                return netSensi.mult(riskWeight).mult(concentrationRiskFactor);
                
            }
        }
        else
            return null;
    }

    public RandomVariableInterface     getParameterG(String bucketKey1, String bucketKey2, double atTime){
        RandomVariableInterface CR1 = this.getConcentrationRiskFactor(bucketKey1,atTime);
        RandomVariableInterface CR2 = this.getConcentrationRiskFactor(bucketKey2,atTime);
        RandomVariableInterface min = CR1.barrier(CR1.sub(CR2),CR2,CR1);
        RandomVariableInterface max = CR1.barrier(CR1.sub(CR2),CR1,CR2);
        return min.div(max);

    }



    public RandomVariableInterface getFactorS(String bucketKey, RandomVariableInterface K, double atTime){
        RandomVariableInterface sum = this.getWeightedSensitivitySum(bucketKey, atTime);
        RandomVariableInterface S1 = K.barrier(sum.sub(K),K,sum);
        RandomVariableInterface KNegative = K.mult(-1);
        S1 = S1.barrier(S1.sub(KNegative),S1,KNegative);
        return S1;
    }

    private RandomVariableInterface   getWeightedSensitivitySum(String bucketKey, double atTime){
        RandomVariableInterface aggregatedSensi = null;

        for (int iIndex = 0; iIndex < calculationSchemeInitialMarginISDA.getIRCurveIndexNames().length; iIndex++) {
            for (int iTenor = 0; iTenor < calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length; iTenor++) {
                String key = calculationSchemeInitialMarginISDA.getIRCurveIndexNames()[iIndex];
                RandomVariableInterface summand = getWeightedNetSensitivity(iTenor, key, bucketKey, atTime);
                aggregatedSensi = aggregatedSensi == null ? aggregatedSensi = summand : aggregatedSensi.add(summand);
            }
        }
        RandomVariableInterface inflationSensi = this.getWeightedNetSensitivity(0,"inflation",bucketKey,atTime);
        RandomVariableInterface ccyBasisSensi = this.getWeightedNetSensitivity(0,"ccybasis",bucketKey,atTime);
        return aggregatedSensi.add(inflationSensi.add(ccyBasisSensi));


    }

    // Inserted 
    public RandomVariableInterface getConcentrationRiskFactor(String bucketKey, double atTime){
    	if(!concentrationRiskFactorMap.containsKey(bucketKey)){
    		doCalculateConcentrationRiskFactor(bucketKey,atTime);
    	}
    	return concentrationRiskFactorMap.get(bucketKey);
    }
    
    private void doCalculateConcentrationRiskFactor(String bucketKey,double atTime){
        RandomVariableInterface sensitivitySum = null;
        for (int iIndex = 0; iIndex < calculationSchemeInitialMarginISDA.getIRCurveIndexNames().length; iIndex++) {
            for (int iTenor = 0; iTenor < calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets.length; iTenor++) {
                String key = calculationSchemeInitialMarginISDA.getIRCurveIndexNames()[iIndex];
                String maturityBucket = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[iTenor];
                RandomVariableInterface summand = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,maturityBucket,key,bucketKey,"delta",atTime);//"ccybasis",bucketKey,"delta",atTime);//getWeightedNetSensitivity(iTenor, key, bucketKey, atTime);
                sensitivitySum = sensitivitySum == null ? sensitivitySum = summand : sensitivitySum.add(summand);
            }
        }
        String firstBucket = calculationSchemeInitialMarginISDA.getParameterCollection().IRMaturityBuckets[0];
        RandomVariableInterface inflationSensi = calculationSchemeInitialMarginISDA.getNetSensitivity(this.productClassKey,this.riskClassKey,firstBucket,"inflation",bucketKey,"delta",atTime);
        sensitivitySum = sensitivitySum.add(inflationSensi); // Inflation Sensi are included in Sum, CCYBasis not

//        Optional<Map.Entry<String,String> > optional = calculationSchemeInitialMarginISDA.getParameterCollection().IRCurrencyMap.entrySet().stream().filter(entry->entry.getKey().contains(bucketKey)).findAny();
//        String currencyMapKey;
//        if (!optional.isPresent())
//            currencyMapKey="High_Volatility_Currencies";
//        else
//            currencyMapKey = optional.get().getValue();

        //double concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(this.riskTypeKey).get(riskClassKey).get(currencyMapKey)[0][0];
        double concentrationThreshold = calculationSchemeInitialMarginISDA.getParameterCollection().MapRiskClassThresholdMap.get(this.riskTypeKey).get(riskClassKey).get(bucketKey)[0][0];
        RandomVariableInterface CR = (sensitivitySum.abs().div(concentrationThreshold)).sqrt();
        CR = CR.barrier(CR.sub(1.0), CR, 1.0);
        concentrationRiskFactorMap.put(bucketKey, CR);
    }



}
