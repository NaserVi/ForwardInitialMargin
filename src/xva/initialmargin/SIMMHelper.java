package xva.initialmargin;


import net.finmath.stochastic.RandomVariableInterface;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SIMMHelper {

    Set<SIMMTradeSpecification> tradeSet;

    public  SIMMHelper(Set<SIMMTradeSpecification> tradeSet){
        this.tradeSet = tradeSet;
    }



    public  static RandomVariableInterface getVarianceCovarianceAggregation(RandomVariableInterface[] contributions, Double[][] correlationMatrix){
        int i = 0;
        RandomVariableInterface value = null;
        for (RandomVariableInterface contribution1 : contributions) {
            if ( contribution1!=null) {
                value = value == null ? contribution1.squared() : value.add(contribution1.squared());
                int j = 0;
                for (RandomVariableInterface contribution2 : contributions) {
                    if (contribution2 != null && i != j) {
                        double correlation = correlationMatrix.length==1 ? correlationMatrix[0][0] : correlationMatrix[i][j];
                        RandomVariableInterface contribution = contribution1.mult(contribution2).mult(correlation);
                        value = value == null ? contribution : value.add(contribution);
                    }
                    j++;
                }
            }
            i++;
        }
        if ( value==null)
            return null;
        value = value.sqrt();
        return value;
    }


    public Set<SIMMTradeSpecification> getTadeSelection(String productClassKey, String riskClassKey, double evaluationTime){
        return tradeSet.stream().filter(trade->trade.getSensitivityKeySet(evaluationTime).stream()
                .filter(sensitivityKey -> sensitivityKey.getProductClass().equals(productClassKey) && sensitivityKey.getRiskClass().equals(riskClassKey)).findAny().isPresent()).collect(Collectors.toSet());
    }

    public Set<String>   getProductClassKeys(double evaluationTime)
    {
        Set<String> riskClasses = tradeSet.stream()
                .flatMap(
                        trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k->k!=null).map(k-> k.getProductClass().name())
                ).distinct().collect(Collectors.toSet());
        return riskClasses;
    }

    public Set<String>   getRiskClassKeysForProductClass(String productClassKey, double evaluationTime)
    {
        Set<String> riskClasses = tradeSet.stream()
                .flatMap(
                        trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k->k!=null).filter(k->k.getProductClass().equals(productClassKey)).map(k-> k.getRiskClass().name())
                ).distinct().collect(Collectors.toSet());
        return riskClasses;
    }

    public Set<String>   getRiskClassKeys(double evaluationTime)
    {
        Set<String> riskClasses = tradeSet.stream()
                .flatMap(
                        trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k->k!=null).map(k-> k.getRiskClass().name())
                ).distinct().collect(Collectors.toSet());
        return riskClasses;
    }

    public Map<String,Set<String> > getRiskClassBucketKeyMap(String riskTypeString, double evaluationTime){

        Set<String> riskClassKeys = getRiskClassKeys(evaluationTime);
        Map<String,Set<String> >     mapRiskClassBucketKeys = new HashMap<>();
        riskClassKeys.stream().forEach(riskClass ->{
            Set<String> riskFactors = tradeSet.stream()
                    .flatMap(
                            trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(  k->  k!=null && k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString)).map(k->k.getBucketKey())
                    ).distinct().collect(Collectors.toSet());
            mapRiskClassBucketKeys.put(riskClass,riskFactors);
        });
        return mapRiskClassBucketKeys;

    }

    public Map<String,Set<String>  >     getRiskClassRiskFactorMap(String riskTypeString, String bucketKey, double evaluationTime){

        Set<String> riskClassKeys = getRiskClassKeys(evaluationTime);
        if (!riskTypeString.equals("vega") ) {
            Map<String, Set<String> > mapRiskClassRiskFactorKeys = new HashMap<>();
//        if ( riskTypeString.equals("delta")) {
            riskClassKeys.stream().forEach(riskClass -> {
                Set<String> riskFactors = tradeSet.stream()
                        .flatMap(
                                trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k -> k!=null &&  k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
                        ).distinct().collect(Collectors.toSet());
                mapRiskClassRiskFactorKeys.put(riskClass, riskFactors);
            });
            return mapRiskClassRiskFactorKeys;
//        }
        }
        else{
            Map<String,Set<String> >     mapRiskClassRiskFactorKeys = new HashMap<>();
            riskClassKeys.stream().forEach(riskClass -> {
                if ( riskClass.equals("InterestRate") ) {
                    Set<String> riskFactors = tradeSet.stream()
                            .flatMap(
                                    trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k -> k!=null && k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
                            ).distinct().collect(Collectors.toSet());
                    mapRiskClassRiskFactorKeys.put(riskClass, riskFactors);
                }
                else{
                    Set<String> riskFactors = tradeSet.stream()
                            .flatMap(
                                    trade -> trade.getSensitivityKeySet(evaluationTime).stream().filter(k ->  k!=null && k.getRiskClass().equals(riskClass) && k.getRiskType().equals(riskTypeString) && k.getBucketKey().equals(bucketKey)).map(k -> k.getRiskFactorKey())
                            ).distinct().collect(Collectors.toSet());
                    mapRiskClassRiskFactorKeys.put(riskClass, riskFactors);
                }
            });
            return mapRiskClassRiskFactorKeys;
        }

    }
}
