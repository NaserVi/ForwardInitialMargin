package xva.initialmargin.simmsensitivityproviders;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;

import net.finmath.stochastic.RandomVariableInterface;
import xva.forwardsensitivityproviders.ForwardSensitivityProviderInterface;
import xva.initialmargin.SIMMParameter;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class SIMMSwapBasePointDeltaSensitivityProvider implements SIMMSensitivityProviderInterface {

    final SIMMParameter.ProductClass productClass = SIMMParameter.ProductClass.RatesFX;
    final SIMMParameter.RiskClass riskClass = SIMMParameter.RiskClass.InterestRate;
    final SIMMParameter.RiskType riskType = SIMMParameter.RiskType.Delta;

    Map<SIMMTradeSpecification.SensitivityKey,Double> notionalMap;

    ForwardSensitivityProviderInterface forwardSensitivityProvider;

    public Set<SIMMTradeSpecification> getSIMMTradeSpecs()
    {
        return null;
    }

    public SIMMSwapBasePointDeltaSensitivityProvider(Set<SIMMTradeSpecification> tradeSet, ForwardSensitivityProviderInterface forwardSensiProvider) throws Exception{
        // BUILD UP MAP
        // For Each Key: Build a ForwardSensitivityProvider
        this.forwardSensitivityProvider = forwardSensiProvider;
    }

    public RandomVariableInterface getSIMMSensitivity(  String productClass,
                                                        String riskClass,
                                                        String riskType,
                                                        String bucketKey,      // currency for IR otherwise bucket number
                                                        String maturityBucket, // only for IR and Credit risk class, null otherwise
                                                        String curveIndexName, // null if riskClass is not IR
                                                        double evaluationTime, LIBORModelMonteCarloSimulationInterface model)
    {


        Optional<SIMMTradeSpecification.SensitivityKey> optional = notionalMap.keySet().stream().filter(key->key.getRiskClass().equals(riskClass) && key.getProductClass().equals(productClass) && key.getRiskType().equals(riskType) && key.getBucketKey().equals(bucketKey)).findAny();
        if (optional.isPresent()){
            double effectiveMaturity = 0.0;//maturityMap.get(curveIndexName);
            double notional = notionalMap.get(curveIndexName);
            double bpvProxyAtZero = effectiveMaturity*notional/10000.;
            return new RandomVariable(evaluationTime,model.getNumberOfPaths(),bpvProxyAtZero);
        }
        else
            return model.getRandomVariableForConstant(0.0);

    }

}
