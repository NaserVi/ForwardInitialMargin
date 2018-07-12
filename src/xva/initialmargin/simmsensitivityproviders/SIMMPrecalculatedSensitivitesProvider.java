package xva.initialmargin.simmsensitivityproviders;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;
import xva.forwardsensitivityproviders.ForwardSensitivityProviderInterface;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class SIMMPrecalculatedSensitivitesProvider {


    ForwardSensitivityProviderInterface forwardSensitivityProvider;

    public Set<SIMMTradeSpecification>    getSIMMTradeSpecs()
    {
        return null;
    }

    public SIMMPrecalculatedSensitivitesProvider(Set<SIMMTradeSpecification> tradeSet) {
    }

    Map<SIMMTradeSpecification.SensitivityKey, Double> sensitivityMap;

    public RandomVariableInterface getSIMMSensitivity(String productClass,
                                                      String riskClass,
                                                      String riskType,
                                                      String bucketKey,      // currency for IR otherwise bucket number
                                                      String maturityBucket, // only for IR and Credit risk class, null otherwise
                                                      String curveIndexName, // null if riskClass is not IR
                                                      double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException {

        Optional<SIMMTradeSpecification.SensitivityKey> optional = sensitivityMap.keySet().stream().filter(key -> key.getRiskClass().equals(riskClass) && key.getProductClass().equals(productClass) && key.getRiskType().equals(riskType) && key.getBucketKey().equals(bucketKey)).findAny();
        if (optional.isPresent()) {

            double externalProvidedSensitivity = sensitivityMap.get(optional.get());
            Map<String,RandomVariableInterface> sensitivityMap = null;
            /*if (SIMMParameter.RiskType.valueOf(riskType).equals(SIMMParameter.RiskType.Delta))
                sensitivityMap = forwardSensitivityProvider.getDeltaSensitivity(evaluationTime,model);
            else if (SIMMParameter.RiskType.valueOf(riskType).equals(SIMMParameter.RiskType.Vega))
                sensitivityMap = forwardSensitivityProvider.getVegaSensitivities(evaluationTime,model);*/



            return model.getRandomVariableForConstant(0.0);
        } else
            return model.getRandomVariableForConstant(0.0);

    }
}
