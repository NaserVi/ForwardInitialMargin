package xva.forwardsensitivityproviders;

import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.Map;

public interface ForwardSensitivityProviderInterface {

    public Map<String,RandomVariableInterface> getDeltaSensitivity(double evaluationTime, String curveKey, LIBORModelMonteCarloSimulationInterface model);



}
