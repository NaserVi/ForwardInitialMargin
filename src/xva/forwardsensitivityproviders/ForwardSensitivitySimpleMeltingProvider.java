package xva.forwardsensitivityproviders;

import net.finmath.modelling.ProductInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.HashMap;
import java.util.Map;

public class ForwardSensitivitySimpleMeltingProvider implements ForwardSensitivityProviderInterface {

    private double timeToMaturity;
    private double sensiValue;
    public ForwardSensitivitySimpleMeltingProvider()
    {

    }


    public Map<String,RandomVariableInterface> getDeltaSensitivity(double evaluationTime, String curveIndexName, LIBORModelMonteCarloSimulationInterface model)
    {
        double ttmRatio = (timeToMaturity-evaluationTime) / timeToMaturity;
        double meltedSensiValue = sensiValue*ttmRatio;
        return new HashMap<>();
    }
}
