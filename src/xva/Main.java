package xva;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.stochastic.RandomVariableInterface;
import xva.forwardsensitivityproviders.ForwardSensitivityProviderInterface;
import xva.forwardsensitivityproviders.ForwardSensitivitySimpleMeltingProvider;
import xva.initialmargin.SIMMParameter;
import xva.initialmargin.SIMMProduct;
import xva.initialmargin.simmsensitivityproviders.SIMMSensitivityProviderInterface;
import xva.initialmargin.simmsensitivityproviders.SIMMSwapBasePointDeltaSensitivityProvider;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.HashSet;
import java.util.Set;

public class Main {

    public static void main(String[] args) throws Exception{

        SIMMParameter parameterSet = new SIMMParameter();
        String calculationCCY = "EUR";


        SIMMTradeSpecification trade = new SIMMTradeSpecification(1.0E6,10.0, "Libor6M");
        SIMMTradeSpecification trade2 = new SIMMTradeSpecification(1.0E6,20.0, "Libor3M");
        Set<SIMMTradeSpecification> tradeSet = new HashSet<>();
        tradeSet.add(trade);
        tradeSet.add(trade2);

        ForwardSensitivityProviderInterface forwardSensitivityProviderInterface = new ForwardSensitivitySimpleMeltingProvider();

        SIMMSensitivityProviderInterface simmSensitivityProvider = new SIMMSwapBasePointDeltaSensitivityProvider(tradeSet,forwardSensitivityProviderInterface);

        double marginCalculationTime = 5.0;
        SIMMProduct product = new SIMMProduct(marginCalculationTime,simmSensitivityProvider,parameterSet,calculationCCY,0.0);
        LIBORMarketModel model = null;
        LIBORModelMonteCarloSimulation simulation = null;
        RandomVariableInterface result = product.getValue(4.0,simulation);


    }
}
