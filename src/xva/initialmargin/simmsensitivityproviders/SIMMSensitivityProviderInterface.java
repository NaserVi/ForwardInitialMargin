package xva.initialmargin.simmsensitivityproviders;


import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;
import xva.tradespecifications.SIMMTradeSpecification;

import java.util.Set;

public interface SIMMSensitivityProviderInterface {

    public Set<SIMMTradeSpecification>    getSIMMTradeSpecs();

    public RandomVariableInterface getSIMMSensitivity(String productClass,
                                                      String riskClass,
                                                      String riskType,
                                                      String bucketKey,      // currency for IR otherwise bucket number
                                                      String maturityBucket, // only for IR and Credit risk class, null otherwise
                                                      String curveIndexName, // null if riskClass is not IR
                                                      double evaluationTime, LIBORModelMonteCarloSimulationInterface model) ;//throws SolverException, CloneNotSupportedException, CalculationException;
}

/*
public     SIMMTradeSpecification.SensitivityKey  getSensitivityKey(String productClassKey, String riskClassKey, String riskTypeKey, String bucketKey, String maturityBucketKey, String riskFactorKey  )
    {
        if ( riskClassKey.equals("InterestRate") ) {
            if (riskFactorKey.equals("inflation") || riskFactorKey.equals("ccybasis") )
                return new SIMMTradeSpecification.SensitivityKey("",riskFactorKey,bucketKey,riskClassKey,riskTypeKey,productClassKey);
            else if ( riskTypeKey.equals("vega"))
                return new SIMMTradeSpecification.SensitivityKey(maturityBucketKey,riskFactorKey,bucketKey,riskClassKey,riskTypeKey,productClassKey);
            else
                return new SIMMTradeSpecification.SensitivityKey(maturityBucketKey,riskFactorKey,bucketKey,riskClassKey,riskTypeKey,productClassKey);
        }
        else if (riskClassKey.contains("Credit")){

            return new SIMMTradeSpecification.SensitivityKey(maturityBucketKey,riskFactorKey,bucketKey,riskClassKey,riskTypeKey,productClassKey);
        }
        else if (riskClassKey.contains("FX")){

            return !riskTypeKey.equals("vega") ?
                    new SIMMTradeSpecification.SensitivityKey("", riskFactorKey, "0", riskClassKey, riskTypeKey, productClassKey) :
                    new SIMMTradeSpecification.SensitivityKey(maturityBucketKey, riskFactorKey, "0", riskClassKey, riskTypeKey, productClassKey);
        }
        else {
            if ( !riskTypeKey.equals("vega"))
                return new SIMMTradeSpecification.SensitivityKey("",riskFactorKey,bucketKey,riskClassKey,riskTypeKey,productClassKey);
            else {
                return new SIMMTradeSpecification.SensitivityKey(maturityBucketKey, riskFactorKey, bucketKey, riskClassKey, riskTypeKey, productClassKey);

            }
        }

    }


 */