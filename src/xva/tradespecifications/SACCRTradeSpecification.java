package xva.tradespecifications;

import net.finmath.montecarlo.AbstractMonteCarloProduct;

import java.time.LocalDate;

public class SACCRTradeSpecification {

    AbstractMonteCarloProduct underlyingValuationProduct;


    public enum CreditRating{
        AAA,
        AA,
        A,
        BBB,
        BB,
        B,
        CCC,
        NotDefined
    }



    public enum tRegulatoryAssetClass
    {
        InterestRate,
        FX,
        Credit,
        EquitySingleName,
        EquityIndex,
        CreditSingleName,
        CreditIndex,
        Commodity,
        Default
    }



    private String                  ccy;
    double[]                        notionals;
    double[]                        periodStartTimes;
    double[]                        periodEndTimes;
    private String                  underlyingCreditRating;
    private LocalDate               maturityDateUnderlying;
    private LocalDate               tradeStartDate;
    private LocalDate               tradeMaturityDate;
    private String                  regulatoryRiskFactorKey;
    private tRegulatoryAssetClass   assetClass;


    public SACCRTradeSpecification(){

    }


    public  AbstractMonteCarloProduct   getUnderlyingValuationProduct(){
        return this.underlyingValuationProduct;
    }

    public String getRegulatoryRiskFactorKey() {
        return regulatoryRiskFactorKey;
    }

    public tRegulatoryAssetClass getAssetClass() {
        return assetClass;
    }


    public CreditRating getUnderlyingCreditRating(){
        return CreditRating.AAA;
    }

    public  double  getRegulatoryPosition(double evaluationTime){
        return 1.0;
    }

    public double getTimeToMaturity(double atTime) {
        return 0.0;
    }

    public double getTimeToMaturityUnderlying(double atTime) {
        return 0.0;
    }

    public double getTimeToStartDate(double atTime) {
        return 0.0;
    }

    public double getMaxNotional(double evaluationTime) {
        if (this.getTimeToMaturity(evaluationTime) <= 0.0)
            return 0.0;
        if (this.periodStartTimes != null && this.periodStartTimes.length > 0) {
            double maxNotional = 0.0;
            for (int i = 0; i < this.periodEndTimes.length; i++) {
                if (this.periodStartTimes[i] >= evaluationTime) // For all remaining notionals beyond atTime
                {
                    if (this.notionals[i] > maxNotional)
                        maxNotional = this.notionals[i];
                }
            }
            return maxNotional;
        }
        return this.notionals[0];
    }

    public double getAverageNotional(double evaluationTime) {
        if (this.getTimeToMaturity(evaluationTime) <= 0.0)
            return 0.0;
        if (this.periodStartTimes != null && this.periodStartTimes.length > 0) {
            double averageNotional = 0.0;
            double sumTimeDiff = 0.0;
            for (int i = 0; i < this.periodEndTimes.length; i++) {
                if (this.periodEndTimes[i] >= evaluationTime) // For all remaining notionals beyond atTime
                {
                    double periodLength = 0.0;
                    periodLength =
                            this.periodStartTimes[i] >= evaluationTime ? this.periodEndTimes[i] - this.periodStartTimes[i] : this.periodEndTimes[i] - evaluationTime;
                    periodLength = java.lang.Math.max(periodLength, 1.0E-12);
                    if (periodLength < 0.0)
                        throw new RuntimeException("getAverageNotional: Period Length should be positive");
                    averageNotional = averageNotional + periodLength * this.notionals[i];
                    sumTimeDiff += periodLength;
                }
            }
            if (sumTimeDiff < 1.0E-12) {
                return this.notionals[0];
            }
            if(sumTimeDiff == 0.0) {
                throw new RuntimeException("sumTimeDiff must not be 0.0");
            }
            averageNotional = averageNotional / sumTimeDiff;
            return averageNotional;
        }
        return this.notionals[0];
    }

}
