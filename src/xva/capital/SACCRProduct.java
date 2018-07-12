package xva.capital;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariableInterface;
import xva.tradespecifications.SACCRTradeSpecification;

import java.util.*;
import java.util.stream.Collectors;

import static java.lang.Math.max;

public class SACCRProduct extends AbstractLIBORMonteCarloProduct {

    final public static double      SACCRAlphaMultiplier = 1.4;

    Collection<SACCRTradeSpecification>   tradeCollection;
    boolean     isCollateralized;
    double      marginPeriodOfRiskFractionOfYears;

    double     capitalReferenceTime;


    public  Collection<SACCRTradeSpecification>   getTradeSelection(SACCRTradeSpecification.tRegulatoryAssetClass assetClass, String riskFactor1, String riskFactor2){
        throw new RuntimeException();
    }


    public SACCRProduct(double capitalReferenceTime){

    }


    RandomVariableInterface      getUnderlyingNetValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
        return tradeCollection.stream().map(trade->{
            try{
                return trade.getUnderlyingValuationProduct().getValue(evaluationTime,model);
            }
            catch(Exception e){
                return null;
            }
        }).filter(v->v!=null).reduce((v1,v2)->v1.add(v2)).orElseGet(()->model.getRandomVariableForConstant(0.0));
    }

    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {

        if (evaluationTime> capitalReferenceTime)
            return model.getRandomVariableForConstant(0.0);

        RandomVariableInterface V = this.getUnderlyingNetValue(evaluationTime,model);
        RandomVariableInterface ReplacementCost = V.floor(0.0);
        RandomVariableInterface pfe = this.getRegulatoryPotentialFutureExposure(V,evaluationTime);
        ReplacementCost = ReplacementCost.add( pfe );
        ReplacementCost = ReplacementCost.mult(SACCRAlphaMultiplier);

        /*Discounted ReplacementCost*/
        ReplacementCost = ReplacementCost.mult(model.getNumeraire(evaluationTime)).div(model.getNumeraire(capitalReferenceTime));

        return ReplacementCost;

    }

    // Paragraph 146
    public		RandomVariableInterface		getRegulatoryPotentialFutureExposure(RandomVariableInterface marketValue, double atTime)
    {
        double addOn = getAggregateNotionalAddOn(atTime);
        RandomVariableInterface V = marketValue;    /* May be adjusted for Collateral */

        RandomVariableInterface multiplier = null;
        multiplier = this.getMultiplier(marketValue,addOn,atTime);

        RandomVariableInterface pfe = multiplier.mult(addOn);
        return pfe;

    }

    public 		RandomVariableInterface		getMultiplier(RandomVariableInterface marketValue, double addOn, double atTime)
    {
        double floor = 0.05;
        double divisor = 2.0 * (1.0-floor)*addOn;
        if( divisor < 1.0E-9 ) divisor = 1.0;
        RandomVariableInterface multiplier = marketValue.div(divisor);
        multiplier = multiplier.exp();
        multiplier = multiplier.mult(1.0-floor);
        multiplier = multiplier.add(floor);
        multiplier = multiplier.cap(1.0);	// Floored at 1.0
        return multiplier;
    }


    // Paragrap 150
    public 		double		getAggregateNotionalAddOn(double atTime)
    {

        double addOn = this.getAddOnForAssetClass(atTime, SACCRTradeSpecification.tRegulatoryAssetClass.InterestRate);
        addOn +=  this.getAddOnForAssetClass(atTime,SACCRTradeSpecification.tRegulatoryAssetClass.FX);
        addOn +=  this.getAddOnForAssetClass(atTime,SACCRTradeSpecification.tRegulatoryAssetClass.Credit);
        addOn +=  this.getAddOnForAssetClass(atTime,SACCRTradeSpecification.tRegulatoryAssetClass.EquityIndex);
        addOn +=  this.getAddOnForAssetClass(atTime,SACCRTradeSpecification.tRegulatoryAssetClass.Commodity);
        return addOn;
    }


    public 		double		getAddOnForAssetClass(double atTime, SACCRTradeSpecification.tRegulatoryAssetClass assetClass)
    {
        double addOn =0.0;

        /* All trades for asset class*/
        Collection<SACCRTradeSpecification>    tradeCollection = this.tradeCollection.stream().filter(trade->trade.getAssetClass()==assetClass).collect(Collectors.toList());

        /*All risk factors*/
        List<String> regulatoryRiskFactors = tradeCollection.stream().map(trade->trade.getRegulatoryRiskFactorKey()).collect(Collectors.toList());


        return addOn;
    }



    // Paragraph 169 / 171 /172
    public 		double		getAddOnForHedgeSet(SACCRTradeSpecification.tRegulatoryAssetClass assetClass,  double evaluationTime)
    {

        double supervisoryFactor = 1.0E12;
        double effectiveNotional = 1.0E12;


        Set<SACCRTradeSpecification>    tradeSelectionAssetClass = tradeCollection.stream().filter(trade->trade.getAssetClass().equals(assetClass)).collect(Collectors.toSet());
        if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.InterestRate )
        {
            supervisoryFactor = 0.005;
            // Paragraph 169: Sum over all Buckets in relevant hedgesets
            double B1 = this.getEffectiveNotional(tradeSelectionAssetClass.stream().filter(trade->trade.getTimeToMaturity(evaluationTime)<1.0).collect(Collectors.toList()),evaluationTime );
            double B2 = this.getEffectiveNotional(tradeSelectionAssetClass.stream().filter(trade->trade.getTimeToMaturity(evaluationTime)>1.0 && trade.getTimeToMaturity(evaluationTime)<5.0).collect(Collectors.toList()),evaluationTime );
            double B3 = this.getEffectiveNotional(tradeSelectionAssetClass.stream().filter(trade->trade.getTimeToMaturity(evaluationTime)>5.0).collect(Collectors.toList()),evaluationTime);

            effectiveNotional = B1*B1 +  B2*B2 + B3*B3 + 1.4*B1*B2 + 1.4*B2*B3 + 0.6*B1*B3;
            effectiveNotional = Math.sqrt(effectiveNotional);

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.FX )
        {
            supervisoryFactor = 0.04;
            // Paragraph 171: FX: Take absolute amount
            effectiveNotional =  Math.abs(this.getEffectiveNotional(tradeSelectionAssetClass,evaluationTime));

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.EquitySingleName  )
        {
            supervisoryFactor = 0.32;
            effectiveNotional = this.getEffectiveNotional(tradeSelectionAssetClass,evaluationTime);

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.EquityIndex  )
        {
            supervisoryFactor = 0.2;
            effectiveNotional =  this.getEffectiveNotional(tradeSelectionAssetClass,evaluationTime);

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.Commodity  )
        {
            supervisoryFactor = 0.18;
            effectiveNotional =  this.getEffectiveNotional(tradeSelectionAssetClass,evaluationTime);

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.CreditSingleName )
        {
            for ( SACCRTradeSpecification trade : tradeSelectionAssetClass )
            {
                String creditRating = trade.getUnderlyingCreditRating().name();
                if ( creditRating.equals("AAA"))
                    supervisoryFactor = 0.0038;
                else if ( creditRating.equals("AA"))
                    supervisoryFactor = 0.0042;
                else if ( creditRating.equals("A"))
                    supervisoryFactor = 0.0054;
                else if ( creditRating.equals("BBB") )
                    supervisoryFactor = 0.0106;
                else if ( creditRating.equals("BB"))
                    supervisoryFactor = 0.016;
                else if ( creditRating.equals("N/A") )
                    supervisoryFactor = 0.06;
                else
                    throw new RuntimeException("Underlying Credit Rating not def!");
            }
            effectiveNotional =  this.getEffectiveNotional(tradeCollection,evaluationTime);

        }
        else if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.CreditIndex )
        {
            supervisoryFactor = 0.0106; // Missing InvestmentGradeFlag
            effectiveNotional =  this.getEffectiveNotional(tradeCollection,evaluationTime);

        }
        else
        {
            supervisoryFactor = 0; // Missing InvestmentGradeFlag
            effectiveNotional =  0;
        }

        return effectiveNotional * supervisoryFactor;
    }


    // 183
    public   	double		getSupervisoryCorrelation(SACCRTradeSpecification.tRegulatoryAssetClass assetClass)
    {
        if(assetClass == SACCRTradeSpecification.tRegulatoryAssetClass .CreditSingleName )
            return 0.5;
        else if(assetClass == SACCRTradeSpecification.tRegulatoryAssetClass .CreditIndex )
            return 0.8;
        else if(assetClass == SACCRTradeSpecification.tRegulatoryAssetClass .EquitySingleName )
            return 0.5;
        else if(assetClass == SACCRTradeSpecification.tRegulatoryAssetClass .EquityIndex )
            return 0.8;
        else if (assetClass == SACCRTradeSpecification.tRegulatoryAssetClass .Commodity )
            return 0.4;
        else
            throw new RuntimeException("getSupervisoryCreditCorrelation: Not defined");
    }


    // Paragraph 168
    private		double		getEffectiveNotional(Collection<SACCRTradeSpecification> tradeCollection, double evaluationTime )
    {

        double sum = 0.0;
        for ( SACCRTradeSpecification trade : tradeCollection )
        {
            double addOn = this.getTradeAddon(trade,evaluationTime);
            sum = sum + addOn;
        }
        return sum;
    }

    public 		double		getSACCRNotional(SACCRTradeSpecification trade, double evaluationTime)
    {
        return trade.getAverageNotional(evaluationTime);
    }

    public		double		getTradeAddon(SACCRTradeSpecification trade, double evaluationTime)
    {

        double SD = this.getSupervisoryDuration(trade.getAssetClass(), trade.getTimeToStartDate(evaluationTime),trade.getTimeToMaturity(evaluationTime));
        double MF = this.getRegulatoryMaturityScaleFactor(trade, evaluationTime);
        double position = trade.getRegulatoryPosition(evaluationTime);
        double SACCRNotional = getSACCRNotional(trade,evaluationTime);
        double adddon = SD * MF * position  * SACCRNotional;
        return adddon;

    }

    public		double	 getRegulatoryMaturityScaleFactor(SACCRTradeSpecification trade, double evaluationTime)
    {
        if ( isCollateralized)
        {
            double MPORMIN = 10.0 / 250.0;
            double MPOR  = marginPeriodOfRiskFractionOfYears;
            MPOR = max(MPOR,MPORMIN);
            return 1.5 * Math.sqrt(MPOR);
        }
        else
        {
            double maturity = max(10.0/250.0, trade.getTimeToMaturity(evaluationTime));
            double MF = Math.min(maturity, 1.0);
            MF = MF / 1.0;
            return Math.sqrt(MF);
        }
    }


    public		double	 getSupervisoryDuration(SACCRTradeSpecification.tRegulatoryAssetClass assetClass, double timeToStartDate, double timeToEndDate)
    {
        if ( assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.InterestRate || assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.CreditIndex || assetClass == SACCRTradeSpecification.tRegulatoryAssetClass.CreditSingleName ) // Nr. 157: Duration: Only applied to CDS and IR
        {
            double startTime = timeToStartDate;
            double endTime = timeToEndDate ;
            double SD = Math.exp(-0.05 * startTime) - Math.exp(-0.05 * endTime );
            SD = SD / 0.05;
            return SD;
        }
        else
            return 1.0;
    }



}
