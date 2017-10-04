/*
 * Created on 22.11.2009
 */
package initialmargin.regression.changedfinmath.products.components;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.MonteCarloSimulationInterface;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.products.indices.LIBORIndex;
import net.finmath.stochastic.RandomVariableInterface;

/**
 * A period. A period has references to the index (coupon) and the notional.
 * It provides the fixing date for the index, the period length, and the payment date.
 * 
 * @author Christian Fries
 * @version 1.1
 */
public class Period extends AbstractPeriod {

	private static final long serialVersionUID = -7107623461781510475L;
	private final boolean couponFlow;
	private final boolean notionalFlow;
	private final boolean payer;
	private final boolean isExcludeAccruedInterest;

	/**
	 * Create a simple period with notional and index (coupon) flow.
	 * 
	 * @param periodStart The period start.
	 * @param periodEnd The period end.
	 * @param fixingDate The fixing date (as double).
	 * @param paymentDate The payment date (as double).
	 * @param notional The notional object relevant for this period.
	 * @param index The index (used for coupon calculation) associated with this period.
	 * @param daycountFraction The daycount fraction (<code>coupon = index(fixingDate) * daycountFraction</code>).
	 * @param couponFlow If true, the coupon will be payed. Otherwise there will be not coupon flow.
	 * @param notionalFlow If true, there will be a positive notional flow at period start (but only if peirodStart &gt; evaluationTime) and a negative notional flow at period end (but only if periodEnd &gt; evaluationTime). Otherwise there will be no notional flows.
	 * @param payer If true, the period will be a payer period, i.e. notional and coupon at period end are payed (negative). Otherwise it is a receiver period.
	 * @param isExcludeAccruedInterest If the true, the valuation will exclude accrued interest, if any.
	 */
	public Period(double periodStart, double periodEnd, double fixingDate,
			double paymentDate, AbstractNotional notional, AbstractProductComponent index, double daycountFraction,
			boolean couponFlow, boolean notionalFlow, boolean payer, boolean isExcludeAccruedInterest) {
		super(periodStart, periodEnd, fixingDate, paymentDate, notional, index, daycountFraction);
		this.couponFlow = couponFlow;
		this.notionalFlow = notionalFlow;
		this.payer = payer;
		this.isExcludeAccruedInterest = isExcludeAccruedInterest;
	}

	/**
	 * Create a simple period with notional and index (coupon) flow.
	 * 
	 * The valuation does not exclude the accrued interest, i.e., the valuation reports a so called dirty price.
	 * 
	 * @param periodStart The period start.
	 * @param periodEnd The period end.
	 * @param fixingDate The fixing date (as double).
	 * @param paymentDate The payment date (as double).
	 * @param notional The notional object relevant for this period.
	 * @param index The index (used for coupon calculation) associated with this period.
	 * @param daycountFraction The daycount fraction (<code>coupon = index(fixingDate) * daycountFraction</code>).
	 * @param couponFlow If true, the coupon will be payed. Otherwise there will be not coupon flow.
	 * @param notionalFlow If true, there will be a positive notional flow at period start (but only if peirodStart &gt; evaluationTime) and a negative notional flow at period end (but only if periodEnd &gt; evaluationTime). Otherwise there will be no notional flows.
	 * @param payer If true, the period will be a payer period, i.e. notional and coupon at period end are payed (negative). Otherwise it is a receiver period.
	 */
	public Period(double periodStart, double periodEnd, double fixingDate,
			double paymentDate, AbstractNotional notional, AbstractProductComponent index, double daycountFraction,
			boolean couponFlow, boolean notionalFlow, boolean payer) {
		this(periodStart, periodEnd, fixingDate, paymentDate, notional, index, daycountFraction, couponFlow, notionalFlow, payer, false);
	}

	/**
	 * Create a simple period with notional and index (coupon) flow.
	 * 
	 * The valuation does not exclude the accrued interest, i.e., the valuation reports a so called dirty price.
	 * 
	 * @param periodStart The period start.
	 * @param periodEnd The period end.
	 * @param fixingDate The fixing date (as double).
	 * @param paymentDate The payment date (as double).
	 * @param notional The notional object relevant for this period.
	 * @param index The index (coupon) associated with this period.
	 * @param couponFlow If true, the coupon will be payed. Otherwise there will be not coupon flow.
	 * @param notionalFlow If true, there will be a positive notional flow at period start (but only if peirodStart &gt; evaluationTime) and a negative notional flow at period end (but only if periodEnd &gt; evaluationTime). Otherwise there will be no notional flows.
	 * @param payer If true, the period will be a payer period, i.e. notional and coupon at period end are payed (negative). Otherwise it is a receiver period.
	 */
	public Period(double periodStart, double periodEnd, double fixingDate,
			double paymentDate, AbstractNotional notional, AbstractProductComponent index,
			boolean couponFlow, boolean notionalFlow, boolean payer) {
		this(periodStart, periodEnd, fixingDate, paymentDate, notional, index, periodEnd-periodStart, couponFlow, notionalFlow, payer);
	}
	/**
	 * This method returns the value random variable of the product within the specified model, evaluated at a given evalutationTime.
	 * Note: For a lattice this is often the value conditional to evalutationTime, for a Monte-Carlo simulation this is the (sum of) value discounted to evaluation time.
	 * Cashflows prior evaluationTime are not considered.
	 * 
	 * @param evaluationTime The time on which this products value should be observed.
	 * @param model The model used to price the product.
	 * @return The random variable representing the value of the product discounted to evaluation time
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method. 
	 */
	@Override
	public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {        

		if(evaluationTime >= this.getPaymentDate()) return new RandomVariable(0.0);

		// Get random variables
		RandomVariableInterface	notionalAtPeriodStart	= getNotional().getNotionalAtPeriodStart(this, model);
		RandomVariableInterface	numeraireAtEval			= model.getNumeraire(evaluationTime);   //1.1301009311819583, 1.2422123463369106, 1.161488600768824, 1.2356093606697227, 1.302681217720419, 1.5392075870444337, 1.3977779276528706, 1.0590818316766346
		RandomVariableInterface	numeraire				= model.getNumeraire(getPaymentDate()); //1.1306018341640918, 1.2428440472798723, 1.1628577880118596, 1.235946759351182, 1.3039588907242465, 1.5450215130402047, 1.39840277997152, 1.0588383837282285, 1.1038289696165438, 1.2473983074952235, 1.0756038020344902, 1.1868661380408725, 1.1177126906537413, 1.2432190229131752, 1.0855080139683173, 1.227599930869064, 1.1338831159673621, 1.1440739948317777
		// @TODO: Add support for weighted Monte-Carlo.
		//        RandomVariableInterface	monteCarloProbabilities	= model.getMonteCarloWeights(getPaymentDate());

		RandomVariableInterface values;

		// Calculate numeraire relative value of coupon flows
		if(couponFlow) {
			values = getCoupon(evaluationTime, model); //0.027467596506918757
			values = values.mult(notionalAtPeriodStart);
			values = values.div(numeraire);
			if(isExcludeAccruedInterest && evaluationTime >= getPeriodStart() && evaluationTime < getPeriodEnd()) {
				double nonAccruedInterestRatio = (getPeriodEnd() - evaluationTime) / (getPeriodEnd() - getPeriodStart());
				values = values.mult(nonAccruedInterestRatio);
			}
		}
		else {
			values = new RandomVariable(0.0,0.0);
		}

		// Apply notional exchange
		if(notionalFlow) {
			RandomVariableInterface	nationalAtPeriodEnd		= getNotional().getNotionalAtPeriodEnd(this, model);

			if(getPeriodStart() > evaluationTime) {
				RandomVariableInterface	numeraireAtPeriodStart	= model.getNumeraire(getPeriodStart());
				values = values.subRatio(notionalAtPeriodStart, numeraireAtPeriodStart);
			}

			if(getPeriodEnd() > evaluationTime) {
				RandomVariableInterface	numeraireAtPeriodEnd	= model.getNumeraire(getPeriodEnd());
				values = values.addRatio(nationalAtPeriodEnd, numeraireAtPeriodEnd);
			}
		}

		if(payer) values = values.mult(-1.0);

		values = values.mult(numeraireAtEval);

		// Return values
		return values;	
	}
	
	//INSERTED
	@Override
	public RandomVariableInterface getValue(double evaluationTime, double fixingTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {        

		if(evaluationTime >= this.getPaymentDate()) return new RandomVariable(0.0);

		// Get random variables
		RandomVariableInterface	notionalAtPeriodStart	= getNotional().getNotionalAtPeriodStart(this, model);
		RandomVariableInterface	numeraireAtEval			= model.getNumeraire(evaluationTime);
		RandomVariableInterface	numeraire				= model.getNumeraire(getPaymentDate());
		// @TODO: Add support for weighted Monte-Carlo.
		//        RandomVariableInterface	monteCarloProbabilities	= model.getMonteCarloWeights(getPaymentDate());

		RandomVariableInterface values;

		// Calculate numeraire relative value of coupon flows
		if(couponFlow) {
			values = getCoupon(evaluationTime, model);
			values = values.mult(notionalAtPeriodStart);
			values = values.div(numeraire);
			if(isExcludeAccruedInterest && evaluationTime >= getPeriodStart() && evaluationTime < getPeriodEnd()) {
				double nonAccruedInterestRatio = (getPeriodEnd() - evaluationTime) / (getPeriodEnd() - getPeriodStart());
				values = values.mult(nonAccruedInterestRatio);
			}
		}
		else {
			values = new RandomVariable(0.0,0.0);
		}

		// Apply notional exchange
		if(notionalFlow) {
			RandomVariableInterface	nationalAtPeriodEnd		= getNotional().getNotionalAtPeriodEnd(this, model);

			if(getPeriodStart() > evaluationTime) {
				RandomVariableInterface	numeraireAtPeriodStart	= model.getNumeraire(getPeriodStart());
				values = values.subRatio(notionalAtPeriodStart, numeraireAtPeriodStart);
			}

			if(getPeriodEnd() > evaluationTime) {
				RandomVariableInterface	numeraireAtPeriodEnd	= model.getNumeraire(getPeriodEnd());
				values = values.addRatio(nationalAtPeriodEnd, numeraireAtPeriodEnd);
			}
		}

		if(payer) values = values.mult(-1.0);

		values = values.mult(numeraireAtEval);

		// Return values
		return values;	
	}
	
	
	

	
	public RandomVariableInterface getCoupon(LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		// Calculate percentage value of coupon (not multiplied with notional, not discounted)
		RandomVariableInterface values = getIndex().getValue(getFixingDate(), model);

		// Apply daycount fraction
		double periodDaycountFraction = getDaycountFraction();
		values = values.mult(periodDaycountFraction);

		return values; // hitherto not discounted
	}
	
	//INSERTED
	public RandomVariableInterface getCoupon(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		// Calculate percentage value of coupon (not multiplied with notional, not discounted)
		RandomVariableInterface values = getIndex().getValue(evaluationTime, getFixingDate(), model);

		// Apply daycount fraction
		double periodDaycountFraction = getDaycountFraction(); //0.5027777777777778
		values = values.mult(periodDaycountFraction);

		return values; // hitherto not discounted
	}

	@Override
	public String toString() {
		return "Period [couponFlow=" + couponFlow + ", notionalFlow="
				+ notionalFlow + ", payer=" + payer + ", toString()="
				+ super.toString() + "]";
	}

	@Override
	public RandomVariableInterface getCF(double initialTime, double finalTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		RandomVariableInterface values = new RandomVariable(0.0);
		if(initialTime >= this.getPaymentDate() || finalTime <this.getPaymentDate()) {
		   // Apply notional exchange
		   if(notionalFlow && finalTime>=getPeriodEnd()) {
			   RandomVariableInterface	notionalAtPeriodEnd		= getNotional().getNotionalAtPeriodEnd(this, model);
			   RandomVariableInterface	numeraireAtPeriodEnd	= model.getNumeraire(getPeriodEnd());
			   values = values.addRatio(notionalAtPeriodEnd, numeraireAtPeriodEnd);	
	    	} else return values;
		} else {
		
		// Get random variables
		RandomVariableInterface	notionalAtPeriodStart	= getNotional().getNotionalAtPeriodStart(this, model);
		RandomVariableInterface	numeraire				= model.getNumeraire(getPaymentDate());
		// @TODO: Add support for weighted Monte-Carlo.
		//        RandomVariableInterface	monteCarloProbabilities	= model.getMonteCarloWeights(getPaymentDate());

		
		// Calculate numeraire relative value of coupon flows
		if(couponFlow) {
			values = getCoupon(finalTime, model); //not discounted
			values = values.mult(notionalAtPeriodStart);
			values = values.div(numeraire);
			if(isExcludeAccruedInterest && finalTime >= getPeriodStart() && finalTime < getPeriodEnd()) {
				double nonAccruedInterestRatio = (getPeriodEnd() - finalTime) / (getPeriodEnd() - getPeriodStart());
				values = values.mult(nonAccruedInterestRatio);
			}
		}
		else {
			values = new RandomVariable(0.0,0.0);
		}

		// Apply notional exchange
		if(notionalFlow && finalTime>=getPeriodEnd()) {
			RandomVariableInterface	notionalAtPeriodEnd		= getNotional().getNotionalAtPeriodEnd(this, model);
			RandomVariableInterface	numeraireAtPeriodEnd	= model.getNumeraire(getPeriodEnd());
			values = values.addRatio(notionalAtPeriodEnd, numeraireAtPeriodEnd);
			
		}
		}

		if(payer) values = values.mult(-1.0);
		RandomVariableInterface	numeraireAtEval			= model.getNumeraire(finalTime); // was finalTime
		values = values.mult(numeraireAtEval);

		// Return values
		return values;	
	}
	
}
