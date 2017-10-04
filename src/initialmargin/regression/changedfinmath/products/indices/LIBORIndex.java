/*
 * Created on 15.09.2006
 *
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 */
package initialmargin.regression.changedfinmath.products.indices;

import java.util.HashSet;
import java.util.Set;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;

/**
 * A forward rate index for a given period start offset (offset from fixing) and period length.
 * 
 * @author Christian Fries
 */
public class LIBORIndex extends AbstractIndex {

	private static final long serialVersionUID = 1L;
	
	private final double periodStartOffset;
	private final double periodLength;
    

	/**
	 * Creates a forward rate index for a given period start offset (offset from fixing) and period length.
	 * 
	 * @param name The name of an index. Used to map an index on a curve.
	 * @param periodStartOffset An offset added to the fixing to define the period start.
	 * @param periodLength The period length
	 */
	public LIBORIndex(String name, double periodStartOffset, double periodLength) {
		super(name);
		this.periodStartOffset = periodStartOffset;
		this.periodLength = periodLength;
	}

	/**
	 * Creates a forward rate index for a given period start offset (offset from fixing) and period length.
	 * 
	 * @param periodStartOffset An offset added to the fixing to define the period start.
	 * @param periodLength The period length
	 */
	public LIBORIndex(double periodStartOffset, double periodLength) {
		this(null, periodStartOffset, periodLength);
	}

	@Override
	public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {

		if(model.getModel().getForwardRateCurve().getName() != null && getName() != null && !model.getModel().getForwardRateCurve().getName().contains(getName())) throw new IllegalArgumentException("No curve for index " + getName() + " found in model.");
		RandomVariableInterface forwardRate = model.getLIBOR(evaluationTime, evaluationTime+periodStartOffset, evaluationTime+periodStartOffset+periodLength);
       //0.05349159667186987
		return forwardRate;
	}
	
	// INSERTED
	@Override
	public RandomVariableInterface getValue(double evaluationTime, double fixingTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {

		if(model.getModel().getForwardRateCurve().getName() != null && getName() != null && !model.getModel().getForwardRateCurve().getName().contains(getName())) throw new IllegalArgumentException("No curve for index " + getName() + " found in model.");
		RandomVariableInterface forwardRate = model.getLIBOR(Math.min(evaluationTime,fixingTime), fixingTime+periodStartOffset, fixingTime+periodStartOffset+periodLength);
//0.054631683660169905     0.027385548994700137, 0.04135428764593474, 0.05645265762282747
		return forwardRate; // 0.0272624742247487, 0.04292296554944075, 0.05854137226845246
	}

	/**
	 * Returns the periodStartOffset as an act/365 daycount.
	 * 
	 * @return the periodStartOffset
	 */
	public double getPeriodStartOffset() {
		return periodStartOffset;
	}

	/**
	 * Returns the tenor encoded as an pseudo act/365 daycount fraction.
	 * 
	 * @return the periodLength The tenor as an act/365 daycount fraction.
	 */
	public double getPeriodLength() {
		return periodLength;
	}

	@Override
	public Set<String> queryUnderlyings() {
		Set<String> underlyingNames = new HashSet<String>();
		underlyingNames.add(getName());
		return underlyingNames;
	}

	@Override
	public String toString() {
		return "LIBORIndex [periodStartOffset=" + periodStartOffset
				+ ", periodLength=" + periodLength + ", toString()="
				+ super.toString() + "]";
	}

	@Override
	public RandomVariableInterface getCF(double initialTime, double finalTime,
			LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		// TODO Auto-generated method stub
		return null;
	}
}
