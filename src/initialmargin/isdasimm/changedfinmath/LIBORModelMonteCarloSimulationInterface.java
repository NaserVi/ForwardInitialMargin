/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 01.03.2008
 */
package initialmargin.isdasimm.changedfinmath;

import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.TermStructureModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * Basic interface which has to be implemented by Monte Carlo models for LIBOR processes.
 * 
 * @author Christian Fries
 * @version 1.0
 */

public interface LIBORModelMonteCarloSimulationInterface extends TermStructureModelMonteCarloSimulationInterface {

	/**
	 * @return Returns the numberOfFactors.
	 */
	int getNumberOfFactors();

	/**
	 * Returns the libor period discretization as time discretization representing start and end dates of periods.
	 * @return Returns the libor period discretization
	 */
	TimeDiscretizationInterface getLiborPeriodDiscretization();

	/**
	 * @return The number of LIBORs in the LIBOR discretization
	 */
	int getNumberOfLibors();

	/**
	 * Returns the period start of the specified forward rate period.
	 * 
	 * @param timeIndex The index corresponding to a given time (interpretation is start of period)
	 * @return The period start of the specified forward rate period.
	 */
	double getLiborPeriod(int timeIndex);

	/**
	 * Same as java.util.Arrays.binarySearch(liborPeriodDiscretization,time).
	 * Will return a negative value if the time is not found, but then -index-1 corresponds to the index of the smallest time greater than the given one.
	 * 
	 * @param time The tenor time (fixing of the forward rate) for which the index is requested.
	 * @return The index corresponding to a given time (interpretation is start of period)
	 */
	int getLiborPeriodIndex(double time);

	/**
	 * Return the forward rate for a given simulation time index and a given forward rate index.
	 * 
	 * @param timeIndex Simulation time index.
	 * @param liborIndex Tenor time index (index corresponding to the fixing of the forward rate).
	 * @return The forward rate as a random variable.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	RandomVariableInterface getLIBOR(int timeIndex, int liborIndex) throws CalculationException;

	/**
	 * Return the forward rate curve for a given simulation time index.
	 * 
	 * @param timeIndex Simulation time index.
	 * @return The forward rate curve for a given simulation time index.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	RandomVariableInterface[] getLIBORs(int timeIndex) throws CalculationException;

	/**
	 * Returns the Brownian motion used to simulate the curve.
	 * 
	 * @return The Brownian motion used to simulate the curve.
	 */
	BrownianMotionInterface getBrownianMotion();

	/**
	 * Returns the underlying model.
	 * 
	 * The model specifies the measure, the initial value, the drift, the factor loadings (covariance model), etc.
	 * 
	 * @return The underlying model
	 */
	TermStructureModelInterface getModel();
	
	/**
	 * Clears the cache containing the numeraire adjustments
	 */
	void clearNumeraireAdjustmentCache();
	
	/**
	 * Returns the map of <code> Double <code> (time) and <code> RandomVariableInterface <code> (numeraire Adjustment)
	 * @return The numeraire adjustment map
	 */
	public Map<Double, RandomVariableInterface> getNumeraireAdjustmentMap();
	
	/**
	 * Returns the numeraire adjustment
	 * @param time The time
	 * @return
	 * @throws CalculationException 
	 */
	RandomVariableInterface getNumeraireAdjustment(double time) throws CalculationException;
	
	/** 
	 * Returns the forward bond P(T;t) on the forward curve. Calculated directly from Libors without using conditional expectation
	 * @param T final time
	 * @param t initial time 
	 * @return P(T;t)
	 * @throws CalculationException 
	 */
	RandomVariableInterface getForwardBondLibor(double T, double t) throws CalculationException;

	/**
	 * Returns the forward bond P(T;t) from on the OIS curve for a given Libor market model
	 * @param T The maturity of the forward bond 
	 * @param t The inception of the forward bond
	 * @return The forward bond P(T;t) on the OIS curve
	 * @throws CalculationException
	 */
	RandomVariableInterface getForwardBondOIS(double T, double t) throws CalculationException;

	/**
	 * Return a clone of this model with a modified Brownian motion using a different seed.
	 * 
	 * @param seed The seed
	 * @return Clone of this object, but having a different seed.
	 * @deprecated
	 */
	@Deprecated
	Object getCloneWithModifiedSeed(int seed);

}