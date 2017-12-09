package initialmargin.isdasimm.changedfinmath.modelplugins;
/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 08.08.2005
 */
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import net.finmath.montecarlo.AbstractRandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * @author Christian Fries
 */
public class LIBORVolatilityModelPiecewiseConstant extends LIBORVolatilityModel {

	private final AbstractRandomVariableFactory	randomVariableFactory;

	private final TimeDiscretizationInterface	simulationTimeDiscretization;
	private final TimeDiscretizationInterface	timeToMaturityDiscretization;

	private Map<Integer, HashMap<Integer, Integer>> 	indexMap = new HashMap<Integer, HashMap<Integer, Integer>>();
	private RandomVariableInterface[] volatility;
	private final	boolean		isCalibrateable;

	public LIBORVolatilityModelPiecewiseConstant(AbstractRandomVariableFactory randomVariableFactory, TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, TimeDiscretizationInterface simulationTimeDiscretization, TimeDiscretizationInterface timeToMaturityDiscretization, double[] volatility, boolean isCalibrateable) {
		super(timeDiscretization, liborPeriodDiscretization);

		this.randomVariableFactory = randomVariableFactory;

		/*
		 * Build index map
		 */
		double maxMaturity = timeToMaturityDiscretization.getTime(timeToMaturityDiscretization.getNumberOfTimes()-1);
		int volatilityIndex = 0;
		for(int simulationTime=0; simulationTime<simulationTimeDiscretization.getNumberOfTimes(); simulationTime++) {
			HashMap<Integer, Integer> timeToMaturityIndexing = new HashMap<Integer, Integer>();
			for(int timeToMaturity=0; timeToMaturity<timeToMaturityDiscretization.getNumberOfTimes(); timeToMaturity++) {
				if(simulationTimeDiscretization.getTime(simulationTime)+timeToMaturityDiscretization.getTime(timeToMaturity) > maxMaturity) 
					continue;

				timeToMaturityIndexing.put(timeToMaturity,volatilityIndex++);
			}
			indexMap.put(simulationTime, timeToMaturityIndexing);
		}

		// catch case where all values are the same
		if(volatility.length == 1) {
			double volatilityValue = volatility[0];
			volatility = new double[volatilityIndex];
			Arrays.fill(volatility, volatilityValue);
		}
		
		this.volatility = new RandomVariableInterface[volatilityIndex];
		setParameter(volatility, true);

		if(volatilityIndex != this.volatility.length) throw new IllegalArgumentException("volatility.length should equal simulationTimeDiscretization.getNumberOfTimes()*timeToMaturityDiscretization.getNumberOfTimes().");
		this.simulationTimeDiscretization = simulationTimeDiscretization;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.isCalibrateable = isCalibrateable;
	}

	public LIBORVolatilityModelPiecewiseConstant(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, TimeDiscretizationInterface simulationTimeDiscretization, TimeDiscretizationInterface timeToMaturityDiscretization, double[] volatility, boolean isCalibrateable) {
		this(new RandomVariableFactory(), timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, volatility, isCalibrateable);
	}

	public LIBORVolatilityModelPiecewiseConstant(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, TimeDiscretizationInterface simulationTimeDiscretization, TimeDiscretizationInterface timeToMaturityDiscretization, double volatility, boolean isCalibrateable) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, new double[] { volatility }, isCalibrateable);
	}

	public LIBORVolatilityModelPiecewiseConstant(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, TimeDiscretizationInterface simulationTimeDiscretization, TimeDiscretizationInterface timeToMaturityDiscretization, double[] volatility) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, volatility, true);
	}

	public LIBORVolatilityModelPiecewiseConstant(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, TimeDiscretizationInterface simulationTimeDiscretization, TimeDiscretizationInterface timeToMaturityDiscretization, double volatility) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, new double[] { volatility });
	}

	@Override
	public RandomVariableInterface[] getParameterAsRandomVariable() {
		if(!isCalibrateable) return null;
		
		return volatility;
	}

	@Override
	public void setParameter(double[] parameter) {
		setParameter(parameter, isCalibrateable);
	}

	private void setParameter(double[] parameter, boolean executeFunction) {
		if(executeFunction) {
			int numberOfParameters = parameter.length;
			if(numberOfParameters != volatility.length) throw new IllegalArgumentException("parameter.length has to conincide with volatility.length!");
			for(int i = 0; i < parameter.length; i++)
				this.volatility[i] = randomVariableFactory.createRandomVariable(parameter[i]);
			}
	}
	
	@Override
	public RandomVariableInterface getVolatility(int timeIndex, int liborIndex) {
		// Create a very simple volatility model here
		double time             = getTimeDiscretization().getTime(timeIndex);
		double maturity         = getLiborPeriodDiscretization().getTime(liborIndex);
		double timeToMaturity   = maturity-time;


		// initialize with zero and time 
		RandomVariableInterface volatilityInstanteaneous = randomVariableFactory.createRandomVariable(time, 0.0); 
		
		if(timeToMaturity > 0)
		{
			int timeIndexSimulationTime = simulationTimeDiscretization.getTimeIndex(time);
			if(timeIndexSimulationTime < 0) timeIndexSimulationTime = -timeIndexSimulationTime-1-1;
			if(timeIndexSimulationTime < 0) timeIndexSimulationTime = 0;
			if(timeIndexSimulationTime >= simulationTimeDiscretization.getNumberOfTimes()) timeIndexSimulationTime--;

			int timeIndexTimeToMaturity = timeToMaturityDiscretization.getTimeIndex(timeToMaturity);
			if(timeIndexTimeToMaturity < 0) timeIndexTimeToMaturity = -timeIndexTimeToMaturity-1-1;
			if(timeIndexTimeToMaturity < 0) timeIndexTimeToMaturity = 0;
			if(timeIndexTimeToMaturity >= timeToMaturityDiscretization.getNumberOfTimes()) timeIndexTimeToMaturity--;

			volatilityInstanteaneous = volatilityInstanteaneous.add(volatility[indexMap.get(timeIndexSimulationTime).get(timeIndexTimeToMaturity)]);
		}
	
		return volatilityInstanteaneous;
	}

	@Override
	public Object clone() {
		return new LIBORVolatilityModelPiecewiseConstant(
				randomVariableFactory,
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				this.simulationTimeDiscretization,
				this.timeToMaturityDiscretization,
				this.getParameter(),
				this.isCalibrateable
				);
	}

}