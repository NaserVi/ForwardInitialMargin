package initialmargin.isdasimm.sensitivity;

import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import initialmargin.isdasimm.changedfinmath.LIBORModelMonteCarloSimulationInterface;
import initialmargin.isdasimm.products.AbstractSIMMProduct;
import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariableInterface;

/** This class contains some functions and methods which we need to calculate forward initial margin.
 * 
 * @author Mario Viehmann
 *
 */
public abstract class AbstractSIMMSensitivityCalculation {
 
    public boolean isUseTimeGridAdjustment; 
    public boolean isUseAnalyticSwapSensitivities;
    public boolean isConsiderOISSensitivities;
    
	public enum SensitivityMode{
    	LinearMelting,
    	Interpolation,
    	Exact, // AAD or Analytic (for Swaps)
    	InterpolationOIS
    }
    
    protected SensitivityMode sensitivityMode;
    
    public enum WeightMode{
		Constant,  //Sets dL/dS(t=0) for all forward IM times, i.e. leave the weight adjustment dL/dS constant
		Stochastic //Calculate dL/dS(t) for all forward IM times, i.e. (weakly) stochastic weight adjustment 
	}
    
    private WeightMode liborWeightMethod;
    private HashMap<Double /*time*/, RandomVariableInterface[][]> riskWeightMap = new HashMap<>();
    
    
    /**
     * 
     * @param sensitivityMode
     * @param liborWeightMode
     * @param isUseTimeGridAdjustment
     * @param isUseAnalyticSwapSensitivities
     * @param isConsiderOISSensitivities
     */
    public AbstractSIMMSensitivityCalculation(SensitivityMode sensitivityMode, WeightMode liborWeightMode, 
    		        boolean isUseTimeGridAdjustment, boolean isUseAnalyticSwapSensitivities, boolean isConsiderOISSensitivities){
    	this.sensitivityMode = sensitivityMode;
    	this.liborWeightMethod = liborWeightMode;
    	this.isUseTimeGridAdjustment = isUseTimeGridAdjustment;
    	this.isUseAnalyticSwapSensitivities = isUseAnalyticSwapSensitivities;
    	this.isConsiderOISSensitivities = isConsiderOISSensitivities;   	
    }
    
    /**
     * 
     * @param sensitivityMode
     * @param liborWeightMode
     * @param interpolationStep
     */
    public AbstractSIMMSensitivityCalculation(SensitivityMode sensitivityMode, WeightMode liborWeightMode, double interpolationStep){
        this(sensitivityMode, liborWeightMode, false, true, true);
    }
    
    
    /** Calculate the delta SIMM sensitivities for a given risk class and index curve at a given evaluation time with the specified Libor market model.
	 *  The sensitivities are calculated by interpolation or melting (a particular case of interpolation).
	 * 		
	 * @param product The product 
	 * @param riskClass The risk class of the product 
	 * @param curveIndexName The name of the index curve
	 * @param evaluationTime The time at which the sensitivities should be calculated
	 * @param model The Libor market model
	 * @return The sensitivities at evaluationTime on the SIMM buckets (or on Libor Buckets)
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	public abstract RandomVariableInterface[] getDeltaSensitivitiesIR(AbstractSIMMProduct product,
															 		  String riskClass, 
															 		  String curveIndexName,
															 		  double evaluationTime, 
															 		  LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException;
		
		
	/**
	 * 
	 * @param product The product
	 * @param riskClass The risk class of the product
	 * @param curveIndexName The name of the curve
	 * @param evaluationTime The time at which the sensitivities should be calculated
	 * @param model The Libor market model
	 * @return The delta sensitivities calculated by AAD
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	public abstract RandomVariableInterface[] getExactDeltaSensitivitiesIR(AbstractSIMMProduct product,
			 															   String riskClass, 
			 															   String curveIndexName,
			 															   double evaluationTime, 
			 															   LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException;
    
	/** Get the sensitivities using sensitivity melting on SIMM Buckets. 
	 * 
	 * @param product
	 * @param sensitivities
	 * @param meltingZeroTime
	 * @param evaluationTime
	 * @param curveIndexName
	 * @param riskClass
	 * @return
	 * @throws SolverException
	 * @throws CloneNotSupportedException
	 * @throws CalculationException
	 */
	public abstract RandomVariableInterface[] getMeltedSensitivities(AbstractSIMMProduct product, RandomVariableInterface[] sensitivities, double meltingZeroTime,
			double evaluationTime, String curveIndexName, String riskClass) throws SolverException, CloneNotSupportedException, CalculationException;
	
	
	/**Calculate the sensitivities dV/dS with respect to all swap rates for given product and curve. This applies to the risk class Interest Rates only.
 	 * 
 	 * @param product The SIMM product
 	 * @param curveIndexName The name of the curve to be considered (OIS, LiborXm)
 	 * @param evaluationTime The time at which the initial margin is calculated
 	 * @param model The Libor market model
 	 * @return The sensitivities dV/dS i.e. with respect to swap rates.
 	 * @throws SolverException
 	 * @throws CloneNotSupportedException
 	 * @throws CalculationException
 	 * 
 	 */
 	public RandomVariableInterface[] doCalculateDeltaSensitivitiesIR(AbstractSIMMProduct product,
 																	 String curveIndexName, 
                                                                     double evaluationTime,
                                                                     LIBORModelMonteCarloSimulationInterface model) throws SolverException, CloneNotSupportedException, CalculationException{
 		
 		RandomVariableInterface[] delta = null; // The vector of delta sensitivities on all SIMM Buckets
 		
 		switch(curveIndexName){
 		    
 		   case("Libor6m"):
 			   
 		        RandomVariableInterface[] dVdL = product.getValueLiborSensitivities(evaluationTime, model);
 		  
 		        // Calculate dV/dS = dV/dL * dL/dS
 		        delta = getValueSwapSensitivities(evaluationTime, dVdL, model);  
 		       
 		        // Map Sensitivities on SIMM Buckets
 		        delta = mapSensitivitiesOnBuckets(delta, "InterestRate" /*riskClass*/, null, model);		       
 		        
 		        break;
 		
 		   case("OIS"):
 			   
 		        if(isConsiderOISSensitivities) {
 			    
 		          // Calculate dV/dS = dV/dP * dP/dS. These Sensis are already on SIMM Buckets
 		          delta = product.getDiscountCurveSensitivities("InterestRate" /*riskClass*/,evaluationTime, model); 
 			    		    
 		        } else {
 		          // Set OIS sensitivities to zero on all SIMM buckets
 		          delta = new RandomVariableInterface[12]; Arrays.fill(delta, new RandomVariable(0.0));
 		        }
 		   
 		        break;
 		        
 		    default:
 		    		throw new IllegalArgumentException("Unknow curve: " + curveIndexName);
 		}
 	
 		return delta;
 	}
 	
 	
 	/**Calculate the sensitivities of the value of a product w.r.t. swap rates given the Libor sensitivities dV/dL
	 * 
	 * @param evaluationTime The time of evaluation
	 * @param dVdL The vector of derivatives dV/dL = dV/dL_0,...,dV/dL_n
	 * @param model The Libor Market Model
	 * @return The derivatives dV/dS 
	 * @throws CalculationException
	 */
	public RandomVariableInterface[] getValueSwapSensitivities(double evaluationTime, 
																RandomVariableInterface[] dVdL,
																LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		// the following line will be removed later. Just checking how timeGridAdjustment affects the result
		int timeGridIndicator = 0; if(!isUseTimeGridAdjustment && !onLiborPeriodDiscretization(evaluationTime,model)) timeGridIndicator = 1;
		
		RandomVariableInterface[] delta = new RandomVariableInterface[dVdL.length-timeGridIndicator];
		RandomVariableInterface[][] dLdS;
		if(this.liborWeightMethod == WeightMode.Stochastic){
			   dLdS = getLiborSwapSensitivities(evaluationTime, model);
		} else dLdS = getLiborSwapSensitivities(0.0, model);
		// Calculate Sensitivities wrt Swaps
		// return multiply(dVdL,dLdS);
		for(int swapIndex = 0; swapIndex<dVdL.length-timeGridIndicator; swapIndex++){
			RandomVariableInterface dVdS  =new RandomVariable(0.0);
			RandomVariableInterface factor;
			for(int liborIndex=0;liborIndex<dVdL.length-timeGridIndicator;liborIndex++){
			    factor = dLdS[liborIndex][swapIndex]==null ?  new RandomVariable(0.0) : dLdS[liborIndex][swapIndex];
			    dVdS = dVdS.addProduct(dVdL[liborIndex+timeGridIndicator], factor);
		    }
			delta[swapIndex]=dVdS;
		}
		return delta;
	 }
	
	/**Performs rebucketing of sensitivities to the SIMM buckets by linear interpolation (Source: Master Thesis of Jamal Issa, modified).
	 * 
	 * @param sensitivities The sensitivities wrt swap rates dV/dS
	 * @param riskClass The risk class
	 * @param riskFactorDays The number of days corresponding to the sensitivities
	 * @return The sensitivities on the SIMM maturity buckets
	 */
	public static RandomVariableInterface[] mapSensitivitiesOnBuckets(RandomVariableInterface[] sensitivities, String riskClass, int[] riskFactorDays, LIBORModelMonteCarloSimulationInterface model){
		//rebucketing to SIMM structure(buckets: 2w, 1m, 3m, 6m, 1y, 2y, 3y, 5y, 10y, 15y, 20y, 30y)	
		int[] riskFactorsSIMM = riskClass=="InterestRate" ? new int[] {14, 30, 90, 180, 365, 730, 1095, 1825, 3650, 5475, 7300, 10950} : /*Credit*/ new int[] {365, 730, 1095, 1825, 3650};	
		RandomVariableInterface[] deltaSIMM = new RandomVariableInterface[riskFactorsSIMM.length];
		for(int i = 0;i<deltaSIMM.length;i++) deltaSIMM[i] = new RandomVariable(0.0);
		
		if(riskFactorDays==null) riskFactorDays = riskFactorDaysLibor(sensitivities, model);
		int counter = 0;
		for(int simmFactor =0; simmFactor<riskFactorsSIMM.length;simmFactor++){
			for(int i = counter; i<sensitivities.length; i++){
				
							    
					if(riskFactorDays[i] < riskFactorsSIMM[0]){
						deltaSIMM[0] = deltaSIMM[0].add(sensitivities[i]);
						counter++;
					}
					else{
						if(riskFactorDays[i] >= riskFactorsSIMM[riskFactorsSIMM.length-1]){
							deltaSIMM[deltaSIMM.length-1] = deltaSIMM[deltaSIMM.length-1].add(sensitivities[i]);
						}
					
						else{
							if(riskFactorDays[i] >= riskFactorsSIMM[simmFactor] && riskFactorDays[i] < riskFactorsSIMM[simmFactor+1]){
					
							deltaSIMM[simmFactor] = deltaSIMM[simmFactor].addProduct(sensitivities[i],((double)(riskFactorsSIMM[simmFactor+1] - riskFactorDays[i]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							deltaSIMM[simmFactor+1] = deltaSIMM[simmFactor+1].addProduct(sensitivities[i],((double)(riskFactorDays[i]-riskFactorsSIMM[simmFactor]) / (riskFactorsSIMM[simmFactor+1]-riskFactorsSIMM[simmFactor])));
							counter++;
							}							
							else{
							break;
							}
						}
					
					}
			}
			
		}
		
	return deltaSIMM;		
			
	}
	
	
	/**Calculates dL/dS: The risk weights used to apply Libor Sensitivities to SIMM: dV/dS = dV/dL * dL/dS.
	 * 
	 * @param evaluationTime The time at which the sensitivity is calculated
	 * @param model The Libor market model
	 * @return The matrix dL/dS 
	 * @throws CalculationException
	 */
	private RandomVariableInterface[][] getLiborSwapSensitivities(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{

		if(riskWeightMap.containsKey(evaluationTime)) return riskWeightMap.get(evaluationTime);
		
		RandomVariableInterface[][] dLdS=null;
		double liborPeriodLength = model.getLiborPeriodDiscretization().getTimeStep(0);
		
	    // Get index of first Libor starting >= evaluationTime
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int numberOfRemainingLibors = model.getNumberOfLibors()-nextLiborIndex;
		dLdS = new RandomVariableInterface [numberOfRemainingLibors][numberOfRemainingLibors];
					
		// Calculate dLdS directly  
		dLdS[0][0]=new RandomVariable(1.0);
		double discountTime = evaluationTime+liborPeriodLength;
		RandomVariableInterface sumDf = model.getForwardBondOIS(discountTime, evaluationTime);
		for(int liborIndex = 1; liborIndex<dLdS.length;liborIndex++){
		    discountTime +=liborPeriodLength;
		    RandomVariableInterface denominator = model.getForwardBondOIS(discountTime, evaluationTime);
		    dLdS[liborIndex][liborIndex-1]=sumDf.div(denominator).mult(-1.0);
		    sumDf = sumDf.add(denominator);
		    dLdS[liborIndex][liborIndex] = sumDf.div(denominator);
		  
		}
		
		riskWeightMap.put(evaluationTime, dLdS);		
		return dLdS;
	}
	
	
	/**Calculates dPdS in a single curve context. Used for calculating sensis with respect to discount curve.
	 * 
	 * @param evaluationTime The time at which the initial margin is calculated
	 * @param model The Libor market model
	 * @param timeSteps The time steps of the swap (between payment dates)
	 * @return The sensitivity of the discount curve (bonds) wrt to swap rates of the same curve.
	 * @throws CalculationException 
	 */
	public static RandomVariableInterface[][] getBondSwapSensitivity(double evaluationTime, 
															   LIBORModelMonteCarloSimulationInterface model, 
															   double[] timeSteps) throws CalculationException{
		int numberOfBonds = timeSteps.length;
		double bondTime = evaluationTime;
		RandomVariableInterface sum= new RandomVariable(0.0);
		RandomVariableInterface[][] dSdP = new RandomVariableInterface[numberOfBonds][numberOfBonds];
		
		for(int swapIndex=0;swapIndex<dSdP.length;swapIndex++){
			bondTime += timeSteps[swapIndex];
			RandomVariableInterface bondOIS = model.getForwardBondOIS(bondTime /*T*/,evaluationTime /*t*/); //P^OIS(T;t)
		    sum = sum.addProduct(bondOIS,timeSteps[swapIndex]);
		    for(int bondIndex=0;bondIndex<dSdP.length;bondIndex++){
		    	if(swapIndex<bondIndex) dSdP[swapIndex][bondIndex] = new RandomVariable(0.0);
		    	else if(swapIndex==bondIndex) dSdP[swapIndex][bondIndex] = sum.mult(-1.0).sub(timeSteps[swapIndex]).addProduct(bondOIS,timeSteps[swapIndex]).div(sum.squared());
		    	else dSdP[swapIndex][bondIndex] = bondOIS.sub(1.0).mult(timeSteps[bondIndex]).div(sum.squared());    			    	
		    }
		} 

		return getPseudoInverse(dSdP, model.getNumberOfPaths()); // PseudoInverse == Inverse for n x n matrix.
	}
	
	
	/**Since dV/dL is wrt the incorrect Libor times this function provides a matrix dL/dL to be multiplied with dV/dL in order to 
	 * have the correct libor times starting at evaluationTime. 
	 * @param evaluationTime The time at which the adjustment should be calculated.
	 * @param model The Libor market model
	 * @return Pseudo Inverse of derivative band matrix; Identity matrix in case of evaluationTime on LiborPeriodDiscretization; 
	 * @throws CalculationException
	 */
	public static RandomVariableInterface[][] getLiborTimeGridAdjustment(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException{
		int numberOfRemainingLibors = getNumberOfRemainingLibors(evaluationTime,model);
		
		// If evaluationTime lies on Libor Time Grid - return identity matrix
		if (onLiborPeriodDiscretization(evaluationTime,model)) {
			RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors];
			for(int i=0;i<dLdL.length;i++) dLdL[i][i]=new RandomVariable(1.0);
		    return dLdL;
		}
		
		// Calculate dLdL. It is a (n-1)x n Matrix!
		RandomVariableInterface[][] dLdL = new RandomVariableInterface[numberOfRemainingLibors][numberOfRemainingLibors+1];
		double swapTenorLength = model.getLiborPeriodDiscretization().getTimeStep(0); // Model must have same tenor as swap!
		double timeOfFirstLiborPriorToEval = getPreviousLiborTime(evaluationTime,model);
		int timeIndexAtEvaluationTime = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		int timeIndexAtFirstLiborPriorToEval = model.getTimeDiscretization().getTimeIndexNearestGreaterOrEqual(timeOfFirstLiborPriorToEval);
		
		for(int liborIndex = 0; liborIndex <numberOfRemainingLibors; liborIndex++){
			double liborTime = evaluationTime+liborIndex*swapTenorLength; // t+j*\Delta T
		    int    previousLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestLessOrEqual(liborTime);
		    double previousLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex);
		    double firstNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+1);
		    double secondNextLiborTime = model.getLiborPeriodDiscretization().getTime(previousLiborIndex+2);
		    double factor1 = (secondNextLiborTime-(liborTime+swapTenorLength))/(secondNextLiborTime-firstNextLiborTime);
		    double factor2 = (liborTime-previousLiborTime)/(firstNextLiborTime-previousLiborTime);
		    int    timeIndex = liborIndex==0 ? timeIndexAtFirstLiborPriorToEval : timeIndexAtEvaluationTime;
		    // Get Libors
		    RandomVariableInterface previousLibor = model.getLIBOR(timeIndex, previousLiborIndex);     
		    RandomVariableInterface nextLibor     = model.getLIBOR(timeIndex, previousLiborIndex + 1); 
		    RandomVariableInterface logInterpol = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).log().mult(-factor1);
		                            logInterpol = logInterpol.add(previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).log().mult(-factor2)).exp();
		    // Set derivatives
		    dLdL[liborIndex][liborIndex]   = nextLibor.mult(secondNextLiborTime-firstNextLiborTime).add(1.0).mult(logInterpol).mult(1-factor2);// dLdL_i-1
		    dLdL[liborIndex][liborIndex+1] = previousLibor.mult(firstNextLiborTime-previousLiborTime).add(1.0).mult(logInterpol).mult(1-factor1);
		}
		
		// dLdL is (n-1) x n matrix. Get PseudoInverse for all paths and then put it back together as RV
		return getPseudoInverse(dLdL, model.getNumberOfPaths());
	}
	
	public static int[] riskFactorDaysLibor(RandomVariableInterface[] sensis, LIBORModelMonteCarloSimulationInterface model){
		
	    int[] riskFactorDays = new int[sensis.length];
		// act/365 as default daycount convention
		for(int i=0;i<sensis.length;i++) riskFactorDays[i] = (int)Math.round(365 * model.getLiborPeriodDiscretization().getTime(i+1));	
		
		return riskFactorDays;
	}
    
   
	
	public static boolean onLiborPeriodDiscretization(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		return (evaluationTime == getNextLiborTime(evaluationTime,model));
	}
	
	public static int getNumberOfRemainingLibors(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getNumberOfLibors()-nextLiborIndex;
	}
	
	public static double getNextLiborTime(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex);
	}
	
	public static double getPreviousLiborTime(double evaluationTime, LIBORModelMonteCarloSimulationInterface model){
		if(evaluationTime==0) return 0.0;
		int nextLiborIndex = model.getLiborPeriodDiscretization().getTimeIndexNearestGreaterOrEqual(evaluationTime);
		return model.getLiborPeriodDiscretization().getTime(nextLiborIndex-1);
	}
	
	
	public SensitivityMode getSensitivityMode(){
		return this.sensitivityMode;
	}
	
	public WeightMode getWeightMode(){
		return this.liborWeightMethod;
	}
	
	public void clearRiskWeights(){
		if(riskWeightMap!=null) riskWeightMap.clear();
	}
	
	public void setWeightMode(WeightMode mode){
		this.liborWeightMethod = mode;
	}
	
	public void setSensitivityMode(SensitivityMode mode){
		this.sensitivityMode = mode;
	}
	//----------------------------------------------------------------------------------------------------------------------------------
    // Some auxiliary functions
	//----------------------------------------------------------------------------------------------------------------------------------

		/**Calculate Pseudo Inverse of matrix of type RandomVariableInterface[][]
		 * 
		 * @param matrix The matrix for which the pseudo inverse is calculated
		 * @return The pseudo inverse of the matrix
		 */
	    public static RandomVariableInterface[][] getPseudoInverse(RandomVariableInterface[][] matrix, int numberOfPaths){
	    	double[][][] inv = new double[matrix[0].length][matrix.length][numberOfPaths];
			double[][] matrixOnPath = new double[matrix.length][matrix[0].length];
			for(int pathIndex=0; pathIndex<numberOfPaths; pathIndex++){
				// Get double[][] matrix on path
				for(int i=0;i<matrixOnPath.length;i++){
					for(int j=0;j<matrixOnPath[0].length;j++){
						matrixOnPath[i][j]=matrix[i][j]==null ? 0 : matrix[i][j].get(pathIndex);
					}
				}
			    // Get Pseudo Inverse 
			    RealMatrix pseudoInverse = new SingularValueDecomposition(MatrixUtils.createRealMatrix(matrixOnPath)).getSolver().getInverse();
			    for(int j=0;j<pseudoInverse.getColumnDimension();j++){
				    double[] columnValues = pseudoInverse.getColumn(j);
				    for(int i=0;i<pseudoInverse.getRowDimension();i++){
					    inv[i][j][pathIndex]= columnValues[i];
				    }
			    }
			}
			// Wrap to RandomVariableInterface[][]
			RandomVariableInterface[][] pseudoInverse = new RandomVariableInterface[matrix[0].length][matrix.length];
			for(int i=0;i<pseudoInverse.length; i++){
				for(int j=0;j<pseudoInverse[0].length; j++){
					pseudoInverse[i][j] = new RandomVariable(0.0 /*should be evaluationTime*/,inv[i][j]);
				}
			}
			return pseudoInverse;
	    }
	   

		public static RandomVariableInterface[][] multiply(RandomVariableInterface[][] A,RandomVariableInterface[][] B){
			RandomVariableInterface[][] AB = new RandomVariableInterface[A.length][B.length];
			RandomVariableInterface ABproduct;
			for(int i=0;i<A.length;i++){
				for(int j=0; j<B.length; j++){
					AB[i][j] = new RandomVariable(0.0);
					for(int k=0;k<B.length;k++) {
						if(A[i][k]==null || B[k][j]==null) {ABproduct = new RandomVariable(0.0);}
						else {ABproduct = A[i][k].mult(B[k][j]);}
						AB[i][j]=AB[i][j].add(ABproduct);
					}
				}
			}
			return AB;
		}
		
		public static RandomVariableInterface[] multiply(RandomVariableInterface[] A,RandomVariableInterface[][] B){
			RandomVariableInterface[] AB = new RandomVariableInterface[B[0].length];
			RandomVariableInterface ABproduct;
			for(int i=0;i<B[0].length;i++){
					AB[i] = new RandomVariable(0.0);
					for(int k=0;k<A.length;k++) {
						if(A[k]==null || B[k][i]==null) {ABproduct = new RandomVariable(0.0);}
						else {ABproduct = A[k].mult(B[k][i]);}
						AB[i]=AB[i].add(ABproduct);
					}
			}
			return AB;
		}
}
