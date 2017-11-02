package initialmargin.isdasimm;

import initialmargin.isdasimm.changedfinmath.products.AbstractLIBORMonteCarloProduct;

/** This class wraps an <code> AbstractLIBORMonteCarloProduct </code> into a product classified according 
 *  to the SIMM methodology requirement. All products to be considered for MVA calculation with the SIMM model
 *  must be classified by generating an object of this class.
 * 
 * @author Mario Viehmann
 *
 */
public class SIMMClassifiedProduct {
	
	    private AbstractLIBORMonteCarloProduct product; 
	    private boolean isCancelable;
	    
	    // Product classification within ISDA SIM
	    private String   productClass;      // RatesFX, Credit, 
	    private String[] riskClass;         // InterestRate, CreditQ, CreditNonQ, Equity, Commodity
	    private String[] curveIndexNames;   // e.g. OIS & Libor6m
	    private String   currency;
	    private boolean  hasOptionality;    // determines the relevance of vega and curvature risk (e.g. Swap has no curvature risk)
	    private String   bucketKey;         // can be null (e.g. in risk class InterestRate it is null because the bucket is given by the currency
	    
	    /**Wraps an <code> AbstractLIBORMonteCarloProduct </code> into a product classified according to the SIMM methodology requirement.
	     * 
	     * @param product The product to be classified
	     * @param productClass The SIMM product class of this product (RatesFx etc.)
	     * @param riskClass The SIMM risk class of this product (InterestRate etc.)
	     * @param curveIndexNames The name of the relevant curves for this product (OIS, Libor6m etc.)
	     * @param currency The currency of this product
	     * @param bucketKey The SIMM bucket key of this product (null for risk class InterestRate)
	     * @param hasOptionality True if this product is not linear
	     * @param isCancelable True if the product may be canceled before its maximum maturity
	     */
	    public SIMMClassifiedProduct(AbstractLIBORMonteCarloProduct product,
			                 String   productClass,
			                 String[] riskClass,     // One product may contribute to several risk Classes
			                 String[] curveIndexNames,
			                 String   currency,
			                 String   bucketKey,
			                 boolean  hasOptionality,
			                 boolean  isCancelable){
	    	
		   this.product=product;
		   this.productClass = productClass; 
		   this.riskClass = riskClass;
		   this.curveIndexNames = curveIndexNames;
		   this.currency=currency;
		   this.hasOptionality = hasOptionality;
		   this.bucketKey = bucketKey;
		   this.isCancelable = isCancelable;
	    }
	    
	    /*
	     * Getters
	     */
	    
	    public AbstractLIBORMonteCarloProduct getProduct(){
	    	return this.product;
	    }
	    public String getProductClass(){
	    	return this.productClass;
	    }
	    public String[] getRiskClasses(){
	    	return this.riskClass;
	    }
	    public String[] getCurveIndexNames(){
	    	return this.curveIndexNames;
	    }
	    public String getCurrency(){
	    	return this.currency;
	    }
	    public boolean getHasOptionality(){
	    	return this.hasOptionality;
	    }
	    public String getBucketKey(){
	    	return this.bucketKey;
	    }
	    public boolean getIsCancelable(){
	    	return this.isCancelable;
	    }
}
