package initialmargin.isdasimm;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;

public class SIMMSchemeSensitivitySet {

    public static class Key {
        public String maturityBucketKey;
        public String riskFactorKey;
        public String bucketKey;
        public String riskClass;
        public String riskType;
        public String productClass;


        /* FOR IRRiskType */
        public Key(String maturityBucket, String riskFactorID, String bucketID, String riskClass, String riskType, String productClass) {
            this.maturityBucketKey = maturityBucket;
            this.riskFactorKey = riskFactorID;
            this.bucketKey = bucketID;
            this.riskClass = riskClass;
            this.riskType = riskType;
            this.productClass = productClass;
        }

        public double getMaturityBucket() {
            return Key.getMaturityBucket(this.maturityBucketKey);
        }

        static public double getMaturityBucket(String key) {
            if (key.contains("y"))
                return Double.parseDouble(key.replace("y", ""));
            if (key.contains("m"))
                return Double.parseDouble(key.replace("m", "")) / 12.0;
            if (key.contains("w"))
                return Double.parseDouble(key.replace("w", "")) / 52.0;
            else
                return 0.0;
        }

        public String getMaturityBucketKey() {
            return maturityBucketKey;
        }

        public String getRiskFactorKey() {
            return riskFactorKey;
        }

        public String getBucketKey() {
            return bucketKey;
        }

        public String getRiskClass() {
            return riskClass;
        }

        public String getRiskType() {
            return riskType;
        }

        public String getProductClass() {
            return productClass;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Key key = (Key) o;

            if (maturityBucketKey != null ? !maturityBucketKey.equals(key.maturityBucketKey) : key.maturityBucketKey != null)
                return false;
            if (riskFactorKey != null ? !riskFactorKey.equals(key.riskFactorKey) : key.riskFactorKey != null)
                return false;
            if (bucketKey != null ? !bucketKey.equals(key.bucketKey) : key.bucketKey != null) return false;
            if (riskClass != null ? !riskClass.equals(key.riskClass) : key.riskClass != null) return false;
            if (riskType != null ? !riskType.equals(key.riskType) : key.riskType != null) return false;
            return productClass != null ? productClass.equals(key.productClass) : key.productClass == null;

        }

        @Override
        public int hashCode() {
            int result = maturityBucketKey != null ? maturityBucketKey.hashCode() : 0;
            result = 31 * result + (riskFactorKey != null ? riskFactorKey.hashCode() : 0);
            result = 31 * result + (bucketKey != null ? bucketKey.hashCode() : 0);
            result = 31 * result + (riskClass != null ? riskClass.hashCode() : 0);
            result = 31 * result + (riskType != null ? riskType.hashCode() : 0);
            result = 31 * result + (productClass != null ? productClass.hashCode() : 0);
            return result;
        }
    }

    Map<Key, RandomVariableInterface> sensitivityMap;

    Map<String,String[]> maturityBucketMap;

    final static double simmVolatilityFactor = Math.sqrt(365./14.) / (new NormalDistribution(0,1).inverseCumulativeProbability(.99));

    int nPathDimension;

    public SIMMSchemeSensitivitySet() {
        this.sensitivityMap = new HashMap<>();
        maturityBucketMap = new HashMap<>();
        nPathDimension = 1;
    }

    public SIMMSchemeSensitivitySet(Map<Key, RandomVariableInterface> sensitivityMap) {
        this.sensitivityMap = sensitivityMap;
        maturityBucketMap = new HashMap<>();
        nPathDimension = 1;
    }

    public  SIMMSchemeSensitivitySet getCloned(){
        SIMMSchemeSensitivitySet clone = new SIMMSchemeSensitivitySet();
        clone.sensitivityMap.putAll(this.sensitivityMap);
        clone.maturityBucketMap.putAll(this.maturityBucketMap);
        clone.nPathDimension = this.nPathDimension;
        return clone;
    }


    public Set<Key> getKeySet(){
        return this.sensitivityMap.keySet();
    }

    public int getPathDimension()
    {
        return nPathDimension;
    }


    public void    setMaturityBuckets(String riskClass, String[] bucketKeys){
        this.maturityBucketMap.put(riskClass,bucketKeys);
    }



    public RandomVariableInterface  getSensitivity(SIMMSchemeSensitivitySet.Key key, double atTime)
    {
    	return new RandomVariable(atTime,nPathDimension,0.0);
       /* try{
            AbstractExposure exposure = this.sensitivityMap.get(key);
            RandomVariableInterface sensiValue = exposure.getExposureValue(atTime);
            if (key.riskType.equals("vega"))
            {
                if ( key.riskClass.equals("Equity") || key.riskClass.equals("Commodity") || key.riskClass.equals("FX"))
                {
                    sensiValue = sensiValue.mult(simmVolatilityFactor);
                }
            }
            if(sensiValue.size() == this.nPathDimension)
                return sensiValue;
            else
                return new RandomVariable(atTime,nPathDimension,sensiValue.get(0));

        }
        catch(Exception e){
            return new RandomVariable(atTime,nPathDimension,0.0);
        }*/
    }




}