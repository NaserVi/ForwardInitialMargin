package xva.tradespecifications;

import net.finmath.modelling.ProductInterface;
import xva.initialmargin.SIMMParameter;


import java.util.*;
import java.util.stream.Collectors;

public class SIMMTradeSpecification  {

    ProductInterface      underlyingValuationProduct;

    public static class SensitivityKey {
        private String maturityBucketKey;
        private String riskFactorKey;
        private String bucketKey;
        private SIMMParameter.RiskClass riskClass;
        private SIMMParameter.RiskType riskType;
        private SIMMParameter.ProductClass productClass;

        public SensitivityKey(String maturityBucket, String riskFactorID, String bucketID, String riskClass, String riskType, String productClass) {
            this.maturityBucketKey = maturityBucket;
            this.riskFactorKey = riskFactorID;
            this.bucketKey = bucketID;
            this.riskClass = SIMMParameter.RiskClass.valueOf(riskClass);
            this.riskType = SIMMParameter.RiskType.valueOf(riskType);
            this.productClass = SIMMParameter.ProductClass.valueOf(productClass);
        }

        public double getMaturityBucket() {
            return SensitivityKey.getMaturityBucket(this.maturityBucketKey);
        }

        static public double getMaturityBucket(String key) {
            if (key.contains("y"))
                return Double.parseDouble(key.replace("y", ""));
            if (key.contains("m"))
                return Double.parseDouble(key.replace("m", "")) / 12.0;
            return key.contains("w") ? Double.parseDouble(key.replace("w", "")) / 52.0 : 0.0;
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

        public SIMMParameter.RiskClass getRiskClass() {
            return riskClass;
        }

        public SIMMParameter.RiskType getRiskType() {
            return riskType;
        }

        public SIMMParameter.ProductClass getProductClass() {
            return productClass;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            SensitivityKey key = (SensitivityKey) o;

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




    Set<SensitivityKey> sensitivityKeySet;

    public SIMMTradeSpecification(double notional, double maturity, String IRCurveKey){

    }

    public double   getMaxTimeToMaturity(){
        return 0.0;
    }

    public double    getNotional(){
        return 0.0;
    }

    public SIMMParameter.ProductClass getProductClass(){
        return sensitivityKeySet.stream().map(key->key.getProductClass()).distinct().findAny().get();
    }

    public Set<SIMMParameter.RiskClass> getRiskClasses(){
        return sensitivityKeySet.stream().map(key->key.getRiskClass()).collect(Collectors.toSet());
    }


    public Set<String>  getRiskfactors(){
        return this.sensitivityKeySet.stream().map(key->key.getRiskFactorKey()).collect(Collectors.toSet());
    }

    public String    getTradeID(){
        return "";
    }

    public Set<SensitivityKey> getSensitivityKeySet(double evaluationTime) {
        return sensitivityKeySet;
    }



}
