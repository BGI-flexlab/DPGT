package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2;

public class VCFQualityControlArgumentCollection {
	 public enum Mode {
	        SNP,
	        INDEL,
	        BOTH
	    }

	    /**
	     * Use either SNP for recalibrating only SNPs (emitting indels untouched in the output VCF) or INDEL for indels (emitting SNPs untouched in the output VCF). There is also a BOTH option for recalibrating both SNPs and indels simultaneously, but this is meant for testing purposes only and should not be used in actual analyses.
	     */
	    public VCFQualityControlArgumentCollection.Mode MODE = VCFQualityControlArgumentCollection.Mode.SNP;

	    /**
	     * Generate a VQSR model using per-allele data instead of the default per-site data, assuming that the input VCF contains allele-specific annotations.
	     * Annotations should be specified using their full names with AS_ prefix. Non-allele-specific (scalar) annotations will be applied to all alleles.
	     */
	    public boolean useASannotations = false;

	    /**
	     * This parameter determines the maximum number of Gaussians that should be used when building a positive model
	     * using the variational Bayes algorithm.
	     */
	    public int MAX_GAUSSIANS = 8;

	    /**
	     * This parameter determines the maximum number of Gaussians that should be used when building a negative model
	     * using the variational Bayes algorithm. The actual maximum used is the smaller value between the mG and mNG
	     * arguments, meaning that if -mG is smaller than -mNG, -mG will be used for both. Note that this number should
	     * be small (e.g. 4) to achieve the best results.
	     */
	    public int MAX_GAUSSIANS_FOR_NEGATIVE_MODEL = 2;

	    /**
	     * This parameter determines the maximum number of VBEM iterations to be performed in the variational Bayes algorithm.
	     * The procedure will normally end when convergence is detected.
	     */
	    public int MAX_ITERATIONS = 150;

	    /**
	     * This parameter determines the number of k-means iterations to perform in order to initialize the means of
	     * the Gaussians in the Gaussian mixture model.
	     */
	    public int NUM_KMEANS_ITERATIONS = 100;

	    /**
	     * If a variant has annotations more than -std standard deviations away from mean, it won't be used for building
	     * the Gaussian mixture model.
	     */
	    public double STD_THRESHOLD = 10.0;

	    public double SHRINKAGE = 1.0;

	    public double DIRICHLET_PARAMETER = 0.001;

	    public double PRIOR_COUNTS = 20.0;

	    /**
	     * The number of variants to use in building the Gaussian mixture model. Training sets larger than this will be randomly downsampled.
	     */
	    protected int MAX_NUM_TRAINING_DATA = 2500000;

	    /**
	     * This parameter determines the minimum number of variants that will be selected from the list of worst scoring
	     * variants to use for building the Gaussian mixture model of bad variants.
	     */
	    public int MIN_NUM_BAD_VARIANTS = 1000;

	    /**
	     * Variants scoring lower than this threshold will be used to build the Gaussian model of bad variants.
	     */
	    public double BAD_LOD_CUTOFF = -5.0;

	    /**
	     * MQ is capped at a "max" value (60 for bwa-mem) when the alignment is considered perfect. Typically, a huge
	     * proportion of the reads in a dataset are perfectly mapped, which yields a distribution of MQ values with a
	     * blob below the max value and a huge peak at the max value. This does not conform to the expectations of the
	     * Gaussian mixture model of VQSR and has been observed to yield a ROC curve with a jump.
	     *
	     * This argument aims to mitigate this problem. Using MQCap = X has 2 effects:  (1) MQs are transformed by a scaled
	     * logit on [0,X] (+ epsilon to avoid division by zero) to make the blob more Gaussian-like and (2) the transformed
	     * MQ=X are jittered to break the peak into a narrow Gaussian.
	     *
	     * Beware that IndelRealigner, if used, adds 10 to MQ for successfully realigned indels. We recommend to either use
	     * --read-filter ReassignOriginalMQAfterIndelRealignment with HaplotypeCaller or use a MQCap=max+10 to take that
	     * into account.
	     *
	     * If this option is not used, or if MQCap is set to 0, MQ will not be transformed.
	     */
	    public int MQ_CAP = 0;

	    /**
	     * The following 2 arguments are hidden because they are only for testing different jitter amounts with and without logit transform.
	     * Once this will have been tested, and the correct jitter amount chosen (perhaps as a function of the logit range [0,max]) they can be removed.
	     */
	    public boolean NO_MQ_LOGIT = false;
	    
	    public double MQ_JITTER = 0.05;
}
