package org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection;

import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pairhmm.PairHMM;


public class UnifiedArgumentCollection extends StandardCallerArgumentCollection {
	
	private static final long serialVersionUID = 1L;

    public UnifiedGenotypingEngine.Model GLmodel = UnifiedGenotypingEngine.Model.SNP;

    /**
     * The PCR error rate is independent of the sequencing error rate, which is necessary because we cannot necessarily
     * distinguish between PCR errors vs. sequencing errors.  The practical implication for this value is that it
     * effectively acts as a cap on the base qualities.
     */
    public Double PCR_error = 1e-4;

    /**
     * Note that calculating the SLOD increases the runtime by an appreciable amount.
     */
    public boolean COMPUTE_SLOD = false;

    /**
     * The PairHMM implementation to use for -glm INDEL genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime.
     */
    public PairHMM.Implementation pairHMM = PairHMM.Implementation.FASTEST_AVAILABLE;

    /**
     * The minimum confidence needed in a given base for it to be used in variant calling.  Note that the base quality of a base
     * is capped by the mapping quality so that bases on reads with low mapping quality may get filtered out depending on this value.
     * Note too that this argument is ignored in indel calling. In indel calling, low-quality ends of reads are clipped off (with fixed threshold of Q20).
     */
    public int MIN_BASE_QUALTY_SCORE = 17;

    /**
     * If the fraction of reads with deletions spanning a locus is greater than this value, the site will not be considered callable and will be skipped.
     * To disable the use of this parameter, set its value to >1.
     */
    public Double MAX_DELETION_FRACTION = 0.05;

    // indel-related arguments
    /**
     * A candidate indel is genotyped (and potentially called) if there are this number of reads with a consensus indel at a site.
     * Decreasing this value will increase sensitivity but at the cost of larger calling time and a larger number of false positives.
     */
    public int MIN_INDEL_COUNT_FOR_GENOTYPING = 5;

    /**
     * Complementary argument to minIndelCnt.  Only samples with at least this fraction of indel-containing reads will contribute
     * to counting and overcoming the threshold minIndelCnt.  This parameter ensures that in deep data you don't end
     * up summing lots of super rare errors up to overcome the 5 read default threshold.  Should work equally well for
     * low-coverage and high-coverage samples, as low coverage samples with any indel containing reads should easily over
     * come this threshold.
     */
    public double MIN_INDEL_FRACTION_PER_SAMPLE = 0.25;

    public byte INDEL_GAP_CONTINUATION_PENALTY = 10;

    public byte INDEL_GAP_OPEN_PENALTY = 45;

    public int INDEL_HAPLOTYPE_SIZE = 80;

    public boolean OUTPUT_DEBUG_INDEL_INFO = false;

    public boolean IGNORE_SNP_ALLELES = false;

    /*
        Generalized ploidy argument (debug only): squash all reads into a single pileup without considering sample info
     */
    public boolean TREAT_ALL_READS_AS_SINGLE_POOL = false;

    /*
       Generalized ploidy argument (debug only): When building site error models, ignore lane information and build only
       sample-level error model
     */
    public boolean IGNORE_LANE_INFO = false;

    /*
        Reference sample name: if included, a site-specific error model will be built in order to improve calling quality. This requires ideally
        that a bar-coded reference sample be included with the polyploid/pooled data in a sequencing experimental design.
        If argument is absent, no per-site error model is included and calling is done with a generalization of traditional statistical calling.
     */
    public String referenceSampleName;

    /**
     * The following argument are for debug-only tweaks when running generalized ploidy with a reference sample
     */
    public byte minQualityScore= 1;

    public byte maxQualityScore= 40;

    public byte phredScaledPrior = 20;

    public double minPower = 0.95;
}
