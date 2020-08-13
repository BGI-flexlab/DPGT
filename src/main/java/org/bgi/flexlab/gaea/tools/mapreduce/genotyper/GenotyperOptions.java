/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.genotyper;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.SNPGenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.AFCalcFactory;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;
import org.bgi.flexlab.gaea.util.pairhmm.PairHMM;
import org.seqdoop.hadoop_bam.SAMFormat;
import org.seqdoop.hadoop_bam.VCFFormat;

import java.util.ArrayList;
import java.util.List;

import static org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY;

/**
 * Created by zhangyong on 2016/12/20.
 */
public class GenotyperOptions extends GaeaOptions implements HadoopOptions {
    private final static String SOFTWARE_NAME = "Genotyper";
    private final static String SOFTWARE_VERSION = "1.0";

    /**
     * input alignment data
     */
    private String input = null;

    private boolean samFormat = false;

    /**
     * output directory
     */
    private String output = null;

    private boolean bcfFormat = false;

    /**
     * Gaea indexed reference
     */
    private String reference = null;

    /**
     * bed region file
     */
    private String bedRegionFile = null;

    /**
     * Which annotations to add to the output VCF file
     */
    private List<String> annotations = new ArrayList<String>();

    /**
     * annotation groups to add to the output VCF file
     */
    private List<String> annotationGroups = new ArrayList<String>();

    /**
     * models for genotype likelihood calculator
     */
    private GenotypeLikelihoodCalculator.Model gtlcalculators = GenotypeLikelihoodCalculator.Model.SNP;

    /**
     * min depth for genotype likelihood calculation
     */
    private int minDepth = 4;

    /**
     * is cap base quality at mapping quality
     */
    private boolean noCapBaseQualsAtMappingQual = false;

    /**
     * minimum base quality
     */
    private byte minBaseQuality = 17;

    /**
     * minimum mapping quality
     */
    private short minMappingQuality = 17;

    /**
     * output mode
     */
    private VariantCallingEngine.OUTPUT_MODE outputMode = EMIT_VARIANTS_ONLY;

    /**
     * Maximum fraction of reads with deletions spanning this locus for it to be callable
     */
    private double maxDeletionFraction = 0.05;

    /**
     * A candidate indel is genotyped (and potentially called) if there are this number of reads with a consensus indel at a site.
     * Decreasing this value will increase sensitivity but at the cost of larger calling time and a larger number of false positives.
     */
    private int minIndelCountForGenotyping = 5;

    /**
     * Complementary argument to minIndelCnt.  Only samples with at least this fraction of indel-containing reads will contribute
     * to counting and overcoming the threshold minIndelCnt.  This parameter ensures that in deep data you don't end
     * up summing lots of super rare errors up to overcome the 5 read default threshold.  Should work equally well for
     * low-coverage and high-coverage samples, as low coverage samples with any indel containing reads should easily over
     * come this threshold.
     */
    private double minIndelFractionPerSample = 0.25;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    private double contaminationFraction = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.05;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     */
    private int maxAlternateAlleles = 6;

    /**
     * Maximum number of genotypes to consider at any site
     */
    private int maxGenotypeCount = 1024;

    /**
     * Maximum number of PL values to output
     */
    private int maxNumPLValues = 100;

    /**
     * Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    private int samplePloidy = GaeaVariantContextUtils.DEFAULT_PLOIDY;

    /**
     * allele frequency calculation model
     */
    private AFCalcFactory.Calculation AFmodel = AFCalcFactory.Calculation.getDefaultModel();

    /**
     *
     */
    private PairHMM.HMM_IMPLEMENTATION pairHmmImplementation = PairHMM.HMM_IMPLEMENTATION.ORIGINAL;

    /**
     * The minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. Only genotypes with
     * confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 for high-pass calling (this
     * is the default).
     */
    private double standardConfidenceForCalling = 10.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    private double standardConfidenceForEmitting = 30.0;

    /**
     * The expected heterozygosity value used to compute prior likelihoods for any locus. The default priors are:
     * het = 1e-3, P(hom-ref genotype) = 1 - 3 * het / 2, P(het genotype) = het, P(hom-var genotype) = het / 2
     */
    private Double heterozygosity = VariantCallingEngine.HUMAN_SNP_HETEROZYGOSITY;

    /**
     * Standard deviation of eterozygosity for SNP and indel calling.
     */
    private double heterozygosityStandardDeviation = 0.01;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    private double indelHeterozygosity = 1.0/8000;

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    private boolean annotateNumberOfAllelesDiscovered = false;

    /**
     * The PCR error rate to be used for computing fragment-based likelihoods
     */
    private double pcr_error = SNPGenotypeLikelihoodCalculator.DEFAULT_PCR_ERROR_RATE;

    /**
     * ndel gap continuation penalty, as Phred-scaled probability. I.e., 30 => 10^-30/10
     */
    private byte indelGapContinuationPenalty = 10;

    /**
     * Indel gap open penalty, as Phred-scaled probability. I.e., 30 => 10^-30/10
     */
    private byte indelGapOpenPenalty = 45;

    /**
     * multi sample mode
     */
    private boolean singleSampleMode;

    /**
     * reducer number
     */
    private int reducerNumber = 30;

    /**
     * window size
     */
    private int windowSize = 10000;

    public GenotyperOptions() {
        addOption("i", "input", true, "Input file containing sequence data (BAM or CRAM)");
        addOption("I", "is_sam_input", false, "the input is in SAM format.");
        addOption("o", "output", true, "directory to which variants should be written in HDFS written format.");
        addOption("O", "outputMode", true, "output mode :EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES");
        addOption("B","is_bcf_output", false, "output variant in BCF format." );
        addOption("b", "bed_file", true, "only variant in this region will be called.");
        addOption("r", "reference", true, "Gaea Indexed reference list for memory sharing.");
        addOption("A", "annotation", true, "One or more specific annotations to apply to variant calls, the tag is separated by \',\'");
        addOption("contamination", "contamination_fraction_to_filter", true, "Fraction of contamination to aggressively remove.");
        addOption("glm", "genotype_likelihoods_model", true, "Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together.");
        addOption("G", "group", true, "One or more classes/groups of annotations to apply to variant calls. The single value 'none' removes the default group, the group tag is separated by \',\'");
        addOption("hets", "heterozygosity", true, "Heterozygosity value used to compute prior likelihoods for any locus.");
        addOption("heterozygosityStandardDeviation", "heterozygosity_stdev", true, "\tStandard deviation of eterozygosity for SNP and indel calling.");
        addOption("indelHeterozygosity", "indel_heterozygosity", true, "Heterozygosity for indel calling.");
        addOption("deletions", "max_deletion_fraction", true, "Maximum fraction of reads with deletions spanning this locus for it to be callable.");
        addOption("mbq", "min_base_quality_score", true, "Minimum base quality required to consider a base for calling.");
        addOption("mmq", "min_mapping_quality_score", true, "Minimum mapping quality required to consider a read for calling.");
        addOption("minIndelCnt", "min_indel_count_for_genotyping", true, "Minimum number of consensus indels required to trigger genotyping run.");
        addOption("minIndelFrac", "min_indel_fraction_per_sample", true, "Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles.");
        addOption("pairHMM", "pair_hmm_implementation", true, "The PairHMM implementation to use for -glm INDEL genotype likelihood calculations");
        addOption("pcrError", "pcr_error_rate", true, "The PCR error rate to be used for computing fragment-based likelihoods.");
        addOption("ploidy", "sample_ploidy", true, "Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).");
        addOption("standCallConf", "standardConfidenceForCalling", true, "standard confidence for calling");
        addOption("standEmitConf", "standardConfidenceForEmitting", true, "standard confidence for emitting");
        addOption("AFmodel", "Allele_frequency_model", true, "allele frequency calculation model");
        addOption("indelGCP", "indelGapContinuationPenalty", true, "Indel gap continuation penalty, as Phred-scaled probability. I.e., 30 => 10^-30/10");
        addOption("indelGOP", "indelGapOpenPenalty", true, "Indel gap open penalty, as Phred-scaled probability. I.e., 30 => 10^-30/10.");
        addOption("maxAltAlleles", "max_alternate_alleles", true, "max alternate alleles");
        addOption("maxGT", "max_genotype_count", true, "Maximum number of genotypes to consider at any site.");
        addOption("maxNumPLValues", "max_num_PL_values", true, "Maximum number of PL values to output");
        addOption("S", "single_sample_mode", false, "will call genotype and variant for each sample separately");
        addOption("D", "minDepth", true, "minimum depth for variant calling.");
        addOption("C", "noCapBaseQualsAtMappingQual", false, "do not cap base quality at mapping quality");
        addOption("numAlleleDis", "annotateNumberOfAllelesDiscovered", false, "annotate Number Of Alleles Discovered");
        addOption("R", "reducer", true, "reducer numbers");
        addOption("W", "window_size", true, "window size that sharding the data.");
        addOption("h", "help", false, "print help information.");

        FormatHelpInfo(SOFTWARE_NAME, SOFTWARE_VERSION);
    }

    @Override
    public void parse(String[] args) {
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            printHelpInfotmation(SOFTWARE_NAME);
            throw new RuntimeException(e);
        }

        if(args.length == 0 || getOptionBooleanValue("h", false)) {
            printHelpInfotmation(SOFTWARE_NAME);
            System.exit(1);
        }
        input = getOptionValue("i", null);
        samFormat = getOptionBooleanValue("I", false);
        output = getOptionValue("o", null);
        bcfFormat = getOptionBooleanValue("B", false);
        bedRegionFile = getOptionValue("b", null);
        reference = getOptionValue("r", null);
        String annotationTags = getOptionValue("A", null);
        if(annotationTags != null) {
            for(String tag : annotationTags.split(","))
            annotations.add(tag);
        }
        annotationGroups.add("Standard");
        String annotationGroupTags = getOptionValue("G", null);
        if(annotationGroupTags != null) {
            if(annotationGroupTags.equals("none")) {
                annotationGroups.clear();
                annotationGroups.add("none");
            } else {
                for (String tag : annotationGroupTags.split(","))
                    annotationGroups.add(tag);
            }
        }
        contaminationFraction = getOptionDoubleValue("C1", DEFAULT_CONTAMINATION_FRACTION);
        heterozygosity = getOptionDoubleValue("hets", VariantCallingEngine.HUMAN_SNP_HETEROZYGOSITY);
        heterozygosityStandardDeviation = getOptionDoubleValue("heterozygosityStandardDeviation", 0.01);
        indelHeterozygosity = getOptionDoubleValue("indelHeterozygosity", 1.0/8000);
        maxDeletionFraction = getOptionDoubleValue("deletions", 0.05);
        noCapBaseQualsAtMappingQual = getOptionBooleanValue("C", false);
        minBaseQuality = getOptionByteValue("mbq", (byte)17);
        minMappingQuality = getOptionShortValue("mmq", (short)17);
        minIndelCountForGenotyping = getOptionIntValue("minIndelCnt", 5);
        minIndelFractionPerSample = getOptionDoubleValue("minIndelFrac", 0.25);
        samplePloidy = getOptionIntValue("ploidy", GaeaVariantContextUtils.DEFAULT_PLOIDY);
        standardConfidenceForCalling = getOptionDoubleValue("standCallConf", 10.0);
        standardConfidenceForEmitting = getOptionDoubleValue("standEmitConf", 30);
        indelGapContinuationPenalty = getOptionByteValue("indelGCP", (byte) 10);
        indelGapOpenPenalty = getOptionByteValue("indelGOP", (byte) 45);
        maxAlternateAlleles = getOptionIntValue("maxAltAlleles", 6);
        maxGenotypeCount = getOptionIntValue("maxGT", 1024);
        maxNumPLValues = getOptionIntValue("maxNumPLValues", 100);
        singleSampleMode = getOptionBooleanValue("S", false);
        reducerNumber = getOptionIntValue("R", 30);
        windowSize = getOptionIntValue("W", 100000);
        pcr_error = getOptionDoubleValue("pcrError", SNPGenotypeLikelihoodCalculator.DEFAULT_PCR_ERROR_RATE);

        try {
            gtlcalculators = GenotypeLikelihoodCalculator.Model.valueOf(getOptionValue("glm", "SNP"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("glm", e.getMessage());
        }
        try {
            outputMode = VariantCallingEngine.OUTPUT_MODE.valueOf(getOptionValue("O", "EMIT_VARIANTS_ONLY"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("O", e.getMessage());
        }
        try {
            AFmodel = AFCalcFactory.Calculation.valueOf(getOptionValue("AFmodel", "EXACT_INDEPENDENT"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("AFmodel", e.getMessage());
        }
        try {
            pairHmmImplementation = PairHMM.HMM_IMPLEMENTATION.valueOf(getOptionValue("pairHMM", "ORIGINAL"));
        } catch (Exception e) {
            throw new UserException.BadArgumentValueException("pairHMM", e.getMessage());
        }

        check();

        minDepth = getOptionIntValue("D", 4);
        annotateNumberOfAllelesDiscovered = getOptionBooleanValue("numAlleleDis", false);
        noCapBaseQualsAtMappingQual = getOptionBooleanValue("C", false);

    }

    private void check() {
        if(input == null)
            throw new UserException.BadArgumentValueException("i", "input directory or file is not assigned.");

        if(output == null)
            throw new UserException.BadArgumentValueException("o", "output directory is not assigned.");

        if(reference == null)
            throw new UserException.BadArgumentValueException("r", "reference can not be null.");

        if(reducerNumber <= 0 ) {
            throw new UserException.BadArgumentValueException("R", "reducer number can not be less than 1.");
        }
    }

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        conf.setStrings("args", args);
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    public GenotypeLikelihoodCalculator.Model getGtlcalculators() {
        return gtlcalculators;
    }

    public void setGtlcalculators(GenotypeLikelihoodCalculator.Model gtlcalculators) {
        this.gtlcalculators = gtlcalculators;
    }

    public SAMFormat getInputFormat(){
        if(samFormat)
            return SAMFormat.SAM;
        return SAMFormat.BAM;
    }

    public VCFFormat getOuptputFormat() {
        if(bcfFormat)
            return VCFFormat.BCF;
        return VCFFormat.VCF;
    }

    public int getMinDepth() {
        return minDepth;
    }

    public void setMinDepth(int minDepth) {
        this.minDepth = minDepth;
    }

    public boolean isCapBaseQualsAtMappingQual() {
        return !noCapBaseQualsAtMappingQual;
    }

    public byte getMinBaseQuality() {
        return minBaseQuality;
    }

    public short getMinMappingQuality() {
        return minMappingQuality;
    }

    public VariantCallingEngine.OUTPUT_MODE getOutputMode() {
        return outputMode;
    }

    public int getMinIndelCountForGenotyping() {
        return minIndelCountForGenotyping;
    }

    public double getMinIndelFractionPerSample() {
        return minIndelFractionPerSample;
    }

    public double getContaminationFraction() {
        return contaminationFraction;
    }

    public int getMaxAlternateAlleles() {
        return maxAlternateAlleles;
    }

    public int getSamplePloidy() {
        return samplePloidy;
    }

    public AFCalcFactory.Calculation getAFmodel() {
        return AFmodel;
    }

    public void setAFmodel(AFCalcFactory.Calculation AFmodel) {
        this.AFmodel = AFmodel;
    }

    public double getStandardConfidenceForCalling() {
        return standardConfidenceForCalling;
    }

    public double getStandardConfidenceForEmitting() {
        return standardConfidenceForEmitting;
    }

    public double getHeterozygosity() {
        return heterozygosity;
    }

    public double getIndelHeterozygosity() {
        return indelHeterozygosity;
    }

    public boolean isAnnotateNumberOfAllelesDiscovered() {
        return annotateNumberOfAllelesDiscovered;
    }

    public boolean isSingleSampleMode() {
        return singleSampleMode;
    }

    public String getInput() {
        return input;
    }

    public String getOutput() {
        return output;
    }

    public String getBAMHeaderOutput(){
        return output+"bamHeader";
    }

    public String getVCFHeaderOutput() {
        return output + "vcfHeader";
    }

    public String getReference() {
        return  reference;
    }

    public int getReducerNumber() {
        return reducerNumber;
    }

    public int getWindowSize() {
        return windowSize;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    public List<String> getAnnotationGroups() {
        return annotationGroups;
    }

    public String getBedRegionFile() {
        return bedRegionFile;
    }

    public double getPcr_error() {
        return pcr_error;
    }

    public byte getIndelGapOpenPenalty() {
        return indelGapOpenPenalty;
    }

    public byte getIndelGapContinuationPenalty() {
        return indelGapContinuationPenalty;
    }

    public PairHMM.HMM_IMPLEMENTATION getPairHmmImplementation() {
        return pairHmmImplementation;
    }
}

