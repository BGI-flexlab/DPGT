package org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections4.map.DefaultedMap;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.HomoSapiensConstants;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.GenotypingOutputMode;
import org.bgi.flexlab.gaea.tools.jointcalling.UnifiedGenotypingEngine.OutputMode;
import org.bgi.flexlab.gaea.tools.jointcalling.afcalculator.AFCalculatorImplementation;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.VariantContext;

public class StandardCallerArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    
    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    public AFCalculatorImplementation requestedAlleleFrequencyCalculationModel;

    public OutputMode outputMode = OutputMode.EMIT_VARIANTS_ONLY;

    /**
     * Advanced, experimental argument: if SNP likelihood model is specified, and if EMIT_ALL_SITES output mode is set, when we set this argument then we will also emit PLs at all sites.
     * This will give a measure of reference confidence and a measure of which alt alleles are more plausible (if any).
     * WARNINGS:
     * - This feature will inflate VCF file size considerably.
     * - All SNP ALT alleles will be emitted with corresponding 10 PL values.
     * - An error will be emitted if EMIT_ALL_SITES is not set, or if anything other than diploid SNP model is used
     */
    public boolean annotateAllSitesWithPLs = false;

    public GenotypingOutputMode genotypingOutputMode = GenotypingOutputMode.DISCOVERY;

    /**
     * When the caller is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    public List<VariantContext> alleles;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.0;

    /**
     *  This argument specifies a file with two columns "sample" and "contamination" specifying the contamination level for those samples.
     *  Samples that do not appear in this file will be processed with CONTAMINATION_FRACTION.
     **/
    public File CONTAMINATION_FRACTION_FILE = null;
    
    private DefaultedMap<String,Double> sampleContamination;
    
    /**
     * Use the new allele frequency / QUAL score model
     */
    public boolean USE_NEW_AF_CALCULATOR = false;

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * The expected heterozygosity value used to compute prior probability that a locus is non-reference.
     *
     * The default priors are for provided for humans:
     *
     * het = 1e-3
     *
     * which means that the probability of N samples being hom-ref at a site is:
     *
     * 1 - sum_i_2N (het / i)
     */
    public Double snpHeterozygosity = HomoSapiensConstants.SNP_HETEROZYGOSITY;

    /**
     * This argument informs the prior probability of having an indel at a site.
     */
    public double indelHeterozygosity = HomoSapiensConstants.INDEL_HETEROZYGOSITY;

    /**
     * The standard deviation of the distribution of alt allele fractions.  The above heterozygosity parameters give the
     * *mean* of this distribution; this parameter gives its spread.
     */
    public double heterozygosityStandardDeviation = 0.01;

    /**
     * The minimum phred-scaled confidence threshold at which variants should be called. Only variant sites with QUAL equal
     * or greater than this threshold will be called.
     */
    public double STANDARD_CONFIDENCE_FOR_CALLING = 10.0;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     */
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * If there are more than this number of genotypes at a locus presented to the genotyper, then only this many genotypes will be used.
     * The possible genotypes are simply different ways of partitioning alleles given a specific ploidy asumption.
     * Therefore, we remove genotypes from consideration by removing alternate alleles that are the least well supported.
     * The estimate of allele support is based on the ranking of the candidate haplotypes coming out of the graph building step.
     * Note that the reference allele is always kept.
     *
     * Note that genotyping sites with large genotype counts is both CPU and memory intensive.
     * Unless there is a good reason to change the default value, we highly recommend that you not play around with this parameter.
     *
     * The maximum number of alternative alleles used in the genotyping step will be the lesser of the two:
     * 1. the largest number of alt alleles, given ploidy, that yields a genotype count no higher than {@link #MAX_GENOTYPE_COUNT}
     * 2. the value of {@link #MAX_ALTERNATE_ALLELES}
     */
    public int MAX_GENOTYPE_COUNT = 1024;

    /**
     * By default, the prior specified with the argument --heterozygosity/-hets is used for variant discovery at a particular locus, using an infinite sites model,
     * see e.g. Waterson (1975) or Tajima (1996).
     * This model asserts that the probability of having a population of k variant sites in N chromosomes is proportional to theta/k, for 1=1:N
     *
     * There are instances where using this prior might not be desireable, e.g. for population studies where prior might not be appropriate,
     * as for example when the ancestral status of the reference allele is not known.
     * By using this argument, user can manually specify priors to be used for calling as a vector for doubles, with the following restriciotns:
     * a) User must specify 2N values, where N is the number of samples.
     * b) Only diploid calls supported.
     * c) Probability values are specified in double format, in linear space.
     * d) No negative values allowed.
     * e) Values will be added and Pr(AC=0) will be 1-sum, so that they sum up to one.
     * f) If user-defined values add to more than one, an error will be produced.
     *
     * If user wants completely flat priors, then user should specify the same value (=1/(2*N+1)) 2*N times,e.g.
     *   -inputPrior 0.33 -inputPrior 0.33
     * for the single-sample diploid case.
     */
    public List<Double> inputPrior = Collections.emptyList();

    /**
     *   Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    public int samplePloidy = HomoSapiensConstants.DEFAULT_PLOIDY;
    
    public boolean dontTrimActiveRegions = false;

    /**
     * the maximum extent into the full active region extension that we're willing to go in genotyping our events
     */
    public int discoverExtension = 25;

    public int ggaExtension = 300;

    /**
     * Include at least this many bases around an event for calling it
     */
    public int indelPadding = 150;

    public int snpPadding = 20;

    /**
     * Copies the values from other into this StandardCallerArgumentCollection
     *
     * @param other StandardCallerArgumentCollection from which to copy values
     */
    public void copyStandardCallerArgsFrom( final StandardCallerArgumentCollection other ) {
        Utils.nonNull(other);
        this.genotypingOutputMode = other.genotypingOutputMode;
        this.alleles = other.alleles; // FeatureInputs are immutable outside of the engine, so this shallow copy is safe
        this.CONTAMINATION_FRACTION = other.CONTAMINATION_FRACTION;
        this.CONTAMINATION_FRACTION_FILE = other.CONTAMINATION_FRACTION_FILE != null ? new File(other.CONTAMINATION_FRACTION_FILE.getAbsolutePath()) : null;
        if ( other.sampleContamination != null ) {
            setSampleContamination(other.sampleContamination);
        }
        this.requestedAlleleFrequencyCalculationModel = other.requestedAlleleFrequencyCalculationModel;
        this.outputMode = other.outputMode;
        this.annotateAllSitesWithPLs = other.annotateAllSitesWithPLs;
        
        this.USE_NEW_AF_CALCULATOR = other.USE_NEW_AF_CALCULATOR;
        this.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = other.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED;
        this.snpHeterozygosity = other.snpHeterozygosity;
        this.indelHeterozygosity = other.indelHeterozygosity;
        this.STANDARD_CONFIDENCE_FOR_CALLING = other.STANDARD_CONFIDENCE_FOR_CALLING;
        this.MAX_ALTERNATE_ALLELES = other.MAX_ALTERNATE_ALLELES;
        this.inputPrior = new ArrayList<>(other.inputPrior);
        this.samplePloidy = other.samplePloidy;
    }

    

    /**
     * Returns true if there is some sample contamination present, false otherwise.
     * @return {@code true} iff there is some sample contamination
     */
    public boolean isSampleContaminationPresent() {
        return (!Double.isNaN(CONTAMINATION_FRACTION) && CONTAMINATION_FRACTION > 0.0) || (sampleContamination != null && !sampleContamination.isEmpty());
    }

    /**
     * Returns an unmodifiable view of the map of SampleId -> contamination.
     */
    public Map<String,Double> getSampleContamination() {
        return Collections.unmodifiableMap(sampleContamination);
    }

    /**
     * Returns the sample contamination or CONTAMINATION_FRACTION if no contamination level was specified for this sample.
     */
    public Double getSampleContamination(final String sampleId){
        Utils.nonNull(sampleId);
        if (sampleContamination == null){
            setSampleContamination(new DefaultedMap<>(CONTAMINATION_FRACTION));//default to empty map
        }
        return sampleContamination.get(sampleId);
    }

    public void setSampleContamination(final DefaultedMap<String, Double> sampleContamination) {
        this.sampleContamination = new DefaultedMap<>(CONTAMINATION_FRACTION);  //NOTE: a bit weird because it ignores the default from the argument and uses ours
        this.sampleContamination.putAll(sampleContamination);                   //make a copy to be safe
    }    
}

