package org.bgi.flexlab.dpgt.jointcalling;

import java.nio.file.Paths;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;
import java.util.Set;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.HashMap;
import java.util.TreeSet;
import java.io.File;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.OutputMode;
import org.broadinstitute.hellbender.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;


public class GVCFsSyncGenotyper {
    private static final Logger logger = LoggerFactory.getLogger(GVCFsSyncGenotyper.class);
    public static final String PHASED_HOM_VAR_STRING = "1|1";
    private static final String GVCF_BLOCK = "GVCFBlock";

    private FeatureContext emptyFeatureContext = new FeatureContext();
    private ReferenceDataSource reference;
    private MultiVariantSyncReader reader = null;
    private VariantAnnotatorEngine annotationEngine;
    private GenotypingEngine<?> genotypingEngine;
    private ReferenceConfidenceVariantContextMerger merger;
    private File outputFile;
    private VariantContextWriter vcfWriter = null;

    private DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();

    private final boolean includeNonVariants = false;
    private GenotypeCalculationArgumentCollection genotypeArgs = null;

    public GVCFsSyncGenotyper() {}

    public GVCFsSyncGenotyper(final String refpath, final List<String> vcfpaths, final String vcfHeader,
        final SimpleInterval interval, final String outpath, final String dbsnpPath, final GenotypeCalculationArgumentCollection genotypeArgs) {
        
        // init reference
        reference = ReferenceDataSource.of(Paths.get(refpath));
        reader = new MultiVariantSyncReader();
        reader.open(vcfpaths, vcfHeader);
        reader.query(interval);
        final VCFHeader inputVCFHeader = reader.getVCFHeader();

        // precalculate a large log10 value to avoid repeatly expand log10cache
        {
            MathUtils.log10(inputVCFHeader.getNGenotypeSamples()*genotypeArgs.samplePloidy*2);
        }
        
        if (dbsnpPath != null) {
            dbsnp.dbsnp = new FeatureInput<>(dbsnpPath, "dnsnp");
        } else {
            dbsnp.dbsnp = null;
        }
        
        List<Annotation> annotations = makeVariantAnnotations();
        annotationEngine = new VariantAnnotatorEngine(annotations, dbsnp.dbsnp, Collections.emptyList(), false);

        // Request INFO field annotations inheriting from RankSumTest and RMSAnnotation added to remove list
        for ( final InfoFieldAnnotation annotation :  annotationEngine.getInfoAnnotations() ) {
            if ( annotation instanceof RankSumTest ||
                    annotation instanceof AS_RMSMappingQuality ||
                    annotation instanceof RMSMappingQuality) {
                final List<String> keyNames = annotation.getKeyNames();
                if ( !keyNames.isEmpty() ) {
                    infoFieldAnnotationKeyNamesToRemove.add(keyNames.get(0));
                }
            }
        }

        final IndexedSampleList samples = new IndexedSampleList(inputVCFHeader.getSampleNamesInOrder());

        this.genotypeArgs = genotypeArgs;
        // We only want the engine to generate the AS_QUAL key if we are using AlleleSpecific annotations.
        genotypingEngine = new MinimalGenotypingEngine(createUAC(), samples, new GeneralPloidyFailOverAFCalculatorProvider(this.genotypeArgs), annotationEngine.isRequestedReducibleRawKey(GATKVCFConstants.AS_QUAL_KEY));

        merger = new ReferenceConfidenceVariantContextMerger(annotationEngine, inputVCFHeader);

        outputFile = new File(outpath);

        vcfWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);
        vcfWriter.setHeader(inputVCFHeader);
    }

    /**
     * make combined vcf header for genotyping by add annotation and genotyping header lines
     * @param vcfHeaderPath combined vcf header for input vcfs
     * @return combined vcf header for genotyping 
     */
    public VCFHeader makeCombinedHeaderForGenotyping(final String vcfHeaderPath, final GenotypeCalculationArgumentCollection genotypeArgs) {
        List<Annotation> annotations = makeVariantAnnotations();
        annotationEngine = new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false);

        this.genotypeArgs = genotypeArgs;
        VCFFileReader vcfHeaderReader = new VCFFileReader(Paths.get(vcfHeaderPath));
        final VCFHeader inputVCFHeader = vcfHeaderReader.getHeader();
        final IndexedSampleList samples = new IndexedSampleList(inputVCFHeader.getSampleNamesInOrder());
        // We only want the engine to generate the AS_QUAL key if we are using AlleleSpecific annotations.
        genotypingEngine = new MinimalGenotypingEngine(createUAC(), samples, new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs), annotationEngine.isRequestedReducibleRawKey(GATKVCFConstants.AS_QUAL_KEY));

        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader combinedVCFHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));

        vcfHeaderReader.close();

        return combinedVCFHeader;
    }

    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        // vcfWriter.writeHeader(vcfHeader);
    }


    public void run() {
        ArrayList<VariantContext> variantContexts;
        while(!(variantContexts = reader.read()).isEmpty()) {
            final SimpleInterval variantInterval = new SimpleInterval(variantContexts.get(0));
            if (containsTrueAltAllele(variantContexts)) {
                apply(variantInterval, variantContexts, new ReferenceContext(reference, variantInterval));
            }
        }
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    public void stop() {
        if (reader != null) {
            reader.close();
        }
    }

    public void apply(final Locatable loc, List<VariantContext> variants, ReferenceContext ref) {
        final VariantContext mergedVC = merger.merge(variants, loc, includeNonVariants ? ref.getBase() : null, !includeNonVariants, false);
        final VariantContext regenotypedVC = regenotypeVC(mergedVC, ref, emptyFeatureContext, includeNonVariants);
        if (regenotypedVC != null) {
            if (!GATKVariantContextUtils.isSpanningDeletionOnly(regenotypedVC)) {
                vcfWriter.add(regenotypedVC);
            }
        }
    }

    /**
     * Does the given list of VariantContexts contain any with an alternate allele other than <NON_REF>?
     *
     * @param VCs  list of VariantContexts
     * @return true if there are one or more variantContexts that contain a true alternate allele, false otherwise
     */
    private static boolean containsTrueAltAllele(final List<VariantContext> VCs) {

        for ( final VariantContext vc : VCs ) {
            if ( vc.getNAlleles() > 2 ) {
                return true;
            }
        }
        return false;
    }

    private static boolean annotationShouldBeSkippedForHomRefSites(VariantAnnotation annotation) {
        return annotation instanceof RankSumTest || annotation instanceof RMSMappingQuality || annotation instanceof AS_RMSMappingQuality;
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features, boolean includeNonVariants) {
        Utils.nonNull(originalVC);

        final VariantContext result;

        if ( originalVC.isVariant()  && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
            // only re-genotype polymorphic sites
            final VariantContext regenotypedVC = calculateGenotypes(originalVC);
            if (regenotypedVC == null || (!isProperlyPolymorphic(regenotypedVC) && !includeNonVariants)) {
                return null;
            }
            if (isProperlyPolymorphic(regenotypedVC) || includeNonVariants) {
                // Note that reversetrimAlleles must be performed after the annotations are finalized because the reducible annotation data maps
                // were generated and keyed on the un reverseTrimmed alleles from the starting VariantContexts. Thus reversing the order will make
                // it difficult to recover the data mapping due to the keyed alleles no longer being present in the variant context.
                final VariantContext withGenotypingAnnotations = addGenotypingAnnotations(originalVC.getAttributes(), regenotypedVC);
                final VariantContext withAnnotations = annotationEngine.finalizeAnnotations(withGenotypingAnnotations, originalVC);
                result = GATKVariantContextUtils.reverseTrimAlleles(withAnnotations);
            } else if (includeNonVariants) {
                result = originalVC;
            } else {
                return null;
            }
        } else {
            result = originalVC;
        }


        // if it turned monomorphic then we either need to ignore or fix such sites
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if (result.isPolymorphicInSamples()) {
            // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
            final VariantContext reannotated = annotationEngine.annotateContext(result, features, ref, null, a -> true);
            return new VariantContextBuilder(reannotated).genotypes(cleanupGenotypeAnnotations(reannotated, false)).make();
        } else if (includeNonVariants) {
            // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
            VariantContext reannotated = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            reannotated = annotationEngine.annotateContext(reannotated, features, ref, null, GVCFsSyncGenotyper::annotationShouldBeSkippedForHomRefSites);
            return removeNonRefAlleles(reannotated);
        } else {
            return null;
        }
    }

    private VariantContext calculateGenotypes(VariantContext vc){
        /*
         * Query the VariantContext for the appropriate model.  If type == MIXED, one would want to use model = BOTH.
         * However GenotypingEngine.getAlleleFrequencyPriors throws an exception if you give it anything but a SNP or INDEL model.
         */
        final GenotypeLikelihoodsCalculationModel model = vc.getType() == VariantContext.Type.INDEL
                ? GenotypeLikelihoodsCalculationModel.INDEL
                : GenotypeLikelihoodsCalculationModel.SNP;
        return genotypingEngine.calculateGenotypes(vc, model, null);
    }

    /**
     * Remove NON-REF alleles from the variant context
     *
     * @param vc   the variant context
     * @return variant context with the NON-REF alleles removed if multiallelic or replaced with NO-CALL alleles if biallelic
     */
    private VariantContext removeNonRefAlleles(final VariantContext vc) {

        // If NON_REF is the only alt allele, ignore this site
        final List<Allele> newAlleles = new ArrayList<>();
        // Only keep alleles that are not NON-REF
        for ( final Allele allele : vc.getAlleles() ) {
            if ( !allele.equals(Allele.NON_REF_ALLELE) ) {
                newAlleles.add(allele);
            }
        }

        // If no alt allele, then remove INFO fields that require alt alleles
        if ( newAlleles.size() == 1 ) {
            final VariantContextBuilder builder = new VariantContextBuilder(vc).alleles(newAlleles);
            // for ( final String name : infoHeaderAltAllelesLineNames ) {
            //     builder.rmAttributes(Arrays.asList(name));
            // }
            return builder.make();
        } else {
            return vc;
        }
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles.
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    public static boolean isProperlyPolymorphic(final VariantContext vc) {
        //obvious cases
        if (vc == null || vc.getAlternateAlleles().isEmpty()) {
            return false;
        } else if (vc.isBiallelic()) {
            return !(GATKVCFConstants.isSpanningDeletion(vc.getAlternateAllele(0)) || vc.isSymbolic());
        } else {
            return true;
        }
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations and any other annotations generated by the genotyper.
        final Map<String, Object> attrs = new LinkedHashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY)) {
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        }
        if (newVC.hasAttribute(GATKVCFConstants.AS_QUAL_KEY)) {
            attrs.put(GATKVCFConstants.AS_QUAL_KEY, newVC.getAttribute(GATKVCFConstants.AS_QUAL_KEY));
        }
        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }


    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param vc            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    @VisibleForTesting
    static List<Genotype> cleanupGenotypeAnnotations(final VariantContext vc, final boolean createRefGTs) {
        final GenotypesContext oldGTs = vc.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ) {
                depth = parseInt(oldGT.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
                builder.DP(depth);
                attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            }

            attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, PHASED_HOM_VAR_STRING);
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && vc.isVariant() ) {
                final int[] AD = new int[vc.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                // move the GQ to RGQ
                if (oldGT.hasGQ()) {
                    builder.noGQ();
                    attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
                }

                //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
                    final List<Allele> refAlleles = Collections.nCopies(oldGT.getPloidy(), vc.getReference());
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    private static int parseInt(Object attribute){
        if( attribute instanceof String) {
            return Integer.parseInt((String)attribute);
        } else if ( attribute instanceof Number){
            return ((Number) attribute).intValue();
        } else {
            throw new IllegalArgumentException("Expected a Number or a String but found something else.");
        }
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);

        //whether to emit non-variant sites is not contained in genotypeArgs and must be passed to uac separately
        //Note: GATK3 uses OutputMode.EMIT_ALL_CONFIDENT_SITES when includeNonVariants is requested
        //GATK4 uses EMIT_ALL_SITES to ensure LowQual sites are emitted.
        uac.outputMode = includeNonVariants ? OutputMode.EMIT_ALL_SITES : OutputMode.EMIT_VARIANTS_ONLY;
        return uac;
    }

    private List<Annotation> makeVariantAnnotations() {
        GATKAnnotationArgumentCollection userArgs = new DefaultGATKVariantAnnotationArgumentCollection();
        GATKAnnotationPluginDescriptor pluginDescriptor = new GATKAnnotationPluginDescriptor(userArgs, Collections.emptyList(), Arrays.asList(StandardAnnotation.class, AS_StandardAnnotation.class));
        findPluginsForDescriptor(pluginDescriptor);
        pluginDescriptor.validateAndResolvePlugins();
        return pluginDescriptor.getResolvedInstances();
    }

    private void findPluginsForDescriptor(
        GATKAnnotationPluginDescriptor pluginDescriptor) {
        final ClassFinder classFinder = new ClassFinder();
        pluginDescriptor.getPackageNames().forEach(
                pkg -> classFinder.find(pkg, pluginDescriptor.getPluginBaseClass()));
        final Set<Class<?>> pluginClasses = classFinder.getClasses();

        final List<Object> plugins = new ArrayList<>(pluginClasses.size());
        for (Class<?> c : pluginClasses) {
            if (pluginDescriptor.includePluginClass(c)) {
                try {
                    final Object plugin = pluginDescriptor.createInstanceForPlugin(c);
                    plugins.add(plugin);
                } catch (InstantiationException | IllegalAccessException e) {
                    throw new CommandLineException.CommandLineParserInternalException("Problem making an instance of plugin " + c +
                            " Do check that the class has a non-arg constructor", e);
                }
            }
        }
    }
}
