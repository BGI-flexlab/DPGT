package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;

import org.apache.commons.lang3.tuple.Pair;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.StandardCallerArgumentCollection;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;

public final class AssemblyRegionTrimmer {

    /**
     * Holds the extension to be used based on whether GGA mode is on or off.
     */
    private int usableExtension;

    /**
     * Records whether the trimming intervals are going to be used to emit reference confidence, {@code true},
     * or regular HC output {@code false}.
     */
    private boolean emitReferenceConfidence;

    private StandardCallerArgumentCollection standardArgs;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     * @param assemblyRegionTrimmerArgs user arguments for the trimmer
     * @param sequenceDictionary dictionary to determine the bounds of contigs
     * @param debug whether to show extra debug log messages.
     * @param isGGA whether the trimming region calculator should act as if we are in GGA mode or not.
     * @param emitReferenceConfidence indicates whether we plan to use this trimmer to generate trimmed regions
     *                                to be used for emitting reference confidence.
     *
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws CommandLineException.BadArgumentValue if any of the user argument values is invalid.
     */
    public void initialize(final StandardCallerArgumentCollection assemblyRegionTrimmerArgs, final SAMSequenceDictionary sequenceDictionary, final boolean debug, final boolean isGGA, final boolean emitReferenceConfidence) {
        if (this.standardArgs != null) {
            throw new IllegalStateException(getClass().getSimpleName() + " instance initialized twice");
        }

        this.standardArgs = Utils.nonNull(assemblyRegionTrimmerArgs);;
        this.sequenceDictionary = sequenceDictionary;

        checkUserArguments();
        usableExtension = isGGA ? this.standardArgs.ggaExtension : this.standardArgs.discoverExtension;
        this.emitReferenceConfidence = emitReferenceConfidence;
    }

    /**
     * Checks user trimming argument values
     *
     * @throws CommandLineException.BadArgumentValue if there is some problem with any of the arguments values.
     */
    private void checkUserArguments() {
        if ( standardArgs.snpPadding < 0 ) {
            throw new UserException.BadArgumentValueException("paddingAroundSNPs", "" + standardArgs.snpPadding + "< 0");
        }
        if ( standardArgs.indelPadding < 0 ) {
            throw new UserException.BadArgumentValueException("paddingAroundIndels", "" + standardArgs.indelPadding + "< 0");
        }
        if ( standardArgs.discoverExtension < 0) {
            throw new UserException.BadArgumentValueException("maxDiscARExtension", "" + standardArgs.discoverExtension + "< 0");
        }
        if ( standardArgs.ggaExtension < 0) {
            throw new UserException.BadArgumentValueException("maxGGAAREExtension", "" + standardArgs.ggaExtension + "< 0");
        }
    }

    /**
     * Holds the result of trimming.
     */
    public static final class Result {

        /**
         * Indicates whether trimming is required per data and user request.
         */
        protected final boolean needsTrimming;

        /**
         * Holds the input active region.
         */
        protected final AssemblyRegion originalRegion;

        /**
         * Holds the smaller range that contain all relevant callable variants in the
         * input active region (not considering the extension).
         *
         */
        protected final GenomeLocation callableSpan;

        /**
         * Maximum available range for the trimmed variant region.
         */
        protected final GenomeLocation maximumSpan;

        /**
         * The trimmed variant region span including the extension.
         */
        protected final GenomeLocation extendedSpan;


        /**
         * The ideal trimmer variant region span including the extension.
         */
        protected final GenomeLocation idealSpan;

        /**
         * Returns the ideal trimming span.
         *
         * <p/>
         * The ideal span is the one containing all callable variation overlapping the original active region span
         * (without extension) and the applicable padding {@link #getPadding()} in both sides.
         *
         *
         * @return never {@code null}.
         */
        public GenomeLocation getIdealSpan() {
            return idealSpan;
        }

        /**
         * Holds the flanking spans that do not contain the callable variants.
         * <p/>
         * The first element of the pair is the left (up-stream) non-variant flank, whereas the second element is
         * the right (down-stream) non-variant flank.
         */
        protected final Pair<GenomeLocation, GenomeLocation> nonVariantFlanks;

        /**
         * Holds the collection of callable events within the variant trimming region.
         */
        protected final List<VariantContext> callableEvents;

        /**
         * Required padding around the variant trimming region.
         */
        protected final int padding;


        /**
         * Returns the required padding around callable variation.
         *
         * <p/>
         * Notice that due to the limiting span of the original active region (including its extension) it
         * is possible that the resulting final trimmed variant region span does not satisfies the padding. However
         * that should be rare.
         *
         * @return 0 or greater.
         */
        public int getPadding() {
            return padding;
        }

        /**
         * Holds the maximum extension around the original active region span considered for the trimmed
         * variation region.
         */
        protected final int usableExtension;

        /**
         * Returns the maximum extension around the original active region span considered for the trimmed
         * variation region.
         *
         * <p/>
         * From time to time, the trimmed region may require a span beyond the input original active region's.
         * For example when there is a callable event close ot one of its ends and the required padding makes it
         * round beyond that limit.
         *
         * <p/>
         * Notice that due to the limiting span of the original active region (including its extended region) it
         * is possible that the resulting final trimmed variant region span goes beyond this extension including more of
         * the original active region own extension.
         *
         * @return 0 or greater.
         */
        public int getUsableExtension() {
            return usableExtension;
        }

        /**
         * Holds variant-containing callable region.
         * <p/>
         * This is lazy-initialized using {@link #callableSpan}.
         */
        protected AssemblyRegion callableRegion;


        /**
         * Non-variant left flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getFirst()}.
         */
        private AssemblyRegion leftFlankRegion;

        /**
         * Non-variant right flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getSecond()}.
         */
        private AssemblyRegion rightFlankRegion;

        /**
         * Whether the variant trimmed region is going to be used for emitting reference confidence records.
         */
        private final boolean emitReferenceConfidence;

        /**
         * Creates a trimming result given all its properties.
         *
         * @param emitReferenceConfidence whether reference confidence output modes are on.
         * @param needsTrimming whether there is any trimming needed at all.
         * @param originalRegion the original active region.
         * @param padding padding around contained callable variation events.
         * @param extension the extension applied to the trimmed variant span.
         * @param overlappingEvents contained callable variation events.
         * @param nonVariantFlanks pair of non-variant flank spans around the variant containing span.
         * @param extendedSpan final trimmed variant span including the extension.
         * @param idealSpan the ideal span, that contains.
         * @param maximumSpan maximum possible trimmed span based on the input original active region extended span.
         * @param callableSpan variant containing span without padding.
         */
        protected Result(final boolean emitReferenceConfidence,
                         final boolean needsTrimming,
                         final AssemblyRegion originalRegion,
                         final int padding,
                         final int extension,
                         final List<VariantContext> overlappingEvents,
                         final Pair<GenomeLocation, GenomeLocation> nonVariantFlanks,
                         final GenomeLocation extendedSpan,
                         final GenomeLocation idealSpan,
                         final GenomeLocation maximumSpan,
                         final GenomeLocation callableSpan) {
            this.emitReferenceConfidence = emitReferenceConfidence;
            this.needsTrimming = needsTrimming;
            this.originalRegion = originalRegion;
            this.nonVariantFlanks = nonVariantFlanks;
            this.padding = padding;
            usableExtension = extension;
            callableEvents = overlappingEvents;
            this.callableSpan = callableSpan;
            this.idealSpan = idealSpan;
            this.maximumSpan = maximumSpan;
            this.extendedSpan = extendedSpan;

            Utils.validateArg(extendedSpan == null || callableSpan == null || extendedSpan.contains(callableSpan), "the extended callable span must include the callable span");
        }


        /**
         * Checks whether there is any variation present in the target region.
         *
         * @return {@code true} if there is any variant, {@code false} otherwise.
         */
        public boolean isVariationPresent() {
            return ! callableEvents.isEmpty();
        }

        /**
         * Checks whether the active region needs trimming.
         */
        public boolean needsTrimming() {
            return needsTrimming;
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getCallableRegion() {
            if (callableRegion == null && extendedSpan != null) {
                //TODO this conditional is a patch to retain the current standard HC run behaviour
                //TODO we should simply remove this difference between trimming with or without GVCF
                //TODO embracing slight changes in the standard HC output
                callableRegion = emitReferenceConfidence ? originalRegion.trim(callableSpan, extendedSpan) : originalRegion.trim(extendedSpan);
            } else if (extendedSpan == null) {
                throw new IllegalStateException("there is no variation thus no variant region");
            }
            return callableRegion;
        }

        /**
         * Checks whether there is a non-empty left flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial left flank region, {@code false} otherwise.
         */
        public boolean hasLeftFlankingRegion() {
            return nonVariantFlanks.getLeft() != null;
        }

        /**
         * Checks whether there is a non-empty right flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial right flank region, {@code false} otherwise.
         */
        public boolean hasRightFlankingRegion() {
            return nonVariantFlanks.getRight() != null;
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *  <p/>
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         *
         *  @throws IllegalStateException if there is not such as left flanking region.
         */
        public AssemblyRegion nonVariantLeftFlankRegion() {
            if (leftFlankRegion == null && nonVariantFlanks.getLeft() != null) {
                leftFlankRegion = originalRegion.trim(nonVariantFlanks.getLeft(), originalRegion.getExtension());
            } else if (nonVariantFlanks.getLeft() == null) {
                throw new IllegalStateException("there is no left flank non-variant trimmed out region");
            }
            return leftFlankRegion;
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public AssemblyRegion nonVariantRightFlankRegion() {
            if (rightFlankRegion == null && nonVariantFlanks.getRight() != null) {
                rightFlankRegion = originalRegion.trim(nonVariantFlanks.getRight(), originalRegion.getExtension());
            } else if (nonVariantFlanks.getRight() == null) {
                throw new IllegalStateException("there is no right flank non-variant trimmed out region");
            }
            return rightFlankRegion;
        }

        /**
         * Creates a result indicating that there was no trimming to be done.
         */
        protected static Result noTrimming(final boolean emitReferenceConfidence,
                                           final AssemblyRegion targetRegion, final int padding,
                                           final int usableExtension,final List<VariantContext> events) {
            final GenomeLocation targetRegionLoc = targetRegion.getSpan();
            final Result result = new Result(emitReferenceConfidence,false,targetRegion,padding,usableExtension,events,Pair.of(null, null),
                    targetRegionLoc,targetRegionLoc,targetRegionLoc,targetRegionLoc);
            result.callableRegion = targetRegion;
            return result;
        }

        /**
         * Creates a result indicating that no variation was found.
         */
        protected static Result noVariation(final boolean emitReferenceConfidence, final AssemblyRegion targetRegion,
                                            final int padding, final int usableExtension) {
            final Result result = new Result(emitReferenceConfidence,false,targetRegion,padding,usableExtension,
                    Collections.emptyList(), Pair.of(targetRegion.getSpan(), null),
                    null, null, null, null);
            result.leftFlankRegion = targetRegion;
            return result;
        }
    }

    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered latter.
     *
     * @param originalRegion the genome location range to trim.
     * @param allVariantsWithinExtendedRegion list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code originalRegion} are simply ignored.
     * @return never {@code null}.
     */
    public Result trim(final AssemblyRegion originalRegion,
                       final SortedSet<VariantContext> allVariantsWithinExtendedRegion) {


        if ( allVariantsWithinExtendedRegion.isEmpty() ) // no variants,
        {
            return Result.noVariation(emitReferenceConfidence, originalRegion, standardArgs.snpPadding, usableExtension);
        }

        final List<VariantContext> withinActiveRegion = new LinkedList<>();
        final GenomeLocation originalRegionRange = originalRegion.getSpan();
        boolean foundNonSnp = false;
        GenomeLocation variantSpan = null;
        for ( final VariantContext vc : allVariantsWithinExtendedRegion ) {
            final GenomeLocation vcLoc = new GenomeLocation(vc);
            if ( originalRegionRange.overlaps(vcLoc) ) {
                foundNonSnp = foundNonSnp || ! vc.isSNP();
                variantSpan = variantSpan == null ? vcLoc : variantSpan.spanWith(vcLoc);
                withinActiveRegion.add(vc);
            }
        }
        final int padding = foundNonSnp ? standardArgs.indelPadding : standardArgs.snpPadding;

        // we don't actually have anything in the region after skipping out variants that don't overlap
        // the region's full location
        if ( variantSpan == null ) {
            return Result.noVariation(emitReferenceConfidence, originalRegion, padding, usableExtension);
        }

        if ( standardArgs.dontTrimActiveRegions) {
            return Result.noTrimming(emitReferenceConfidence, originalRegion, padding, usableExtension, withinActiveRegion);
        }

        final GenomeLocation maximumSpan = originalRegionRange.expandWithinContig(usableExtension, sequenceDictionary);
        final GenomeLocation idealSpan = variantSpan.expandWithinContig(padding, sequenceDictionary);
        final GenomeLocation finalSpan = maximumSpan.intersect(idealSpan).mergeWithContiguous(variantSpan);

        // Make double sure that, if we are emitting GVCF we won't call non-variable positions beyond the target active region span.
        // In regular call we don't do so so we don't care and we want to maintain behavior, so the conditional.
        final GenomeLocation callableSpan = emitReferenceConfidence ? variantSpan.intersect(originalRegionRange) : variantSpan;

        final Pair<GenomeLocation, GenomeLocation> nonVariantRegions = nonVariantTargetRegions(originalRegion, callableSpan);

        return new Result(emitReferenceConfidence,true,originalRegion,padding, usableExtension,withinActiveRegion,nonVariantRegions,finalSpan,idealSpan,maximumSpan,variantSpan);
    }

    /**
     * Calculates the list of region to trim away.
     * @param targetRegion region for which to generate the flanking regions.
     * @param variantSpan the span of the core region containing relevant variation and required padding.
     * @return never {@code null}; 0, 1 or 2 element list.
     */
    private Pair<GenomeLocation, GenomeLocation> nonVariantTargetRegions(final AssemblyRegion targetRegion, final GenomeLocation variantSpan) {
        final GenomeLocation targetRegionRange = targetRegion.getSpan();
        final int finalStart = variantSpan.getStart();
        final int finalStop = variantSpan.getEnd();

        final int targetStart = targetRegionRange.getStart();
        final int targetStop = targetRegionRange.getEnd();

        final boolean preTrimmingRequired = targetStart < finalStart;
        final boolean postTrimmingRequired = targetStop > finalStop;
        if (preTrimmingRequired) {
            final String contig = targetRegionRange.getContig();
            return postTrimmingRequired ? Pair.of(new GenomeLocation(contig, targetStart, finalStart - 1), new GenomeLocation(contig, finalStop + 1, targetStop)) :
                                          Pair.of(new GenomeLocation(contig, targetStart, finalStart - 1), null);
        } else if (postTrimmingRequired) {
            return Pair.of(null, new GenomeLocation(targetRegionRange.getContig(), finalStop + 1, targetStop));
        } else {
            return Pair.of(null, null);
        }
    }
}
