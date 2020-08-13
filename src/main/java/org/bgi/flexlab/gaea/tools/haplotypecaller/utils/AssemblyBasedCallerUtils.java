package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.Haplotype;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReferenceConfidenceModel;
import org.bgi.flexlab.gaea.tools.haplotypecaller.SampleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.HaplotypeCallerArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.LikelihoodEngineArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection.ReadThreadingAssemblerArgumentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.AssemblyRegion;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.AssemblyResultSet;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.ReadErrorCorrector;
import org.bgi.flexlab.gaea.tools.haplotypecaller.assembly.ReadThreadingAssembler;
import org.bgi.flexlab.gaea.tools.haplotypecaller.engine.HaplotypeCallerEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.engine.PairHMMLikelihoodCalculationEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.engine.RandomLikelihoodCalculationEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.engine.ReadLikelihoodCalculationEngine;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.FragmentCollection;
import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SmithWatermanAligner;
import org.bgi.flexlab.gaea.tools.haplotypecaller.writer.HaplotypeBAMWriter;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.util.QualityUtils;
import org.bgi.flexlab.gaea.util.ReadUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;

public final class AssemblyBasedCallerUtils {

	static final int REFERENCE_PADDING_FOR_ASSEMBLY = 500;

	/**
	 * Returns a map with the original read as a key and the realigned read as
	 * the value.
	 * <p>
	 * Missing keys or equivalent key and value pairs mean that the read was not
	 * realigned.
	 * </p>
	 * 
	 * @return never {@code null}
	 */
	public static Map<GaeaSamRecord, GaeaSamRecord> realignReadsToTheirBestHaplotype(
			final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype,
			final Locatable paddedReferenceLoc, final SmithWatermanAligner aligner) {
		final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
		final Map<GaeaSamRecord, GaeaSamRecord> result = new HashMap<>(bestAlleles.size());

		for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
			final GaeaSamRecord originalRead = bestAllele.read;
			final Haplotype bestHaplotype = bestAllele.allele;
			final boolean isInformative = bestAllele.isInformative();
			final GaeaSamRecord realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype,
					refHaplotype, paddedReferenceLoc.getStart(), isInformative, aligner);
			result.put(originalRead, realignedRead);
		}
		return result;
	}

	public static void finalizeRegion(final AssemblyRegion region, final boolean errorCorrectReads,
			final boolean dontUseSoftClippedBases, final byte minTailQuality, final SAMFileHeader readsHeader,
			final SampleList samplesList) {
		if (region.isFinalized()) {
			return;
		}

		// Loop through the reads hard clipping the adaptor and low quality
		// tails
		final List<GaeaSamRecord> readsToUse = new ArrayList<>(region.getReads().size());
		for (final GaeaSamRecord myRead : region.getReads()) {
			final byte minTailQualityToUse = errorCorrectReads
					? HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION : minTailQuality;
			GaeaSamRecord clippedRead = ReadClipper.hardClipLowQualEnds(myRead, minTailQualityToUse);

			// remove soft clips if we cannot reliably clip off adapter sequence
			// or if the user doesn't want to use soft clips at all
			// otherwie revert soft clips so that we see the alignment start and
			// end assuming the soft clips are all matches
			// TODO -- WARNING -- still possibility that unclipping the soft
			// clips will introduce bases that aren't
			// TODO -- truly in the extended region, as the unclipped bases
			// might actually include a deletion
			// TODO -- w.r.t. the reference. What really needs to happen is that
			// kmers that occur before the
			// TODO -- reference haplotype start must be removed
			clippedRead = dontUseSoftClippedBases || !ReadUtils.hasWellDefinedFragmentSize(clippedRead)
					? ReadClipper.hardClipSoftClippedBases(clippedRead)
					: ReadClipper.revertSoftClippedBases(clippedRead);

			clippedRead = clippedRead.isUnmapped() ? clippedRead : ReadClipper.hardClipAdaptorSequence(clippedRead);
			if (!clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0) {
				clippedRead = ReadClipper.hardClipToRegion(clippedRead, region.getExtendedSpan().getStart(),
						region.getExtendedSpan().getEnd());
				if (region.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0) {
					readsToUse.add(clippedRead);
				}
			}
		}

		// TODO -- Performance optimization: we partition the reads by sample 4
		// times right now; let's unify that code.
		// final List<GaeaSamRecord> downsampledReads =
		// DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse),
		// maxReadsInRegionPerSample, minReadsPerAlignmentStart);
		Collections.sort(readsToUse, new ReadCoordinateComparator(readsHeader)); 

		// handle overlapping read pairs from the same fragment
		cleanOverlappingReadPairs(readsToUse, samplesList, readsHeader);

		region.clearReads();
		region.addAll(readsToUse);
		region.setFinalized(true);
	}

	/**
	 * Clean up reads/bases that overlap within read pairs
	 *
	 * @param reads
	 *            the list of reads to consider
	 */
	private static void cleanOverlappingReadPairs(final List<GaeaSamRecord> reads, final SampleList samplesList,
			final SAMFileHeader readsHeader) {
		for (final List<GaeaSamRecord> perSampleReadList : splitReadsBySample(samplesList, readsHeader, reads)
				.values()) {
			final FragmentCollection<GaeaSamRecord> fragmentCollection = FragmentCollection.create(perSampleReadList);
			for (final List<GaeaSamRecord> overlappingPair : fragmentCollection.getOverlappingPairs()) {
				FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
			}
		}
	}

	public static Map<String, List<GaeaSamRecord>> splitReadsBySample(final SampleList samplesList,
			final SAMFileHeader header, final Collection<GaeaSamRecord> reads) {
		final Map<String, List<GaeaSamRecord>> returnMap = new HashMap<>();
		for (final String sample : samplesList.asListOfSamples()) {
			returnMap.put(sample, new ArrayList<>());
		}

		for (final GaeaSamRecord read : reads) {
			returnMap.get(ReadUtils.getSampleName(read, header)).add(read);
		}

		return returnMap;
	}

	/**
	 * Helper function to create the reference haplotype out of the active
	 * region and a padded loc
	 * 
	 * @param region
	 *            the active region from which to generate the reference
	 *            haplotype
	 * @param paddedReferenceLoc
	 *            the interval which includes padding and shows how big the
	 *            reference haplotype should be
	 * @return a non-null haplotype
	 */
	public static Haplotype createReferenceHaplotype(final AssemblyRegion region,
			final GenomeLocation paddedReferenceLoc, final ChromosomeInformationShare reference) {
		return ReferenceConfidenceModel.createReferenceHaplotype(region,
				region.getAssemblyRegionReference(reference), paddedReferenceLoc);
	}

	public static GenomeLocation getPaddedReferenceLoc(final AssemblyRegion region, final int referencePadding,
			final ChromosomeInformationShare reference) {
		final int padLeft = Math.max(region.getExtendedSpan().getStart() - referencePadding, 1);
		final int padRight = Math.min(region.getExtendedSpan().getEnd() + referencePadding, reference.getLength());
		return new GenomeLocation(region.getExtendedSpan().getContig(), padLeft, padRight);
	}

	/**
	 * Instantiates the appropriate likelihood calculation engine.
	 *
	 * @return never {@code null}.
	 */
	public static ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine(
			final LikelihoodEngineArgumentCollection likelihoodArgs) {
		final double log10GlobalReadMismappingRate = likelihoodArgs.phredScaledGlobalReadMismappingRate < 0
				? -Double.MAX_VALUE
				: QualityUtils.qualToErrorProbLog10(likelihoodArgs.phredScaledGlobalReadMismappingRate);

		switch (likelihoodArgs.likelihoodEngineImplementation) {
		case PairHMM:
			return new PairHMMLikelihoodCalculationEngine((byte) likelihoodArgs.gcpHMM,
					likelihoodArgs.pairHMMNativeArgs.getPairHMMArgs(), likelihoodArgs.pairHMM,
					log10GlobalReadMismappingRate, likelihoodArgs.pcrErrorModel,
					likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD);
		case Random:
			return new RandomLikelihoodCalculationEngine();
		default:
			throw new UserException("Unsupported likelihood calculation engine.");
		}
	}

	public static ReadThreadingAssembler createReadThreadingAssembler(
			final HaplotypeCallerArgumentCollection args) {
		final ReadThreadingAssemblerArgumentCollection rtaac = args.assemblerArgs;
		final ReadThreadingAssembler assemblyEngine = new ReadThreadingAssembler(rtaac.maxNumHaplotypesInPopulation,
				rtaac.kmerSizes, rtaac.dontIncreaseKmerSizesForCycles, rtaac.allowNonUniqueKmersInRef,
				rtaac.numPruningSamples);
		assemblyEngine.setErrorCorrectKmers(rtaac.errorCorrectKmers);
		assemblyEngine.setPruneFactor(rtaac.minPruneFactor);
		assemblyEngine.setDebug(args.debug);
		assemblyEngine.setDebugGraphTransformations(rtaac.debugGraphTransformations);
		assemblyEngine.setRecoverDanglingBranches(!rtaac.doNotRecoverDanglingBranches);
		assemblyEngine.setMinDanglingBranchLength(rtaac.minDanglingBranchLength);
		assemblyEngine.setMinBaseQualityToUseInAssembly(args.minBaseQualityScore);

		if (rtaac.graphOutput != null) {
			assemblyEngine.setGraphWriter(new File(rtaac.graphOutput));
		}

		return assemblyEngine;
	}

	public static Optional<HaplotypeBAMWriter> createBamWriter(final HaplotypeCallerArgumentCollection args,
			final boolean createBamOutIndex, final boolean createBamOutMD5, final SAMFileHeader header) {
		return args.bamOutputPath != null ? Optional.of(HaplotypeBAMWriter.create(args.bamWriterType,
				new File(args.bamOutputPath), createBamOutIndex, createBamOutMD5, header)) : Optional.empty();
	}

	// create the assembly using just high quality reads (eg Q20 or higher). We
	// may want to use lower
	// quality reads in the PairHMM downstream, so we can't use a ReadFilter
	public static AssemblyRegion assemblyRegionWithWellMappedReads(final AssemblyRegion originalAssemblyRegion,
			final int minMappingQuality, final SAMFileHeader readsHeader) {
		final AssemblyRegion result = new AssemblyRegion(originalAssemblyRegion.getSpan(),
				originalAssemblyRegion.getSupportingStates(), originalAssemblyRegion.isActive(),
				originalAssemblyRegion.getExtension(), readsHeader);
		originalAssemblyRegion.getReads().stream().filter(rec -> rec.getMappingQuality() >= minMappingQuality)
				.forEach(result::add);
		return result;
	}

	// Contract: the List<Allele> alleles of the resulting VariantContext is the
	// ref allele followed by alt alleles in the
	// same order as in the input vcs
	public static VariantContext makeMergedVariantContext(final List<VariantContext> vcs) {
		if (vcs.isEmpty()) {
			return null;
		}
		final List<String> haplotypeSources = vcs.stream().map(VariantContext::getSource).collect(Collectors.toList());
		return GaeaGvcfVariantContextUtils.simpleMerge(vcs, haplotypeSources,
				GaeaGvcfVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
				GaeaGvcfVariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, null, false, false);
	}

	/**
	 * High-level function that runs the assembler on the given region's reads,
	 * returning a data structure with the resulting information needed for
	 * further HC steps
	 */
	public static AssemblyResultSet assembleReads(final AssemblyRegion region, final List<VariantContext> givenAlleles,
			final HaplotypeCallerArgumentCollection argumentCollection, final SAMFileHeader header,
			final SampleList sampleList, final ChromosomeInformationShare referenceReader,
			final ReadThreadingAssembler assemblyEngine, final SmithWatermanAligner aligner,final int maxDepthForAssembly) {
		finalizeRegion(region, argumentCollection.errorCorrectReads, argumentCollection.dontUseSoftClippedBases,
				(byte) (argumentCollection.minBaseQualityScore - 1), header, sampleList);

		final byte[] fullReferenceWithPadding = region.getAssemblyRegionReference(referenceReader,
				REFERENCE_PADDING_FOR_ASSEMBLY);
		final GenomeLocation paddedReferenceLoc = getPaddedReferenceLoc(region, REFERENCE_PADDING_FOR_ASSEMBLY,
				referenceReader);
		final Haplotype referenceHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, referenceReader);

		final ReadErrorCorrector readErrorCorrector = argumentCollection.errorCorrectReads
				? new ReadErrorCorrector(argumentCollection.assemblerArgs.kmerLengthForReadErrorCorrection,
						HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
						argumentCollection.assemblerArgs.minObservationsForKmerToBeSolid, argumentCollection.debug,
						fullReferenceWithPadding)
				: null;

		try {
			final AssemblyResultSet assemblyResultSet = assemblyEngine.runLocalAssembly(region, referenceHaplotype,
					fullReferenceWithPadding, paddedReferenceLoc, givenAlleles, readErrorCorrector, header, aligner,maxDepthForAssembly);
			return assemblyResultSet;
		} catch (final Exception e) {
			// Capture any exception that might be thrown, and write out the
			// assembly failure BAM if requested
			if (argumentCollection.captureAssemblyFailureBAM) {
				try (final SAMFileWriter writer = ReadUtils.createCommonSAMWriter(new File("assemblyFailure.bam"), null,
						header, false, false, false)) {
					for (final GaeaSamRecord read : region.getReads()) {
						writer.addAlignment(read);
					}
				}
			}
			throw e;
		}
	}
}
