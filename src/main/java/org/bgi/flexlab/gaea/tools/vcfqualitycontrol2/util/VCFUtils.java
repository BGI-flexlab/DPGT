package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.GaeaVariantContextUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class VCFUtils {
	
	enum VCFIndexType {
	    DYNAMIC_SEEK,       // use DynamicIndexCreator(IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
	    DYNAMIC_SIZE,       // use DynamicIndexCreator(IndexFactory.IndexBalanceApproach.FOR_SIZE)
	    LINEAR,             // use LinearIndexCreator()
	    INTERVAL            // use IntervalIndexCreator()
	}
	
	public static final String VCF_FILE_EXTENSION = "vcf";
    public static final String BCF_FILE_EXTENSION = "bcf";

    // as determined experimentally Nov-Dec 2013
    public final static VCFIndexType DEFAULT_GVCF_INDEX_TYPE = VCFIndexType.LINEAR;
    public final static Integer DEFAULT_GVCF_INDEX_PARAMETER = 128000;

    private VCFUtils(){}

    public static SortedSet<String> getSortedSampleSet(Map<String, VCFHeader> headers, GaeaVariantContextUtils.GenotypeMergeType mergeOption) {
        final SortedSet<String> samples = new TreeSet<>();
        for (final Map.Entry<String, VCFHeader> val : headers.entrySet()) {
            VCFHeader header = val.getValue();
            samples.addAll(header.getGenotypeSamples().stream().map(sample -> GaeaVariantContextUtils.mergedSampleName(val.getKey(), sample,
                    mergeOption == GaeaVariantContextUtils.GenotypeMergeType.UNIQUIFY)).collect(Collectors.toList()));
        }

        return samples;
    }


    /**
     * Check if using the GCVF indexing arguments' values
     *
     * @param variantIndexType variant indexing strategy
     * @param variantIndexParameter variant indexing parameter
     * @return true if the index type and parameter are the default GVCF values, false otherwise
     */
    public static boolean usingGVCFIndexingArguments(final VCFIndexType variantIndexType, final int variantIndexParameter) {
        return variantIndexType == DEFAULT_GVCF_INDEX_TYPE && variantIndexParameter == DEFAULT_GVCF_INDEX_PARAMETER;
    }

    //TODO: these should be refactored/consolidated as part of
    // https://github.com/broadinstitute/gatk/issues/121 and
    // https://github.com/broadinstitute/gatk/issues/1116
    /**
     * Given a set of VCF header lines, update the set with contig
     * lines from the provided reference dictionary.
     * @param oldLines
     * @param referenceFile
     * @param refDict
     * @param referenceNameOnly
     * @return Updated list of VCF header lines.
     */
    public static Set<VCFHeaderLine> updateHeaderContigLines(
            final Set<VCFHeaderLine> oldLines,
            final File referenceFile,
            final SAMSequenceDictionary refDict,
            final boolean referenceNameOnly) {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>(oldLines.size());

        for (final VCFHeaderLine line : oldLines) {
            if (line instanceof VCFContigHeaderLine) {
                continue; // skip old contig lines
            }
            if (line.getKey().equals(VCFHeader.REFERENCE_KEY)) {
                continue; // skip the old reference key
            }
            lines.add(line);
        }

        lines.addAll(makeContigHeaderLines(refDict, referenceFile).stream().collect(Collectors.toList()));

        if (referenceFile != null) {
            final String referenceValue;
            if (referenceNameOnly) {
                final int extensionStart = referenceFile.getName().lastIndexOf(".");
                referenceValue = extensionStart == -1 ? referenceFile.getName() : referenceFile.getName().substring(0, extensionStart);
            }
            else {
                referenceValue = "file://" + referenceFile.getAbsolutePath();
            }
            lines.add(new VCFHeaderLine(VCFHeader.REFERENCE_KEY, referenceValue));
        }
        return lines;
    }

    private static List<VCFContigHeaderLine> makeContigHeaderLines(final SAMSequenceDictionary refDict,
                                                                   final File referenceFile) {
        final List<VCFContigHeaderLine> lines = new ArrayList<>();
        final String assembly = referenceFile != null ? referenceFile.getName() : null;
        lines.addAll(refDict.getSequences().stream().map(contig -> makeContigHeaderLine(contig, assembly)).collect(Collectors.toList()));
        return lines;
    }

    private static VCFContigHeaderLine makeContigHeaderLine(final SAMSequenceRecord contig, final String assembly) {
        final Map<String, String> map = new LinkedHashMap<>(3);
        map.put(GaeaVCFConstants.CONTIG_ID_KEY, contig.getSequenceName());
        map.put(GaeaVCFConstants.CONTIG_LENGTH_KEY, String.valueOf(contig.getSequenceLength()));
        if (assembly != null) {
            map.put(GaeaVCFConstants.ASSEMBLY_NAME_KEY, assembly);
        }
        return new VCFContigHeaderLine(map, contig.getSequenceIndex());
    }
    
    public static void addVQSRStandardHeaderLines(final Set<VCFHeaderLine> hInfo) {
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.VQS_LOD_KEY));
        hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.CULPRIT_KEY));
        hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.NEGATIVE_LABEL_KEY));
    }
}
