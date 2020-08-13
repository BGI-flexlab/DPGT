package org.bgi.flexlab.gaea.tools.haplotypecaller.writer;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.GaeaVCFHeaderLines;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class GVCFWriter implements VariantContextWriter {

    /** Where we'll ultimately write our VCF records */
	private static final int MAX_GENOTYPE_QUAL = VCFConstants.MAX_GENOTYPE_QUAL;

    //
    // Final fields initialized in constructor
    //
    /** Where we'll ultimately write our VCF records */
    private final VariantContextWriter underlyingWriter;

    private final List<HomRefBlock> GQPartitions;

    /** fields updated on the fly during GVCFWriter operation */
    int nextAvailableStart = -1;
    String contigOfNextAvailableStart = null;
    private String sampleName = null;
    private HomRefBlock currentBlock = null;
    private final int defaultPloidy;

    /**
     * Are the proposed GQ partitions well-formed?
     *
     * @param GQPartitions proposed GQ partitions
     * @return a non-null string if something is wrong (string explains issue)
     */
    protected static List<HomRefBlock> parsePartitions(final List<Integer> GQPartitions, final int defaultPloidy) {
        if ( GQPartitions == null ) {
            throw new IllegalArgumentException("The list of GQ partitions cannot be null.");
        }
        if ( GQPartitions.isEmpty() ) {
            throw new IllegalArgumentException("The list of GQ partitions cannot be empty.");
        }

        final List<HomRefBlock> result = new LinkedList<>();
        int lastThreshold = 0;
        for ( final Integer value : GQPartitions ) {
            if ( value == null || value <= 0 ) {
                throw new IllegalArgumentException("The list of GQ partitions contains a null or non-positive integer.");
            }
            if ( value < lastThreshold ) {
                throw new IllegalArgumentException(String.format("The list of GQ partitions is out of order. " +
                        "Previous value is %d but the next is %d.", lastThreshold, value));
            }
            if ( value == lastThreshold ) {
                throw new IllegalArgumentException(String.format("The value %d appears more than once in the list of GQ partitions.", value));
            }
            if ( value > MAX_GENOTYPE_QUAL + 1 ) {
                throw new IllegalArgumentException(String.format("The value %d in the list of GQ partitions is " +
                        "greater than VCFConstants.MAX_GENOTYPE_QUAL + 1 = %d.", value, VCFConstants.MAX_GENOTYPE_QUAL + 1));
            }
            result.add(new HomRefBlock(lastThreshold, value, defaultPloidy));
            lastThreshold = value;
        }
        if ( lastThreshold <= MAX_GENOTYPE_QUAL ) {
            result.add(new HomRefBlock(lastThreshold, MAX_GENOTYPE_QUAL + 1, defaultPloidy));
        }
        return result;
    }

    /**
     * Create a new GVCF writer
     *
     * Should be a non-empty list of boundaries.  For example, suppose this variable is
     *
     * [A, B, C]
     *
     * We would partition our hom-ref sites into the following bands:
     *
     * X < A
     * A <= X < B
     * B <= X < C
     * X >= C
     *
     * @param underlyingWriter the ultimate destination of the GVCF records
     * @param GQPartitions a well-formed list of GQ partitions
     * @param defaultPloidy the assumed ploidy for input variant context without one.
     */
    public GVCFWriter(final VariantContextWriter underlyingWriter, final List<Integer> GQPartitions, final int defaultPloidy) {
        if ( underlyingWriter == null ) throw new IllegalArgumentException("underlyingWriter cannot be null");
        this.underlyingWriter = underlyingWriter;
        this.GQPartitions = parsePartitions(GQPartitions,defaultPloidy);
        this.defaultPloidy = defaultPloidy;
    }

    /**
     * Write the VCF header
     *
     * Adds standard GVCF fields to the header
     *
     * @param header a non-null header
     */
    @Override
    public void writeHeader(VCFHeader header) {
        if ( header == null ) throw new IllegalArgumentException("header cannot be null");
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GaeaVCFHeaderLines.getFormatLine(GaeaVCFConstants.MIN_DP_FORMAT_KEY));

        for ( final HomRefBlock partition : GQPartitions ) {
            header.addMetaDataLine(partition.toVCFHeaderLine());
        }

        underlyingWriter.writeHeader(header);
    }

    /**
     * Close this GVCF writer.  Finalizes any pending hom-ref blocks and emits those to the underlyingWriter as well
     */
    @Override
    public void close() {
        close(true);
    }

    /**
     * Horrible work around because there's no clean way to get our VCFWriter closed by the GATK
     *
     * If closeUnderlyingWriter is true, then we'll close the underlying writer, otherwise we'll leave it open
     * so the GATK closes it later
     *
     * @param closeUnderlyingWriter should we leave the underlying writer open or closed?
     */
    public void close(final boolean closeUnderlyingWriter) {
        emitCurrentBlock();
        if ( closeUnderlyingWriter ) underlyingWriter.close();
    }

    /**
     * Add hom-ref site from vc to this gVCF hom-ref state tracking, emitting any pending states if appropriate
     *
     * @param vc a non-null VariantContext
     * @param g a non-null genotype from VariantContext
     * @return a VariantContext to be emitted, or null if non is appropriate
     */
    protected VariantContext addHomRefSite(final VariantContext vc, final Genotype g) {

        if ( nextAvailableStart != -1 ) {
            // don't create blocks while the hom-ref site falls before nextAvailableStart (for deletions)
            if ( vc.getStart() <= nextAvailableStart && vc.getChr().equals(contigOfNextAvailableStart) )
                return null;
            // otherwise, reset to non-relevant
            nextAvailableStart = -1;
            contigOfNextAvailableStart = null;
        }

        final VariantContext result;
        if (genotypeCanBeMergedInCurrentBlock(g)) {
            currentBlock.add(vc.getStart(), g);
            result = null;
        } else {
            result = blockToVCF(currentBlock);
            currentBlock = createNewBlock(vc, g);
        }
        return result;
    }

    private boolean genotypeCanBeMergedInCurrentBlock(final Genotype g) {
        return currentBlock != null && currentBlock.withinBounds(capToMaxGQ(g.getGQ())) && currentBlock.getPloidy() == g.getPloidy()
                && (currentBlock.getMinPLs() == null || !g.hasPL() || (currentBlock.getMinPLs().length == g.getPL().length));
    }

    private int capToMaxGQ(final int gq) {
        return Math.min(gq, MAX_GENOTYPE_QUAL);
    }

    /**
     * Flush the current hom-ref block, if necessary, to the underlying writer, and reset the currentBlock to null
     */
    private void emitCurrentBlock() {
        if ( currentBlock != null ) {
            // there's actually some work to do
            underlyingWriter.add(blockToVCF(currentBlock));
            currentBlock = null;
        }
    }

    /**
     * Convert a HomRefBlock into a VariantContext
     *
     * @param block the block to convert
     * @return a VariantContext representing the gVCF encoding for this block.
     * It will return {@code null} if input {@code block} is {@code null}, indicating that there
     * is no variant-context to be output into the VCF.
     */
    private VariantContext blockToVCF(final HomRefBlock block) {
        if ( block == null ) return null;

        final VariantContextBuilder vcb = new VariantContextBuilder(block.getStartingVC());
        vcb.attributes(new HashMap<String, Object>(2)); // clear the attributes
        vcb.stop(block.getStop());
        vcb.attribute(VCFConstants.END_KEY, block.getStop());

        // create the single Genotype with GQ and DP annotations
        final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GaeaGvcfVariantContextUtils.homozygousAlleleList(block.getRef(),block.getPloidy()));
        gb.noAD().noPL().noAttributes(); // clear all attributes

        final int[] minPLs = block.getMinPLs();
        gb.PL(minPLs);
        final int gq = GaeaGvcfVariantContextUtils.calculateGQFromPLs(minPLs);
        gb.GQ(gq);
        gb.DP(block.getMedianDP());
        gb.attribute(GaeaVCFConstants.MIN_DP_FORMAT_KEY, block.getMinDP());

        // This annotation is no longer standard
        //gb.attribute(MIN_GQ_FORMAT_FIELD, block.getMinGQ());

        return vcb.genotypes(gb.make()).make();
    }

    /**
     * Helper function to create a new HomRefBlock from a variant context and current genotype
     *
     * @param vc the VariantContext at the site where want to start the band
     * @param g the genotype of the sample from vc that should be used to initialize the block
     * @return a newly allocated and initialized block containing g already
     */
    private HomRefBlock createNewBlock(final VariantContext vc, final Genotype g) {
        // figure out the GQ limits to use based on the GQ of g
        HomRefBlock partition = null;
        for ( final HomRefBlock maybePartition : GQPartitions ) {
            if ( maybePartition.withinBounds(capToMaxGQ(g.getGQ())) ) {
                partition = maybePartition;
                break;
            }
        }

        if ( partition == null )
            throw new IllegalStateException("GQ " + g + " from " + vc + " didn't fit into any partition");

        // create the block, add g to it, and return it for use
        final HomRefBlock block = new HomRefBlock(vc, partition.getGQLowerBound(), partition.getGQUpperBound(), defaultPloidy);
        block.add(vc.getStart(), g);
        return block;
    }

    /**
     * Add a VariantContext to this writer for emission
     *
     * Requires that the VC have exactly one genotype
     *
     * @param vc a non-null VariantContext
     */
    @Override
    public void add(VariantContext vc) {
        if ( vc == null ) throw new IllegalArgumentException("vc cannot be null");

        if ( sampleName == null )
            sampleName = vc.getGenotype(0).getSampleName();

        if ( ! vc.hasGenotypes() ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has genotypes");
        } else if ( vc.getGenotypes().size() != 1 ) {
            throw new IllegalArgumentException("GVCF assumes that the VariantContext has exactly one genotype but saw " + vc.getGenotypes().size());
        } else {
            if ( currentBlock != null && ! currentBlock.isContiguous(vc) ) {
                // we've made a non-contiguous step (across interval, onto another chr), so finalize
                emitCurrentBlock();
            }

            final Genotype g = vc.getGenotype(0);
            if ( g.isHomRef() && vc.hasAlternateAllele(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE) && vc.isBiallelic() ) {
                // create bands
                final VariantContext maybeCompletedBand = addHomRefSite(vc, g);
                if ( maybeCompletedBand != null ) underlyingWriter.add(maybeCompletedBand);
            } else {
                // g is variant, so flush the bands and emit vc
                emitCurrentBlock();
                nextAvailableStart = vc.getEnd();
                contigOfNextAvailableStart = vc.getChr();
                underlyingWriter.add(vc);
            }
        }
    }

    /**
     * Check the return from PrintStream.checkError() if underlying stream is a java.io.PrintStream
     * @return false, no error since the underlying stream is not a java.io.PrintStream
     */
    public boolean checkError(){
        return false;
    }
}

