package org.bgi.flexlab.gaea.tools.haplotypecaller.writer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;

final class HomRefBlock implements Locatable {

	private final VariantContext startingVC;
    private int stop;
    private final int minGQ, maxGQ;
    private int[] minPLs = null;
    final private List<Integer> GQs = new ArrayList<>(100);
    final private List<Integer> DPs = new ArrayList<>(100);
    private final Allele ref;
    private final int ploidy;

    /**
     * Create a new HomRefBlock
     *
     * @param startingVC the VariantContext that starts this band (for starting position information)
     * @param minGQ the minGQ (inclusive) to use in this band
     * @param maxGQ the maxGQ (exclusive) to use in this band
     */
    public HomRefBlock(final VariantContext startingVC, final int minGQ, final int maxGQ, final int defaultPloidy) {
        if ( startingVC == null ) throw new IllegalArgumentException("startingVC cannot be null");
        if ( minGQ > maxGQ ) throw new IllegalArgumentException("bad minGQ " + minGQ + " as its > maxGQ " + maxGQ);

        this.startingVC = startingVC;
        this.stop = getStart() - 1;
        this.ref = startingVC.getReference();
        this.minGQ = minGQ;
        this.maxGQ = maxGQ;
        this.ploidy = startingVC.getMaxPloidy(defaultPloidy);
    }

    /**
     * Create a new HomRefBlock only for doing bounds checking
     *
     * @param minGQ the minGQ (inclusive) to use in this band
     * @param maxGQ the maxGQ (exclusive) to use in this band
     */
    public HomRefBlock(final int minGQ, final int maxGQ, final int ploidy) {
        if ( minGQ > maxGQ ) throw new IllegalArgumentException("bad minGQ " + minGQ + " as its > maxGQ " + maxGQ);

        this.startingVC = null;
        this.stop = -1;
        this.ref = null;
        this.minGQ = minGQ;
        this.maxGQ = maxGQ;
        this.ploidy = ploidy;
    }

    /**
     * Add information from this Genotype to this band
     * @param g a non-null Genotype with GQ and DP attributes
     */
    public void add(final int pos, final Genotype g) {
        if ( g == null ) throw new IllegalArgumentException("g cannot be null");
        if ( ! g.hasGQ() ) throw new IllegalArgumentException("g must have GQ field");
        if ( ! g.hasPL() ) throw new IllegalArgumentException("g must have PL field");
        if ( pos != stop + 1 ) throw new IllegalArgumentException("adding genotype at pos " + pos + " isn't contiguous with previous stop " + stop);
        if ( g.getPloidy() != ploidy)
            throw new IllegalArgumentException("cannot add a genotype with a different ploidy: " + g.getPloidy() + " != " + ploidy);

        if( minPLs == null )
            minPLs = g.getPL();
        else { // otherwise take the min with the provided genotype's PLs
            final int[] PL = g.getPL();
            if (PL.length != minPLs.length)
                throw new IllegalStateException("trying to merge different PL array sizes: " + PL.length + " != " + minPLs.length);
            for (int i = 0; i < PL.length; i++)
                if (minPLs[i] > PL[i])
                    minPLs[i] = PL[i];
        }
        stop = pos;
        GQs.add(Math.min(g.getGQ(), 99)); // cap the GQs by the max. of 99 emission
        DPs.add(Math.max(g.getDP(),0));
    }

    /**
     * Is the GQ value within the bounds of this GQ (GQ >= minGQ && GQ < maxGQ)
     * @param GQ the GQ value to test
     * @return true if within bounds, false otherwise
     */
    public boolean withinBounds(final int GQ) {
        return GQ >= minGQ && GQ < maxGQ;
    }

    /** Get the min GQ observed within this band */
    /*public int getMinGQ() { return GvcfMathUtils.arrayMin(GQs); }

    public int getMedianGQ() { return GvcfMathUtils.median(GQs); }*/

    public int getMinDP() { return GvcfMathUtils.arrayMin(DPs); }

    public int getMedianDP() { return (int) Math.round(MathUtils.median(DPs)); }
    /** Get the min PLs observed within this band, can be null if no PLs have yet been observed */
    public int[] getMinPLs() { return minPLs; }

    protected int getGQUpperBound() { return maxGQ; }
    protected int getGQLowerBound() { return minGQ; }

    public boolean isContiguous(final VariantContext vc) {
        return vc.getEnd() == getStop() + 1 && startingVC.getContig().equals(vc.getContig());
    }

    public VariantContext getStartingVC() { return startingVC; }
    public int getStart() { return startingVC.getStart(); }
    public int getStop() { return stop; }
    public Allele getRef() { return ref; }
    public int getSize() { return getStop() - getStart() + 1; }

    @Override
    public String toString() {
        return "HomRefBlock{" +
                "minGQ=" + minGQ +
                ", maxGQ=" + maxGQ +
                '}';
    }

    public VCFHeaderLine toVCFHeaderLine() {
        // Need to uniquify the key for the header line using the min/max GQ, since
        // VCFHeader does not allow lines with duplicate keys.
        final String key = String.format("GVCFBlock%d-%d", getGQLowerBound(), getGQUpperBound());
        return new VCFHeaderLine(key, "minGQ=" + getGQLowerBound() + "(inclusive),maxGQ=" + getGQUpperBound() + "(exclusive)");
    }

    /**
     * Get the ploidy of this hom-ref block.
     * @return
     */
    public int getPloidy() {
        return ploidy;
    }

	@Override
	public String getContig() {
		return startingVC.getContig();
	}

	@Override
	public int getEnd() {
		return getStop();
	}
}
