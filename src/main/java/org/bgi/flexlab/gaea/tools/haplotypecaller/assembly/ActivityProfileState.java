package org.bgi.flexlab.gaea.tools.haplotypecaller.assembly;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.Locatable;

public final class ActivityProfileState {
    private final GenomeLocation loc;
    private double activeProb;
    private final Type resultState;
    private final Number resultValue;

    public double isActiveProb() {
        return activeProb;
    }

    /**
     * Set the probability that this site is active.
     *
     * Probabilities should be between 0.0 and 1.0, however this is not currently enforced
     * because the {@link BandPassActivityProfile} can sometimes generate probabilities that
     * slightly exceed 1.0 when moving probability mass around. We intend to fix this by
     * capping at 1.0, but first we must evaluate the effects of capping on the HaplotypeCaller.
     *
     * @param activeProb probability (should be between 0.0 and 1.0) that the site is active
     */
    public void setIsActiveProb( final double activeProb ) {
        this.activeProb = activeProb;
    }

    /**
     * @return The type of the value returned by {@link #getResultValue}
     */
    public Type getResultState() {
        return resultState;
    }

    /**
     * @return Numeric value associated with {@link #getResultState}. If {@link #getResultState} is HIGH_QUALITY_SOFT_CLIPS,
     *         this is the number of bp affected by the soft clips
     */
    public Number getResultValue() {
        return resultValue;
    }

    /**
     * The type of the value returned by {@link #getResultValue}
     */
    public enum Type {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLocation loc, final double activeProb) {
        this(loc, activeProb, Type.NONE, null);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HIGH_QUALITY_SOFT_CLIPS, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final GenomeLocation loc, final double activeProb, final Type resultState, final Number resultValue) {
        if ( loc.size() != 1 ) {
            throw new IllegalArgumentException("Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        }
        if ( resultValue != null && resultValue.doubleValue() < 0 ) {
            throw new IllegalArgumentException("Result value isn't null and its < 0, which is illegal: " + resultValue);
        }

        this.loc = loc;
        setIsActiveProb(activeProb);
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a Locatable
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final Locatable regionStartLoc) {
        Utils.nonNull(regionStartLoc);
        return getLoc().getStart() - regionStartLoc.getStart();
    }

    /**
     * Get the locus associated with the ActivityProfileState
     * @return the locus associated with the ActivityProfileState as a GenomeLocation
     */
    public GenomeLocation getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileState{" +
                "loc=" + loc +
                ", activeProb=" + activeProb +
                ", resultState=" + resultState +
                ", resultValue=" + resultValue +
                '}';
    }
}
