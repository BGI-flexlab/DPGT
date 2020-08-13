package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.Haplotype;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class GenotypingGivenAllelesUtils {

    /**
     * Composes the given allele variant-context providing information about the rods and reference location.
     * @param tracker the meta data tracker.
     * @param loc the query location.
     * @param snpsOnly whether we only should consider SNP variation.
     * @param logger where to output warnings.
     * @param allelesBinding the target variation context binding containing the given alleles.
     * @return never {@code null}
     */
    public static VariantContext composeGivenAllelesVariantContextFromRod(final RefMetaDataTracker tracker,
                                                                          final Locatable loc,
                                                                          final boolean snpsOnly,
                                                                          final List<VariantContext> allelesBinding) {
        Utils.nonNull(tracker, "tracker may not be null");
        Utils.nonNull(loc, "location may not be null");
        Utils.nonNull(allelesBinding, "alleles binding may not be null");
        VariantContext vc = null;

        // search for usable record
        for ( final VariantContext rodVc : tracker.getValues(allelesBinding, new GenomeLocation(loc)) ) {
            if ( rodVc != null && ! rodVc.isFiltered() && (! snpsOnly || rodVc.isSNP() )) {
                if ( vc == null ) {
                    vc = rodVc;
                }
            }
        }

        return vc;
    }

    /**
     * Create the list of artificial GGA-mode haplotypes by injecting each of the provided alternate alleles into the reference haplotype
     *
     * @param refHaplotype the reference haplotype
     * @param givenHaplotypes the list of alternate alleles in VariantContexts
     * @param activeRegionWindow the window containing the reference haplotype
     *
     * @return a non-null list of haplotypes
     */
    public static List<Haplotype> composeGivenHaplotypes(final Haplotype refHaplotype, final List<VariantContext> givenHaplotypes, final GenomeLocation activeRegionWindow) {
        Utils.nonNull(refHaplotype, "reference haplotype may not be null");
        Utils.nonNull(givenHaplotypes, "given haplotypes may not be null");
        Utils.nonNull(activeRegionWindow, "active region window may not be null");
        Utils.validateArg(activeRegionWindow.size() == refHaplotype.length(), "inconsistent reference haplotype and active region window");

        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();
        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();

        for( final VariantContext compVC : givenHaplotypes ) {
            Utils.validateArg(GaeaGvcfVariantContextUtils.overlapsRegion(compVC, activeRegionWindow),
                    " some variant provided does not overlap with active region window");
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                if( insertedRefHaplotype != null ) { // can be null if the requested allele can't be inserted into the haplotype
                    returnHaplotypes.add(insertedRefHaplotype);
                }
            }
        }

        return new ArrayList<>(returnHaplotypes);
    }
}
