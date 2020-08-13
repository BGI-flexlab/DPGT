package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;
import org.bgi.flexlab.gaea.tools.jointcalling.util.HeterozygosityUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ExcessHet extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation{

	private final double minNeededValue = 1.0E-16;
    private Set<String> founderIds;
    private final boolean RETURN_ROUNDED = true;
    private int sampleCount = -1;
    
    @Override
    public void initialize ( Set<VCFHeaderLine>headerLines,Set<String> sampleList ) {
        if(founderIds == null && sampleList != null) {
            founderIds = sampleList;
        }

    }
    
	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		return makeEHAnnotation(vc);
	}
	
	protected double calculateEH(final VariantContext vc, final GenotypesContext genotypes) {
        HeterozygosityUtils heterozygosityUtils = new HeterozygosityUtils(RETURN_ROUNDED);
        final double[] genotypeCountsDoubles = heterozygosityUtils.getGenotypeCountsForRefVsAllAlts(vc, genotypes);
        sampleCount = heterozygosityUtils.getSampleCount();
        final int[] genotypeCounts = new int[genotypeCountsDoubles.length];
        for(int i = 0; i < genotypeCountsDoubles.length; i++) {
            genotypeCounts[i] = (int)genotypeCountsDoubles[i];
        }

        double pval = exactTest(genotypeCounts);

        //If the actual phredPval would be infinity we will probably still filter out just a very large number
        if (pval == 0) {
        	return -10.0 * Math.log10(minNeededValue);
        }
        double phredPval = -10.0 * Math.log10(pval);

        return phredPval;
    }
	
	protected Map<String, Object> makeEHAnnotation(final VariantContext vc) {
        final GenotypesContext genotypes = (founderIds == null || founderIds.isEmpty()) ? vc.getGenotypes() : vc.getGenotypes(founderIds);
        if (genotypes == null || !vc.isVariant())
            return null;
        double EH = calculateEH(vc, genotypes);
        if (sampleCount < 1)
            return null;
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", EH));
    }
	
	protected double exactTest(final int[] genotypeCounts) {
        if (genotypeCounts.length != 3) {
            throw new IllegalStateException("Input genotype counts must be length 3 for the number of genotypes with {2, 1, 0} ref alleles.");
        }
        final int REF_INDEX = 0;
        final int HET_INDEX = 1;
        final int VAR_INDEX = 2;

        final int refCount = genotypeCounts[REF_INDEX];
        final int hetCount = genotypeCounts[HET_INDEX];
        final int homCount = genotypeCounts[VAR_INDEX];

        if (hetCount < 0 || refCount < 0 || homCount < 0) {
            throw new IllegalArgumentException("Genotype counts cannot be less than 0");
        }

        //Split into observed common allele and rare allele
        final int obsHomR;
        final int obsHomC;
        if (refCount < homCount) {
            obsHomR = refCount;
            obsHomC = homCount;
        } else {
            obsHomR = homCount;
            obsHomC = refCount;
        }

        final int rareCopies = 2 * obsHomR + hetCount;
        final int N = hetCount + obsHomC + obsHomR;

        //If the probability distribution has only 1 point, then the mid p-value is .5
        if (rareCopies <= 1) {
            return .5;
        }

        double[] probs = new double[rareCopies + 1];

        //Find (something close to the) mode for the midpoint
        int mid = (int) Math.floor(((double) rareCopies * (2.0 * (double) N - (double) rareCopies)) / (2.0 * (double) N - 1.0));
        if ((mid % 2) != (rareCopies % 2)) {
            mid++;
        }

        probs[mid] = 1.0;
        double mysum = 1.0;

        //Calculate probabilities from midpoint down
        int currHets = mid;
        int currHomR = (rareCopies - mid) / 2;
        int currHomC = N - currHets - currHomR;

        while (currHets >= 2) {
            double potentialProb = probs[currHets] * (double) currHets * ((double) currHets - 1.0) / (4.0 * ((double) currHomR + 1.0) * ((double) currHomC + 1.0));
            if (potentialProb < minNeededValue) {
                break;
            }

            probs[currHets - 2] = potentialProb;
            mysum = mysum + probs[currHets - 2];

            //2 fewer hets means one additional homR and homC each
            currHets = currHets - 2;
            currHomR = currHomR + 1;
            currHomC = currHomC + 1;
        }

        //Calculate probabilities from midpoint up
        currHets = mid;
        currHomR = (rareCopies - mid) / 2;
        currHomC = N - currHets - currHomR;

        while (currHets <= rareCopies - 2) {
            double potentialProb = probs[currHets] * 4.0 * (double) currHomR * (double) currHomC / (((double) currHets + 2.0) * ((double) currHets + 1.0));
            if (potentialProb < minNeededValue) {
                break;
            }

            probs[currHets + 2] = potentialProb;
            mysum = mysum + probs[currHets + 2];

            //2 more hets means 1 fewer homR and homC each
            currHets = currHets + 2;
            currHomR = currHomR - 1;
            currHomC = currHomC - 1;
        }

        double rightPval = probs[hetCount] / (2.0 * mysum);
        //Check if we observed the highest possible number of hets
        if (hetCount == rareCopies) {
            return rightPval;
        }
        rightPval = rightPval + StatUtils.sum(Arrays.copyOfRange(probs, hetCount + 1, probs.length)) / mysum;

        return (rightPval);
    }
	
	@Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.EXCESS_HET_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

}
