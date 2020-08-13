package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACcounts;
import org.bgi.flexlab.gaea.tools.genotyer.genotypecaller.ExactACset;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GvcfMathUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.util.MathUtils;

import static java.lang.Math.pow;

public abstract class DiploidExactAFCalculator extends ExactAFCalculator {
	private static final double LOG10_OF_2 = GvcfMathUtils.Log10Cache.get(2);

    public DiploidExactAFCalculator() {
    }

    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                               final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        final int numAlternateAlleles = vc.getNAlleles() - 1;

        final ArrayList<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true, vc.hasAllele(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE));
        //raw code
//        final int numSamples = genotypeLikelihoods.size()-1;
//        final int numChr = 2*numSamples;
//
//        // queue of AC conformations to process
//        final LinkedList<ExactACset> ACqueue = new LinkedList<>();
//
//        // mapping of ExactACset indexes to the objects
//        final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<>(numChr+1);
//
//        // add AC=0 to the queue
//        final int[] zeroCounts = new int[numAlternateAlleles];
//        ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
//        ACqueue.add(zeroSet);
//        indexesToACset.put(zeroSet.getACcounts(), zeroSet);
//
//        while ( !ACqueue.isEmpty() ) {
//            stateTracker.incNEvaluations(); // keep track of the number of evaluations
//
//            // compute log10Likelihoods
//            final ExactACset set = ACqueue.remove();
//
//            calculateAlleleCountConformation(set, genotypeLikelihoods, numChr, ACqueue,
//                    indexesToACset, log10AlleleFrequencyPriors,stateTracker);
//
//            // clean up memory
//            indexesToACset.remove(set.getACcounts());
//            //if ( DEBUG )
//            //    System.out.printf(" *** removing used set=%s%n", set.ACcounts);
//        }
        //raw code ends
        //EM algorithm

        final int numSamples3 = genotypeLikelihoods.size()-1;
        final int numChr3 = 2*numSamples3;
        double af=0.5;
        int iter=0;
        double lastAf=0;
        for(int j=1;j<=numSamples3;j++) {
            double[] gls = genotypeLikelihoods.get(j);
            for (int k = 0; k < 3; k++) {
                gls[k] = pow(10, gls[k]);
            }
        }
        for(int i=1;i<=100;i++){
            double totalSum=0;
            for(int j=1;j<=numSamples3;j++) {
                double[] gls = genotypeLikelihoods.get(j);
                double fenzi=0+1*gls[1]*2*af*(1-af)+2*gls[2]*af*af;
                double fenmu=gls[0]*(1-af)*(1-af)+gls[1]*2*af*(1-af)+gls[2]*af*af;
                totalSum+=fenzi/fenmu;
            }
            af=totalSum/numChr3;
            if(Math.abs(af-lastAf)<1.0/(numChr3*100) || af<1.0/(numChr3*5)){
                iter=i;
                break;
            }
            lastAf=af;
        }
        double dAC=af*numChr3;
        int ac=(int)(dAC+0.5d);
        int acCutoff=100;
        final ArrayList<double[]> genotypeLikelihoods3 = getGLs(vc.getGenotypes(), true, vc.hasAllele(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE));
        final LinkedList<ExactACset> ACqueue3 = new LinkedList<>();

        final HashMap<ExactACcounts, ExactACset> indexesToACset3 = new HashMap<>(numChr3+1);

        final int[] zeroCounts3 = new int[numAlternateAlleles];
        ExactACset zeroSet3 = new ExactACset(numSamples3+1, new ExactACcounts(zeroCounts3));
        ACqueue3.add(zeroSet3);
        indexesToACset3.put(zeroSet3.getACcounts(), zeroSet3);
        while ( !ACqueue3.isEmpty() ) {
            stateTracker.incNEvaluations();
            final ExactACset set = ACqueue3.remove();
            calculateAlleleCountConformation(set, genotypeLikelihoods3, numChr3, ACqueue3,
                    indexesToACset3, log10AlleleFrequencyPriors,stateTracker);
            indexesToACset3.remove(set.getACcounts());
            if(stateTracker.getNEvaluations()>=10 && ac>acCutoff){
                break;
            }
        }
        if(!ACqueue3.isEmpty()) {
            int[] updateAC = new int[1];
            updateAC[0] = stateTracker.getNEvaluations() - 2;
            int startAC=updateAC[0];
            if (ac > acCutoff) {
                double curMle=stateTracker.log10MLE;
                double virtualMLE = curMle > -2.0 ? curMle * 0.5 : -2.0;
                double interval = (virtualMLE - curMle) / (ac - startAC);

                for(int i=0;i<ac-startAC;i++) {
                    updateAC[0]++;
                    double newMle = stateTracker.log10MLE + interval;
                    stateTracker.updateMLEifNeeded(newMle, updateAC);
                    newMle += log10AlleleFrequencyPriors[updateAC[0]];
                    stateTracker.updateMAPifNeeded(newMle, updateAC);
                    stateTracker.incNEvaluations();
                    if(i==ac-startAC-1){
                        for(int j=0;j<3;j++) {
                            curMle=stateTracker.log10MLE;
                            newMle = curMle - 1.5 * interval;
                            updateAC[0]++;
                            if(updateAC[0]>numChr3){
                                break;
                            }
                            stateTracker.updateMLEifNeeded(newMle, updateAC);
                            newMle += log10AlleleFrequencyPriors[updateAC[0]];
                            stateTracker.updateMAPifNeeded(newMle, updateAC);
                            stateTracker.incNEvaluations();
                        }
                    }
                }
            }
        }

        //EM code ends

        //filter sample code start
        /*
        if(ac==-100) {
            final ArrayList<double[]> genotypeLikelihoods2 = getGLs2(vc.getGenotypes(), true, vc.hasAllele(GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE));
            final int numSamples = genotypeLikelihoods.size() - 1;
            final int numChr = 2 * numSamples;

            //reduced
            final int numSamples2 = genotypeLikelihoods2.size() - 1;
            final int numChr2 = 2 * numSamples2;
            // queue of AC conformations to process
            final LinkedList<ExactACset> ACqueue = new LinkedList<>();
            final LinkedList<ExactACset> ACqueue2 = new LinkedList<>();

            // mapping of ExactACset indexes to the objects
            final HashMap<ExactACcounts, ExactACset> indexesToACset = new HashMap<>(numChr + 1);
            final HashMap<ExactACcounts, ExactACset> indexesToACset2 = new HashMap<>(numChr + 1);

            // add AC=0 to the queue
            final int[] zeroCounts = new int[numAlternateAlleles];
            ExactACset zeroSet = new ExactACset(numSamples + 1, new ExactACcounts(zeroCounts));
            // reduced
            ExactACset zeroSet2 = new ExactACset(numSamples2 + 1, new ExactACcounts(zeroCounts));
            ACqueue.add(zeroSet);
            ACqueue2.add(zeroSet2);
            indexesToACset.put(zeroSet.getACcounts(), zeroSet);
            indexesToACset2.put(zeroSet2.getACcounts(), zeroSet2);
            StateTracker stateTracker2 = new StateTracker(stateTracker.getAlleleCountsOfMAP().length);


            // compute log10Likelihoods
            ExactACset set2 = ACqueue.remove();
            computeLofK(set2, genotypeLikelihoods, log10AlleleFrequencyPriors, stateTracker2);
            final double log10LofK = set2.getLog10Likelihoods()[set2.getLog10Likelihoods().length - 1];
            indexesToACset.remove(set2.getACcounts());
            indexesToACset.clear();
            set2 = null;
            genotypeLikelihoods.clear();
            stateTracker2 = null;
            // clean up memory

            //reduced
            while (!ACqueue2.isEmpty()) {
                stateTracker.incNEvaluations(); // keep track of the number of evaluations

                // compute log10Likelihoods
                final ExactACset set = ACqueue2.remove();

                calculateAlleleCountConformation(set, genotypeLikelihoods2, numChr2, ACqueue2,
                        indexesToACset2, log10AlleleFrequencyPriors, stateTracker);

                // clean up memory
                indexesToACset2.remove(set.getACcounts());
                //if ( DEBUG )
                //    System.out.printf(" *** removing used set=%s%n", set.ACcounts);
            }
            int addACs = 0;
            for (HashMap.Entry<Integer, Integer> kv : ExactAFCalculator.excludedACs.entrySet()) {
                if (kv.getKey() == 1) {
                    addACs += kv.getValue();
                } else if (kv.getKey() == 2) {
                    addACs += 2 * kv.getValue();
                }
            }

            stateTracker.setLog10LikelihoodOfAFzero(log10LofK);
            stateTracker.setLog10PosteriorOfAFzero(log10LofK + log10AlleleFrequencyPriors[0]);
            double curMLE = log10LofK;
            int[] updateAC = new int[1];
            updateAC[0] = (int) MathUtils.sum(stateTracker.getAlleleCountsOfMAP()) + addACs;
            stateTracker.updateMLEifNeeded(stateTracker.log10MLE + 1, updateAC);
            excludedACs.clear();
        }
        */
        //filter sample code ends
        return getResultFromFinalState(vc, log10AlleleFrequencyPriors, stateTracker);
    }


    @Override
    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        return GaeaGvcfVariantContextUtils.subsetAlleles(vc, allelesToUse, GaeaGvcfVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }

    @Override
    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        final ArrayList<double[]> GLs = getGLs(vc.getGenotypes(), true);
        for ( final double[] likelihoods : GLs ) {
            final int PLindexOfBestGL = GvcfMathUtils.maxElementIndex(likelihoods);
            if ( PLindexOfBestGL != PL_INDEX_OF_HOM_REF ) {
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindexOfBestGL);
                final int alleleLikelihoodIndex1 = alleles.alleleIndex1 - 1;
                final int alleleLikelihoodIndex2 = alleles.alleleIndex2 - 1;
                if ( alleles.alleleIndex1 != 0 )
                    likelihoodSums[alleleLikelihoodIndex1].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
                // don't double-count it
                if ( alleles.alleleIndex2 != 0 && alleles.alleleIndex2 != alleles.alleleIndex1 )
                    likelihoodSums[alleleLikelihoodIndex2].sum += likelihoods[PLindexOfBestGL] - likelihoods[PL_INDEX_OF_HOM_REF];
            }
        }
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;

        public DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }

    private double calculateAlleleCountConformation(final ExactACset set,
                                                    final ArrayList<double[]> genotypeLikelihoods,
                                                    final int numChr,
                                                    final LinkedList<ExactACset> ACqueue,
                                                    final HashMap<ExactACcounts, ExactACset> indexesToACset,
                                                    final double[] log10AlleleFrequencyPriors,
                                                    final StateTracker stateTracker) {

        //if ( DEBUG )
        //    System.out.printf(" *** computing LofK for set=%s%n", set.ACcounts);

        // compute the log10Likelihoods
        computeLofK(set, genotypeLikelihoods, log10AlleleFrequencyPriors, stateTracker);

        final double log10LofK = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];

        // can we abort early because the log10Likelihoods are so small?
        if ( stateTracker.abort(log10LofK, set.getACcounts(), true, false) ) {
            //if ( DEBUG )
            //    System.out.printf(" *** breaking early set=%s log10L=%.2f maxLog10L=%.2f%n", set.ACcounts, log10LofK, maxLog10L);
            return log10LofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ) // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return log10LofK;

        final int numAltAlleles = set.getACcounts().getCounts().length;

        // add conformations for the k+1 case
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // to get to this conformation, a sample would need to be AB (remember that ref=0)
            final int PLindex = GenotypeLikelihoods.calculatePLindex(0, allele+1);
            updateACset(ACcountsClone, numChr, set, PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final ArrayList<DependentSet> differentAlleles = new ArrayList<>(numAltAlleles * numAltAlleles);
            final ArrayList<DependentSet> sameAlleles = new ArrayList<>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.getACcounts().getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    // to get to this conformation, a sample would need to be BB or BC (remember that ref=0, so add one to the index)
                    final int PLindex = GenotypeLikelihoods.calculatePLindex(allele_i+1, allele_j+1);
                    if ( allele_i == allele_j )
                        sameAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    else
                        differentAlleles.add(new DependentSet(ACcountsClone, PLindex));
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( DependentSet dependent : differentAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            for ( DependentSet dependent : sameAlleles )
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        return log10LofK;
    }

    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also pushes its value to the given callingSetIndex.
    private void updateACset(final int[] newSetCounts,
                             final int numChr,
                             final ExactACset dependentSet,
                             final int PLsetIndex,
                             final Queue<ExactACset> ACqueue,
                             final HashMap<ExactACcounts, ExactACset> indexesToACset,
                             final ArrayList<double[]> genotypeLikelihoods) {
        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // push data from the dependency to the new set
        //if ( DEBUG )
        //    System.out.println(" *** pushing data from " + index + " to " + dependencySet.ACcounts);
        pushData(indexesToACset.get(index), dependentSet, PLsetIndex, genotypeLikelihoods);
    }

    private void computeLofK(final ExactACset set,
                             final ArrayList<double[]> genotypeLikelihoods,
                             final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {

        set.getLog10Likelihoods()[0] = 0.0; // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for ( int j = 1; j < set.getLog10Likelihoods().length; j++ )
                set.getLog10Likelihoods()[j] = set.getLog10Likelihoods()[j-1] + genotypeLikelihoods.get(j)[HOM_REF_INDEX];

            final double log10Lof0 = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];
            stateTracker.setLog10LikelihoodOfAFzero(log10Lof0);
            stateTracker.setLog10PosteriorOfAFzero(log10Lof0 + log10AlleleFrequencyPriors[0]);
            return;
        }

        // if we got here, then k > 0 for at least one k.
        // the non-AA possible conformations were already dealt with by pushes from dependent sets;
        // now deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
        for ( int j = 1; j < set.getLog10Likelihoods().length; j++ ) {

            if ( totalK < 2*j-1 ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue = GvcfMathUtils.Log10Cache.get(2*j-totalK) + GvcfMathUtils.Log10Cache.get(2*j-totalK-1) + set.getLog10Likelihoods()[j-1] + gl[HOM_REF_INDEX];
                set.getLog10Likelihoods()[j] = GvcfMathUtils.approximateLog10SumLog10(set.getLog10Likelihoods()[j], conformationValue);
            }

            final double logDenominator = GvcfMathUtils.Log10Cache.get(2*j) + GvcfMathUtils.Log10Cache.get(2*j-1);
            set.getLog10Likelihoods()[j] = set.getLog10Likelihoods()[j] - logDenominator;
        }

        double log10LofK = set.getLog10Likelihoods()[set.getLog10Likelihoods().length-1];

        // update the MLE if necessary
        stateTracker.updateMLEifNeeded(log10LofK, set.getACcounts().getCounts());

        // apply the priors over each alternate allele
        for ( final int ACcount : set.getACcounts().getCounts() ) {
            if ( ACcount > 0 )
                log10LofK += log10AlleleFrequencyPriors[ACcount];
        }

        stateTracker.updateMAPifNeeded(log10LofK, set.getACcounts().getCounts());
    }

    private void pushData(final ExactACset targetSet,
                          final ExactACset dependentSet,
                          final int PLsetIndex,
                          final ArrayList<double[]> genotypeLikelihoods) {
        final int totalK = targetSet.getACsum();
        final double[] targetLog10Likelihoods = targetSet.getLog10Likelihoods();
        final double[] dependentLog10Likelihoods = dependentSet.getLog10Likelihoods();
        final int[] targetACcounts = targetSet.getACcounts().getCounts();

        // skip impossible conformations, namely those for which there aren't enough possible chromosomes (2*j in the if loop below)
        // to fill up the totalK alternate alleles needed; we do want to ensure that there's at least one sample included (hence the Math.max)
        final int firstIndex = Math.max(1, (totalK + 1) / 2);

        // find the 2 alleles that are represented by this PL index
        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLsetIndex);
        // if neither allele is reference then the coefficient is constant throughout this invocation
        final Double constantCoefficient = (alleles.alleleIndex1 == 0) ? null : determineCoefficient(alleles, 0, targetACcounts, totalK);

        for ( int j = firstIndex; j < targetLog10Likelihoods.length; j++ ) {
            final double[] gl = genotypeLikelihoods.get(j);
            final double coefficient = (constantCoefficient != null) ? constantCoefficient : determineCoefficient(alleles, j, targetACcounts, totalK);

            final double conformationValue = coefficient + dependentLog10Likelihoods[j-1] + gl[PLsetIndex];
            targetLog10Likelihoods[j] = GvcfMathUtils.approximateLog10SumLog10(targetLog10Likelihoods[j], conformationValue);
        }
    }

    private double determineCoefficient(final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles, final int j, final int[] ACcounts, final int totalK) {
        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        // *** note that throughout this method we subtract one from the alleleIndex because ACcounts ***
        // *** doesn't consider the reference allele whereas the GenotypeLikelihoods PL cache does.   ***

        // the AX het case
        if ( alleles.alleleIndex1 == 0 )
            return GvcfMathUtils.Log10Cache.get(2*ACcounts[alleles.alleleIndex2-1]) + GvcfMathUtils.Log10Cache.get(2*j-totalK);

        final int k_i = ACcounts[alleles.alleleIndex1-1];

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles.alleleIndex1 == alleles.alleleIndex2 ) {
            coeff = GvcfMathUtils.Log10Cache.get(k_i) + GvcfMathUtils.Log10Cache.get(k_i - 1);
        }
        // the het non-ref case (e.g. BC, BD, CD)
        else {
            final int k_j = ACcounts[alleles.alleleIndex2-1];
            coeff = LOG10_OF_2 + GvcfMathUtils.Log10Cache.get(k_i) + GvcfMathUtils.Log10Cache.get(k_j);
        }

        return coeff;
    }

    @Override
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                          final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        if (defaultPloidy != 2)
            throw new IllegalArgumentException("cannot support ploidy different than 2 and the default ploidy is " + defaultPloidy);
        return allelesToUse.size() == 1
                ? GaeaGvcfVariantContextUtils.subsetToRefOnly(vc, 2)
                : GaeaGvcfVariantContextUtils.subsetAlleles(vc, allelesToUse,
                     assignGenotypes ? GaeaGvcfVariantContextUtils.GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN : GaeaGvcfVariantContextUtils.GenotypeAssignmentMethod.SET_TO_NO_CALL);
    }
}
