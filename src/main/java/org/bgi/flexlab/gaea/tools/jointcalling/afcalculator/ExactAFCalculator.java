package org.bgi.flexlab.gaea.tools.jointcalling.afcalculator;

import java.util.*;

import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.MathUtils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

abstract class ExactAFCalculator extends AFCalculator {
    public static int excludedSamples=0;
    public static Map<Integer,Integer> excludedACs=new HashMap<Integer,Integer>();
    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    // useful so that we don't keep printing out the same warning messages
    protected static boolean printedMaxAltAllelesWarning = false;

    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {

        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return - Double.compare(o1.sum,o2.sum);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with higher likelihood are first but make sure that
     * NON_REF alleles are place are last.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            if (o1.allele == GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE)
                return 1;
            else if (o2.allele == GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE)
                return -1;
            else
                return o1.compareTo(o2);
        }
    };
    /**
     * Sorts {@link ExactAFCalculator.LikelihoodSum} instances where those with lower alternative allele index are first regardless of
     * the likelihood sum.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_INDEX_COMPARATOR = new Comparator<LikelihoodSum>() {
        @Override
        public int compare(final LikelihoodSum o1, final LikelihoodSum o2) {
            return Integer.compare(o1.index, o2.index);
        }
    };

    protected ExactAFCalculator() {

    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public final Allele allele;
        public final int index;

        public LikelihoodSum(final Allele allele, final int index) { this.allele = allele; this.index = index; }

        public int compareTo(LikelihoodSum other) {
            final double diff = sum - other.sum;
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        return getGLs(GLs, includeDummy, false);
    }

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @param keepUninformative Don't filter out uninformative genotype likelihoods (i.e. all log likelihoods near 0)
     *                          This is useful for VariantContexts with a NON_REF allele
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy, final boolean keepUninformative) {
        final ArrayList<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();
                
                if ( MathUtils.sum(gls) < GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL || keepUninformative )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }
    protected static ArrayList<double[]> getGLs2(final GenotypesContext GLs, final boolean includeDummy, final boolean keepUninformative) {
        final ArrayList<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        int haveAC=0;
        int totalNumber=0;
        int nonAC=0;
        int sampleSize=0;

        int highAC=0;
        int lowAC=0;
        int highRef=0;
        int lowRef=0;
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ){
            sampleSize++;
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();
                double maxGT=-1000;
                int maxIndex=-1;
                for(int i=0;i<gls.length;i++){
                    if(gls[i]>maxGT){
                        maxGT=gls[i];
                        maxIndex=i;
                    }
                }
                if ( MathUtils.sum(gls) < GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL || keepUninformative ){
                    int gq=sample.getGQ();
                    if(sample.getGQ()==-1){
                        int maxNumber=0;
                        double secondMax=-1000;
                        for(int i=0;i<gls.length;i++){
                            if(gls[i]==maxGT){
                                maxNumber++;
                            }else if(gls[i]>secondMax){
                                secondMax=gls[i];
                            }
                        }
                        gq=(int)(maxGT-secondMax)*10;
                    }
                    if(sample.getDP()<12 || gq<60){
                        if(maxIndex>0){
                            lowAC++;
                        }else{
                            lowRef++;
                        }
                    }else {
                        if(maxIndex>0){
                            highAC++;
                        }else{
                            highRef++;
                        }
                    }

                }
            }
        }
        int highAcNum=0;
        int highRefNum=0;
        double acTotalThreshold=sampleSize*2*0.01;
        double acLowThreashold=1000;
        if(highAC+lowAC<acTotalThreshold){
            highAcNum=sampleSize*2+1;
            highRefNum=highRef/2;
        }else{
            if(lowAC<acLowThreashold){
                highAcNum=highAC/10>sampleSize*0.001?highAC/10:highAC;
                highRefNum=highRef/2;
            }
        }
        int addedHighAc=0;
        int addedHighNonAc=0;
        int addedLowAc=0;
        int addedLowNonAc=0;
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();
                double maxGT=-1000;
                int maxIndex=-1;
                for(int i=0;i<gls.length;i++){
                    if(gls[i]>maxGT){
                        maxGT=gls[i];
                        maxIndex=i;
                    }
                }
                if ( MathUtils.sum(gls) < GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL || keepUninformative ){
//                    if(maxIndex>0){
//                        haveAC++;
//                    }else{
//                        nonAC++;
//                    }
                    if(sample.getDP()<12 || sample.getGQ()<60) {
                        boolean exclude=false;
                        if(sample.getGQ()==-1){
                            if(sample.getDP()>15){
                                int maxNumber=0;
                                double secondMax=-1000;
                                for(int i=0;i<gls.length;i++){
                                    if(gls[i]==maxGT){
                                        maxNumber++;
                                    }else if(gls[i]>secondMax){
                                        secondMax=gls[i];
                                    }
                                }
                                if(maxGT-secondMax>=6.0){
                                    exclude=true;
                                }
                            }

                        }
                        if(exclude) {
                            if (maxIndex>0) {
                                haveAC++;
                                if(haveAC<=highAcNum) {
                                    addedHighAc++;
                                    genotypeLikelihoods.add(gls);
                                }
                            }else{
                                nonAC++;
                                if(nonAC<=highRefNum){
                                    addedHighNonAc++;
                                    genotypeLikelihoods.add(gls);
                                }
                            }
                        }else{
                            if(maxIndex>0){
                                addedLowAc++;
                            }else{
                                addedLowNonAc++;
                            }
                            genotypeLikelihoods.add(gls);
                        }
                    }else{
                        if(maxIndex>0){
                            haveAC++;
                            if(haveAC<=highAcNum) {
                                addedHighAc++;
                                genotypeLikelihoods.add(gls);
                                continue;
                            }
                        }else{
                            nonAC++;
                            if(nonAC<=highRefNum){
                                addedHighNonAc++;
                                genotypeLikelihoods.add(gls);
                                continue;
                            }
                        }
                        excludedSamples++;
                        //excludedPLs.add(gls);
                        double minPL=-1000;
                        int minIndex=-1;
                        for(int i=0;i<gls.length;i++){
                            if(gls[i]>minPL){
                                minPL=gls[i];
                                minIndex=i;
                            }
                        }
                        if(minIndex>=0) {
                            if(excludedACs.containsKey(minIndex)){
                                int v=excludedACs.get(minIndex);
                                v++;
                                excludedACs.put(minIndex,v);
                            }else{
                                excludedACs.put(minIndex, 1);
                            }
                        }
                    }
                }

            }
        }
        System.out.println("HighAC\t"+addedHighAc+"\tLowAC\t"+addedLowAc+"\tHighNonAC\t"+addedHighNonAc+"\tLowNonAC\t"+addedLowNonAc);
        return genotypeLikelihoods;
    }
    protected static ArrayList<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy, final boolean keepUninformative,boolean output) {
        final ArrayList<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);

        if ( includeDummy ) genotypeLikelihoods.add(new double[]{0.0,0.0,0.0}); // dummy
        for ( Genotype sample : GLs.iterateInSampleNameOrder() ) {
            if ( sample.hasLikelihoods() ) {
                final double[] gls = sample.getLikelihoods().getAsVector();
                
                if(output){
                	System.err.print("genotypeLikelihoods gt:"+sample+"\t");
                	for(int i = 0 ; i < gls.length ; i++)
                		System.err.print(gls[i]+"\t");
                	System.err.println();
                }
                
                if ( MathUtils.sum(gls) < GaeaGvcfVariantContextUtils.SUM_GL_THRESH_NOCALL || keepUninformative )
                    genotypeLikelihoods.add(gls);
            }
        }

        return genotypeLikelihoods;
    }

    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        // don't try to genotype too many alternate alleles
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        // only if output allele has reduced from the input alt allele set size we should care.
        final int altAlleleReduction = inputAltAlleles.size() - outputAltAlleles.size();

        if (altAlleleReduction == 0)
            return vc;

        if ( !printedMaxAltAllelesWarning ) {
            printedMaxAltAllelesWarning = true;
        }

        final List<Allele> alleles = new ArrayList<>(maximumAlternativeAlleles + 1);
        alleles.add(vc.getReference());
        alleles.addAll(reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles));
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.alleles(alleles);
        builder.genotypes(reduceScopeGenotypes(vc, defaultPloidy, alleles));
        if (altAlleleReduction < 0)
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + - altAlleleReduction + " " + vc + " " + builder.make());
        return builder.make();
    }

    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {

        // Look  for the <NON_REF> allele to exclude it from the pruning if present.
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GaeaGvcfVariantContextUtils.indexOfAltAllele(vc,
                GaeaVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing, so we need to count it out when
        // considering if alt. allele downsizing is required.
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles)
            return vc.getAlternateAlleles();

        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
            final Allele allele = vc.getAlternateAllele(i);
            likelihoodSums[i] = new LikelihoodSum(allele,i);
        }

        // Calculate the allele likelihood sums.
        reduceScopeCalculateLikelihoodSums(vc, defaultPloidy, likelihoodSums);

        // sort them by probability mass and choose the best ones
        // Make sure that the <NON_REF> allele is last if present.
        Collections.sort(Arrays.asList(likelihoodSums), nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // We need to return the best likelihood alleles in the original alternative allele index order.
        // This heap will keep track of that index order.
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);

        for ( int i = 0; i < numAllelesToChoose; i++ )
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[i]);

        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent)
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[nonRefAltAlleleIndex]);

        final ArrayList<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty())
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);

        return orderedBestAlleles;
    }

    protected static final int PL_INDEX_OF_HOM_REF = 0;

    /**
     * Update the likelihood sums with using the variant context genotype likelihoods.
     * @param vc source variant context.
     * @param likelihoodSums where to update the likelihood sums.
     */
    protected abstract void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums);

    /**
     * Transforms the genotypes of the variant context according to the new subset of possible alleles.
     *
     * @param vc original variant-context.
     * @param allelesToUse possible alleles.
     * @return never {@code null}, the new set of genotype calls for the reduced scope.
     */
    protected abstract GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse);
}
