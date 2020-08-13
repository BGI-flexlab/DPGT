package org.bgi.flexlab.gaea.tools.haplotypecaller;

import java.util.ArrayList;
import java.util.List;

import org.bgi.flexlab.gaea.tools.haplotypecaller.allele.AlleleList;
import org.bgi.flexlab.gaea.tools.haplotypecaller.allele.AlleleListPermutation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculator;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood.GenotypeLikelihoodCalculators;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;

public final class IndependentSampleGenotypesModel {
    private static final int DEFAULT_CACHE_PLOIDY_CAPACITY = 10;
    private static final int DEFAULT_CACHE_ALLELE_CAPACITY = 50;

    private final int cacheAlleleCountCapacity;
    private final int cachePloidyCapacity;
    private GenotypeLikelihoodCalculator[][] likelihoodCalculators;
    private final GenotypeLikelihoodCalculators calculators;

    public IndependentSampleGenotypesModel() { this(DEFAULT_CACHE_PLOIDY_CAPACITY, DEFAULT_CACHE_ALLELE_CAPACITY); }

    /**
     *  Initialize model with given maximum allele count and ploidy for caching
     */
    public IndependentSampleGenotypesModel(final int calculatorCachePloidyCapacity, final int calculatorCacheAlleleCapacity) {
        cachePloidyCapacity = calculatorCachePloidyCapacity;
        cacheAlleleCountCapacity = calculatorCacheAlleleCapacity;
        likelihoodCalculators = new GenotypeLikelihoodCalculator[calculatorCachePloidyCapacity][calculatorCacheAlleleCapacity];
        calculators = new GenotypeLikelihoodCalculators();
    }

    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles, final GenotypingData<A> data) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = AlleleLikelihoodMatrixMapper.newInstance(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);
        final int alleleCount = genotypingAlleles.numberOfAlleles();

        GenotypeLikelihoodCalculator likelihoodsCalculator = sampleCount > 0 ? getLikelihoodsCalculator(ploidyModel.samplePloidy(0), alleleCount) : null;
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);

            // get a new likelihoodsCalculator if this sample's ploidy differs from the previous sample's
            if (samplePloidy != likelihoodsCalculator.ploidy()) {
                likelihoodsCalculator = getLikelihoodsCalculator(samplePloidy, alleleCount);
            }

            final LikelihoodMatrix<A> sampleLikelihoods = alleleLikelihoodMatrixMapper.apply(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(likelihoodsCalculator.genotypeLikelihoods(sampleLikelihoods));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }

    private GenotypeLikelihoodCalculator getLikelihoodsCalculator(final int samplePloidy, final int alleleCount) {
        if (samplePloidy >= cachePloidyCapacity || alleleCount >= cacheAlleleCountCapacity) {
            return calculators.getInstance(samplePloidy, alleleCount);
        }
        final GenotypeLikelihoodCalculator result = likelihoodCalculators[samplePloidy][alleleCount];
        if (result != null) {
            return result;
        } else {
            final GenotypeLikelihoodCalculator newOne = calculators.getInstance(samplePloidy, alleleCount);
            likelihoodCalculators[samplePloidy][alleleCount] = newOne;
            return newOne;
        }
    }
}
