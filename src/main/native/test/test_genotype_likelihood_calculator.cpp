#include "genotyper/genotype_likelihood_calculators.hpp"
#include "genotyper/genotype_likelihood_calculator.hpp"



int main() {
    GenotypeLikelihoodCalculators calculators;
    GenotypeLikelihoodCalculator calculator = 
        createGenotypeLikelihoodCalculator(2, 3, calculators);
    
    std::cout << "ploidy: " << calculator.ploidy() << ", "
        << "alleles: " << calculator.alleleCount() << ", "
        << "genotypes: " << calculator.genotypeCount() << std::endl;
    
    for (int i = 0; i < 6; ++i) {
        std::cout << i << ":" << calculator.genotypeAlleleCountsAt(i)
            << std::endl;
    }
    
    return 0;
}

