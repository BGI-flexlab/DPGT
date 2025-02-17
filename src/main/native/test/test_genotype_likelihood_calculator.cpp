/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
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

