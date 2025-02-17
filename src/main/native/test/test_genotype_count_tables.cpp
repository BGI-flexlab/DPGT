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
#include "boost/math/special_functions/binomial.hpp"
#include "common/simple_matrix.hpp"
#include "genotyper/genotype_likelihood_calculators.hpp"


int calcNumLikelihoods(int ploidy, int num_alleles)
{
    return static_cast<int>(boost::math::binomial_coefficient<double>(
            ploidy + num_alleles - 1, ploidy));
}

void printGenotypeCountTable(int max_ploidy, int max_allele) {
    const int rows = max_ploidy + 1;
    const int cols = max_allele + 1;
    SimpleMatrix<int> result(rows, cols);
    result.fill(0);
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            result(i, j) = calcNumLikelihoods(i, j);
        }
    }
    std::cout << result;
}

void printAlleleFirstGenotypeOffsetTable(int max_ploidy, int max_allele) {
    std::cout << GenotypeLikelihoodCalculators::buildAlleleFirstGenotypeOffsetTable(max_ploidy, max_allele);
}


int main() {
    printGenotypeCountTable(10, 10);
    printAlleleFirstGenotypeOffsetTable(10, 10);
}



