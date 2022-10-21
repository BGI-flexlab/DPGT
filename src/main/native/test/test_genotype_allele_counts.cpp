#include <iostream>
#include "genotyper/genotype_allele_counts.hpp"


int main() {

    GenotypeAlleleCounts gac = GenotypeAlleleCounts::first(2);

    for (int i = 0; i < 6; ++i) {
        std::cout << gac << std::endl;
        gac.increase();
    }

    GenotypeAlleleCounts gac1 = GenotypeAlleleCounts();

    gac1 = gac;

    std::cout << gac1 << std::endl;

    return 0;
}

