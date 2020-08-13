/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.vcfstats;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public enum VariantType {
    SNP// Single nucleotide polymorphism (i.e. 1 base is changed)
    , MNP // Multiple nucleotide polymorphism (i.e. several bases are changed)
    , INS // Insertion (i.e. some bases added)
    , DEL // Deletion (some bases removed)
    , InDel
    , MIXED // A mixture of insertion, deletions, SNPs and or MNPs (a.k.a. substitution)
    , NO_VARIATION
    , SYMBOLIC
    , INV // Inversion (structural variant)
    , DUP // Duplication (structural variant)
    , BND; // Break-ends (rearrangement)


    public static VariantType determineType(VariantContext vc, String sample) {
        VariantType type;

        switch ( vc.getNAlleles() ) {
            case 0:
                throw new IllegalStateException("Unexpected error: requested type of VariantContext with no alleles!");
            case 1:
                // note that this doesn't require a reference allele.  You can be monomorphic independent of having a
                // reference allele
                type = NO_VARIATION;
                break;
            default:
                type = determinePolymorphicType(vc, sample);
        }
        return type;
    }

    private static VariantType determinePolymorphicType(VariantContext vc, String sample) {
        Genotype gt = vc.getGenotype(sample);
        VariantType type = null;

        // do a pairwise comparison of all alleles against the reference allele
        for ( Allele allele : gt.getAlleles() ) {
            if ( allele.isReference())
                continue;

            // find the type of this allele relative to the reference
            VariantType biallelicType = typeOfBiallelicVariant(vc.getReference(), allele);

            // for the first alternate allele, set the type to be that one
            if ( type == null ) {
                type = biallelicType;
            }
            // if the type of this allele is different from that of a previous one, assign it the MIXED type and quit
            else if ( biallelicType != type ) {
                type = MIXED;
            }
        }
        return type;
    }

    public static VariantType typeOfBiallelicVariant(Allele ref, Allele allele) {
        if ( ref.isSymbolic() )
            throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");

        if ( allele.isSymbolic() )
            return SYMBOLIC;

        if ( ref.length() == allele.length() ) {
            if ( isMNP(ref, allele))
                return MNP;
            else
                return SNP;
        } else if ( ref.length() > allele.length() ) {
            if(isInDel(ref,allele))
//            if(isInDel2(ref,allele))
                return InDel;
            else
                return DEL;
        } else
            if(isInDel(allele, ref))
//            if(isInDel2(allele, ref))
                return InDel;
            else
                return INS;
    }

    public static boolean isMNP(Allele ref, Allele alt){
        int diffBaseNum = 0;
        for (int i = 0; i < ref.length(); i++) {
            if(ref.getBaseString().charAt(i) != alt.getBaseString().charAt(i))
                diffBaseNum++;
        }
        return diffBaseNum > 1;
    }

    public static boolean isInDel(Allele longAllele, Allele shortAllele){
        return !longAllele.getBaseString().startsWith(shortAllele.getBaseString());
    }

    public static boolean isInDel2(Allele longAllele, Allele shortAllele){
        return shortAllele.length() > 1;
    }
}
