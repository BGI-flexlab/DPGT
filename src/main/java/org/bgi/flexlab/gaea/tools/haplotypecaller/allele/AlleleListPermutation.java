package org.bgi.flexlab.gaea.tools.haplotypecaller.allele;

import htsjdk.variant.variantcontext.Allele;

public interface AlleleListPermutation<A extends Allele> extends Permutation<A>, AlleleList<A> {
}
