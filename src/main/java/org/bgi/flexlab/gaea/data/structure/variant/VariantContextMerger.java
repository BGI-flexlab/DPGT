package org.bgi.flexlab.gaea.data.structure.variant;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextMerger {
    public final static String NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
    public final static Allele NON_REF_SYMBOLIC_ALLELE = Allele.create("<"+NON_REF_SYMBOLIC_ALLELE_NAME+">", false);

	static Allele determineReferenceAllele(final List<VariantContext> VCs, final GenomeLocation location) {
		Allele ref = null;

		for (VariantContext vc : VCs) {
			if (location == null || vc.getStart() == location.getStart()) {
				Allele vcRef = vc.getReference();
				if (ref == null || ref.length() < vcRef.length())
					ref = vcRef;
				else if (ref.length() == vcRef.length() && !ref.equals(vcRef))
					throw new UserException(
							String.format("Inconsistent references for the same position at %s:%d, %s vs. %s",
									vc.getContig(), vc.getStart(), ref.toString(), vcRef.toString()));
			}
		}

		return ref;
	}
	
	static List<Allele> replaceWithNoCallsAndDels(final VariantContext vc) {
	    if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");

	    final List<Allele> result = new ArrayList<>(vc.getNAlleles());

	    // no-call the reference allele
	    result.add(Allele.NO_CALL);

	    for ( final Allele allele : vc.getAlternateAlleles() ) {
	        final Allele replacement;
	        if ( allele.equals(NON_REF_SYMBOLIC_ALLELE) )
	            replacement = allele;
	        else if ( allele.length() < vc.getReference().length() )
	            replacement = Allele.SPAN_DEL;
	        else
	            replacement = Allele.NO_CALL;

	        result.add(replacement);
	    }
	    return result;
	}
	
	static List<Allele> remapAlleles(final VariantContext vc, final Allele refAllele, final LinkedHashSet<Allele> finalAlleles) {

	    final Allele vcRef = vc.getReference();
	    final byte[] refBases = refAllele.getBases();
	    final int extraBaseCount = refBases.length - vcRef.getBases().length;
	    if (extraBaseCount < 0) throw new IllegalStateException("the wrong reference was selected");

	    final List<Allele> result = new ArrayList<>(vc.getNAlleles());
	    result.add(refAllele);

	    for (final Allele a : vc.getAlternateAlleles()) {
	        if (a.isSymbolic()) {
	            result.add(a);
	            // we always skip <NON_REF> when adding to finalAlleles; this is done outside if it applies.
	            // we also skip <*:DEL> if there isn't a real alternate allele.
	            if ( !a.equals(NON_REF_SYMBOLIC_ALLELE) && !vc.isSymbolic() )
	                finalAlleles.add(a);
	        } else if ( a == Allele.SPAN_DEL ) {
	            result.add(a);
	            // we skip * if there isn't a real alternate allele.
	            if ( !vc.isBiallelic() )
	                finalAlleles.add(a);
	        } else if (a.isCalled()) {
	            final Allele newAllele;
	            if (extraBaseCount > 0) {
	                final byte[] oldBases = a.getBases();
	                final byte[] newBases = Arrays.copyOf(oldBases,oldBases.length + extraBaseCount);
	                System.arraycopy(refBases,refBases.length - extraBaseCount,newBases,oldBases.length,extraBaseCount);
	                newAllele = Allele.create(newBases,false);
	            } else
	                newAllele = a;
	            result.add(newAllele);
	            finalAlleles.add(newAllele);
	        } else { // NO_CALL and strange miscellanea
	            result.add(a);
	        }
	    }
	    return result;
	}
}
