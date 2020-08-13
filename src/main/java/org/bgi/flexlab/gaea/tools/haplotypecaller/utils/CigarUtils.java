package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SWParameters;
import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SmithWatermanAligner;
import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SmithWatermanAlignment;
import org.bgi.flexlab.gaea.tools.haplotypecaller.smithwaterman.SmithWatermanJavaAligner.SWOverhangStrategy;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class CigarUtils {
	private static final String SW_PAD = "NNNNNNNNNN";
	
	public static int countRefBasesBasedOnCigar(final GaeaSamRecord read, final int cigarStartIndex, final int cigarEndIndex){
        if (read == null){
            throw new IllegalArgumentException("null read");
        }
        final List<CigarElement> elems = read.getCigar().getCigarElements();
        if (cigarStartIndex < 0 || cigarEndIndex > elems.size() || cigarStartIndex > cigarEndIndex){
            throw new IllegalArgumentException("invalid index:" + 0 + " -" + elems.size());
        }
        int result = 0;
        for(int i = cigarStartIndex; i < cigarEndIndex; i++){
            final CigarElement cigarElement = elems.get(i);
            switch (cigarElement.getOperator()) {
                case M:
                case D:
                case N:
                case EQ:
                case X:
                case S:
                case H:
                    result += cigarElement.getLength();
                    break;
                case I:
                case P:        //for these two, nothing happens.
                    break;
                default:
                    throw new UserException("Unsupported cigar operator: " + cigarElement.getOperator());
            }
        }
        return result;
    }
	
	public static boolean containsNOperator(final List<CigarElement> cigarElements) {
	    return Utils.nonNull(cigarElements).stream().anyMatch(el -> el.getOperator() == CigarOperator.N);
	}
	
	public static boolean isGood(final Cigar c) {
	    Utils.nonNull(c, "cigar is null");

	    //Note: the string comes from the SAMRecord so it must be a wellformed CIGAR (that is, in "\*|([0-9]+[MIDNSHPX=])+" as per SAM spec).
	    //We don't have to check that
	    if (c.isValid(null, -1) != null){  //if it's invalid, then it's not good
	        return false;
	    }
	    final List<CigarElement> elems = c.getCigarElements();
	    if (hasConsecutiveIndels(elems)){
	        return false;
	    }
	    if (startsWithDeletionIgnoringClips(elems)){
	        return false;
	    }
	    //revert the list and check deletions at the end
	    final List<CigarElement> elemsRev = new ArrayList<>(elems);
	    Collections.reverse(elemsRev);
	    return !startsWithDeletionIgnoringClips(elemsRev);
	}
	
	private static boolean hasConsecutiveIndels(final List<CigarElement> elems) {
	    boolean prevIndel = false;
	    for (final CigarElement elem : elems) {
	        final CigarOperator op = elem.getOperator();
	        final boolean isIndel = (op == CigarOperator.INSERTION || op == CigarOperator.DELETION);
	        if (prevIndel && isIndel) {
	            return true;
	        }
	        prevIndel = isIndel;
	    }
	    return false;
	}
	
	private static boolean startsWithDeletionIgnoringClips(final List<CigarElement> elems) {
	    final Iterator<CigarElement> iter = elems.iterator();
	    boolean isClip = true;
	    CigarOperator op = null;
	    while(iter.hasNext() && isClip) { //consume clips at the beginning
	        final CigarElement elem = iter.next();
	        op = elem.getOperator();
	        isClip = (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP);
	    }
	    //once all clips are consumed, is it a deletion or not?
	    return op == CigarOperator.DELETION;
	}
	
	public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq, final SmithWatermanAligner aligner) {
	    Utils.nonNull(refSeq, "refSeq");
	    Utils.nonNull(altSeq, "altSeq");
	    if ( altSeq.length == 0 ) {
	        // horrible edge case from the unit tests, where this path has no bases
	        return new Cigar(Collections.singletonList(new CigarElement(refSeq.length, CigarOperator.D)));
	    }

	    //Note: this is a performance optimization.
	    // If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
	    if (Arrays.equals(refSeq, altSeq)){
	        final Cigar matching = new Cigar();
	        matching.add(new CigarElement(refSeq.length, CigarOperator.MATCH_OR_MISMATCH));
	        return matching;
	    }

	    final Cigar nonStandard;

	    final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
	    final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
	    final SmithWatermanAlignment alignment = aligner.align(paddedRef.getBytes(), paddedPath.getBytes(), new SWParameters(200, -150, -260, -11), SWOverhangStrategy.SOFTCLIP);

	    if ( isSWFailure(alignment) ) {
	        return null;
	    }


	    // cut off the padding bases
	    final int baseStart = SW_PAD.length();
	    final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
	    nonStandard = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

	    if ( nonStandard.getReferenceLength() != refSeq.length ) {
	        nonStandard.add(new CigarElement(refSeq.length - nonStandard.getReferenceLength(), CigarOperator.D));
	    }

	    // finally, return the cigar with all indels left aligned
	    return leftAlignCigarSequentially(nonStandard, refSeq, altSeq, 0, 0);
	}
	
	private static boolean isSWFailure(final SmithWatermanAlignment alignment) {
	    // check that the alignment starts at the first base, which it should given the padding
	    if ( alignment.getAlignmentOffset() > 0 ) {
	        return true;
	    }

	    // check that we aren't getting any S operators (which would be very bad downstream)
	    for ( final CigarElement ce : alignment.getCigar().getCigarElements() ) {
	        if ( ce.getOperator() == CigarOperator.S )
	            return true;
	        // soft clips at the end of the alignment are really insertions
	    }

	    return false;
	}
	
	public static Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
	    Utils.nonNull(cigar, "cigar null");
	    Utils.nonNull(refSeq, "refSeq null");
	    Utils.nonNull(readSeq, "readSeq null");

	    final Cigar cigarToReturn = new Cigar();
	    Cigar cigarToAlign = new Cigar();
	    for (int i = 0; i < cigar.numCigarElements(); i++) {
	        final CigarElement ce = cigar.getCigarElement(i);
	        if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
	            cigarToAlign.add(ce);
	            final Cigar leftAligned = AlignmentUtils.leftAlignSingleIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false);
	            for ( final CigarElement toAdd : leftAligned.getCigarElements() ) { cigarToReturn.add(toAdd); }
	            refIndex += cigarToAlign.getReferenceLength();
	            readIndex += cigarToAlign.getReadLength();
	            cigarToAlign = new Cigar();
	        } else {
	            cigarToAlign.add(ce);
	        }
	    }
	    if( !cigarToAlign.isEmpty() ) {
	        for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
	            cigarToReturn.add(toAdd);
	        }
	    }

	    final Cigar result = AlignmentUtils.consolidateCigar(cigarToReturn);
	    Utils.validateArg(result.getReferenceLength() == cigar.getReferenceLength(),
	            () -> "leftAlignCigarSequentially failed to produce a valid CIGAR.  Reference lengths differ.  Initial cigar " + cigar + " left aligned into " + result);
	    return result;
	}
	
	public static int countLeftHardClippedBases(final Cigar cigar) {
	    Utils.nonNull(cigar, "the input cigar cannot not be null");
	    if (cigar.numCigarElements() < 2) {
	        return 0;
	    } else if (cigar.getCigarElement(0).getOperator() != CigarOperator.H) {
	        return 0;
	    } else {
	        return cigar.getCigarElement(0).getLength();
	    }
	}
	
	public static int countRightHardClippedBases(final Cigar cigar) {
	    Utils.nonNull(cigar, "the input cigar cannot not be null");
	    if (cigar.numCigarElements() < 2) {
	        return 0;
	    } else {
	        final List<CigarElement> elements = cigar.getCigarElements();
	        int lastElementIndex;
	        if (elements.get(lastElementIndex = elements.size() - 1).getOperator() != CigarOperator.H) {
	            return 0;
	        } else {
	            return elements.get(lastElementIndex).getLength();
	        }
	    }
	}
}
