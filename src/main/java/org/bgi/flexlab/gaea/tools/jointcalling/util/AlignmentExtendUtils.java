package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.Haplotype;
import org.bgi.flexlab.gaea.util.AlignmentUtil;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class AlignmentExtendUtils extends AlignmentUtil {

	public final static String HAPLOTYPE_TAG = "HC";
	
	private static class CigarPairTransform {
	    private final EnumSet<CigarOperator> op12, op23;
	    private final CigarOperator op13;
	    private final int advance12, advance23;

	    private CigarPairTransform(CigarOperator op12, CigarOperator op23, CigarOperator op13, int advance12, int advance23) {
	        this.op12 = getCigarSet(op12);
	        this.op23 = getCigarSet(op23);
	        this.op13 = op13;
	        this.advance12 = advance12;
	        this.advance23 = advance23;
	    }

	    private static EnumSet<CigarOperator> getCigarSet(final CigarOperator masterOp) {
	        switch ( masterOp ) {
	            case M: return EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X);
	            case I: return EnumSet.of(CigarOperator.I, CigarOperator.S);
	            case D: return EnumSet.of(CigarOperator.D);
	            default: throw new IllegalStateException("Unexpected state " + masterOp);
	        }
	    }

	    @Override
	    public String toString() {
	        return "CigarPairTransform{" +
	                "op12=" + op12 +
	                ", op23=" + op23 +
	                ", op13=" + op13 +
	                ", advance12=" + advance12 +
	                ", advance23=" + advance23 +
	                '}';
	    }
	}
	
	private final static List<CigarPairTransform> cigarPairTransformers = Arrays.asList(
			//
			// op12 is a match
			//
			// 3: xxx B yyy
			// ^^^^^^^^^^^^
			// 2: xxx M yyy
			// 1: xxx M yyy
			new CigarPairTransform(CigarOperator.M, CigarOperator.M, CigarOperator.M, 1, 1),
			// 3: xxx I yyy
			// ^^^^^^^^^^^^
			// 2: xxx I yyy
			// 1: xxx M yyy
			new CigarPairTransform(CigarOperator.M, CigarOperator.I, CigarOperator.I, 1, 1),
			// 3: xxx D yyy
			// ^^^^^^^^^^^^
			// 2: xxx D yyy
			// 1: xxx M yyy
			new CigarPairTransform(CigarOperator.M, CigarOperator.D, CigarOperator.D, 0, 1),

			//
			// op12 is a deletion
			//
			// 3: xxx D M yyy
			// ^^^^^^^^^^^^
			// 2: xxx M yyy
			// 1: xxx D yyy
			new CigarPairTransform(CigarOperator.D, CigarOperator.M, CigarOperator.D, 1, 1),
			// 3: xxx D1 D2 yyy
			// ^^^^^^^^^^^^
			// 2: xxx D2 yyy
			// 1: xxx D1 yyy
			new CigarPairTransform(CigarOperator.D, CigarOperator.D, CigarOperator.D, 1, 0),
			// 3: xxx X yyy => no-op, we skip emitting anything here
			// ^^^^^^^^^^^^
			// 2: xxx I yyy
			// 1: xxx D yyy
			new CigarPairTransform(CigarOperator.D, CigarOperator.I, null, 1, 1),

			//
			// op12 is a insertion
			//
			// 3: xxx I M yyy
			// ^^^^^^^^^^^^
			// 2: xxx M yyy
			// 1: xxx I yyy
			new CigarPairTransform(CigarOperator.I, CigarOperator.M, CigarOperator.I, 1, 0),
			// 3: xxx I D yyy
			// ^^^^^^^^^^^^
			// 2: xxx D yyy
			// 1: xxx I yyy
			new CigarPairTransform(CigarOperator.I, CigarOperator.D, CigarOperator.I, 1, 0),
			// 3: xxx I1 I2 yyy
			// ^^^^^^^^^^^^
			// 2: xxx I2 yyy
			// 1: xxx I1 yyy
			new CigarPairTransform(CigarOperator.I, CigarOperator.I, CigarOperator.I, 1, 0)
			);

	public static GaeaSamRecord createReadAlignedToRef(final GaeaSamRecord originalRead, final Haplotype haplotype,
			final Haplotype refHaplotype, final int referenceStart, final boolean isInformative) {
		if (originalRead == null)
			throw new IllegalArgumentException("originalRead cannot be null");
		if (haplotype == null)
			throw new IllegalArgumentException("haplotype cannot be null");
		if (refHaplotype == null)
			throw new IllegalArgumentException("ref haplotype cannot be null");
		if (haplotype.getCigar() == null)
			throw new IllegalArgumentException("Haplotype cigar not set " + haplotype);
		if (referenceStart < 1)
			throw new IllegalArgumentException("reference start much be >= 1 but got " + referenceStart);

		// compute the smith-waterman alignment of read -> haplotype
		final SWPairwiseAlignment swPairwiseAlignment = new SWPairwiseAlignment(haplotype.getBases(),
				originalRead.getReadBases());
		if (swPairwiseAlignment.getAlignmentStart2wrt1() == -1)
			// sw can fail (reasons not clear) so if it happens just don't
			// realign the read
			return originalRead;
		final Cigar swCigar = consolidateCigar(swPairwiseAlignment.getCigar());

		// since we're modifying the read we need to clone it
		final GaeaSamRecord read = new GaeaSamRecord(originalRead.getHeader(),originalRead);

		// only informative reads are given the haplotype tag to enhance
		// visualization
		if (isInformative)
			read.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());

		// compute here the read starts w.r.t. the reference from the SW result
		// and the hap -> ref cigar
		final Cigar extendedHaplotypeCigar = getConsolidatedPaddedCigar(haplotype.getCigar(), 1000);
		final int readStartOnHaplotype = calcFirstBaseMatchingReferenceInCigar(extendedHaplotypeCigar,
				swPairwiseAlignment.getAlignmentStart2wrt1());
		final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnHaplotype;
		read.setAlignmentStart(readStartOnReference);
		read.resetSoftStartAndEnd();

		// compute the read -> ref alignment by mapping read -> hap -> ref from
		// the
		// SW of read -> hap mapped through the given by hap -> ref
		final Cigar haplotypeToRef = trimCigarByBases(extendedHaplotypeCigar,
				swPairwiseAlignment.getAlignmentStart2wrt1(), extendedHaplotypeCigar.getReadLength() - 1);
		final Cigar readToRefCigarRaw = applyCigarToCigar(swCigar, haplotypeToRef);
		final Cigar readToRefCigarClean = cleanUpCigar(readToRefCigarRaw);
		final Cigar readToRefCigar = leftAlignIndel(readToRefCigarClean, refHaplotype.getBases(),
				originalRead.getReadBases(), swPairwiseAlignment.getAlignmentStart2wrt1(), 0, true);

		read.setCigar(readToRefCigar);

		if (readToRefCigar.getReadLength() != read.getReadLength())
			throw new IllegalStateException("Cigar " + readToRefCigar + " with read length "
					+ readToRefCigar.getReadLength() + " != read length " + read.getReadLength() + " for read "
					+ read.format() + "\nhapToRef " + haplotypeToRef + " length " + haplotypeToRef.getReadLength() + "/"
					+ haplotypeToRef.getReferenceLength() + "\nreadToHap " + swCigar + " length "
					+ swCigar.getReadLength() + "/" + swCigar.getReferenceLength());

		return read;
	}

	public static Cigar getConsolidatedPaddedCigar(Cigar cigar, final int padSize) {
		if (padSize < 0)
			throw new IllegalArgumentException("padSize must be >= 0 but got " + padSize);
		final Cigar extendedHaplotypeCigar = new Cigar(cigar.getCigarElements());
		if (padSize > 0)
			extendedHaplotypeCigar.add(new CigarElement(padSize, CigarOperator.M));
		return consolidateCigar(extendedHaplotypeCigar);
	}

	public static Cigar consolidateCigar(final Cigar c) {
		if (c == null) {
			throw new IllegalArgumentException("Cigar cannot be null");
		}

		// fast check to determine if there's anything worth doing before we
		// create new Cigar and actually do some work
		if (!needsConsolidation(c))
			return c;

		final Cigar returnCigar = new Cigar();
		int sumLength = 0;
		CigarElement lastElement = null;

		for (final CigarElement cur : c.getCigarElements()) {
			if (cur.getLength() == 0)
				continue; // don't add elements of 0 length

			if (lastElement != null && lastElement.getOperator() != cur.getOperator()) {
				returnCigar.add(new CigarElement(sumLength, lastElement.getOperator()));
				sumLength = 0;
			}

			sumLength += cur.getLength();
			lastElement = cur;
		}

		if (sumLength > 0) {
			returnCigar.add(new CigarElement(sumLength, lastElement.getOperator()));
		}

		return returnCigar;
	}

	private static boolean needsConsolidation(final Cigar c) {
		if (c.numCigarElements() <= 1)
			return false; // fast path for empty or single cigar

		CigarOperator lastOp = null;
		for (final CigarElement cur : c.getCigarElements()) {
			if (cur.getLength() == 0 || lastOp == cur.getOperator())
				return true;
			lastOp = cur.getOperator();
		}

		return false;
	}

	public static int calcFirstBaseMatchingReferenceInCigar(final Cigar cigar, int readStartByBaseOfCigar) {
		if (cigar == null)
			throw new IllegalArgumentException("cigar cannot be null");
		if (readStartByBaseOfCigar >= cigar.getReadLength())
			throw new IllegalArgumentException("readStartByBaseOfCigar " + readStartByBaseOfCigar
					+ " must be <= readLength " + cigar.getReadLength());

		int hapOffset = 0, refOffset = 0;
		for (final CigarElement ce : cigar.getCigarElements()) {
			for (int i = 0; i < ce.getLength(); i++) {
				switch (ce.getOperator()) {
				case M:
				case EQ:
				case X:
					if (hapOffset >= readStartByBaseOfCigar)
						return refOffset;
					hapOffset++;
					refOffset++;
					break;
				case I:
				case S:
					hapOffset++;
					break;
				case D:
					refOffset++;
					break;
				default:
					throw new IllegalStateException("calcFirstBaseMatchingReferenceInCigar does not support cigar "
							+ ce.getOperator() + " in cigar " + cigar);
				}
			}
		}

		throw new IllegalStateException("Never found appropriate matching state for cigar " + cigar + " given start of "
				+ readStartByBaseOfCigar);
	}

	public static Cigar trimCigarByBases(final Cigar cigar, final int start, final int end) {
		if (start < 0)
			throw new IllegalArgumentException("Start must be >= 0 but got " + start);
		if (end < start)
			throw new IllegalArgumentException("End " + end + " is < start = " + start);
		if (end > cigar.getReadLength())
			throw new IllegalArgumentException("End is beyond the cigar's read length " + end + " for cigar " + cigar);

		final Cigar result = trimCigar(cigar, start, end, false);

		final int expectedSize = end - start + 1;
		if (result.getReadLength() != expectedSize)
			throw new IllegalStateException("trimCigarByBases failure: start " + start + " end " + end + " for " + cigar
					+ " resulted in cigar with wrong size " + result + " with size " + result.getReadLength()
					+ " expected " + expectedSize + " for input cigar " + cigar);
		return result;
	}

	private static Cigar trimCigar(final Cigar cigar, final int start, final int end, final boolean byReference) {
		final List<CigarElement> newElements = new LinkedList<CigarElement>();

		int pos = 0;
		for (final CigarElement elt : cigar.getCigarElements()) {
			if (pos > end && (byReference || elt.getOperator() != CigarOperator.D))
				break;

			switch (elt.getOperator()) {
			case D:
				if (!byReference) {
					if (pos >= start)
						newElements.add(elt);
					break;
				}
				// otherwise fall through to the next case
			case EQ:
			case M:
			case X:
				pos = addCigarElements(newElements, pos, start, end, elt);
				break;
			case S:
			case I:
				if (byReference) {
					if (pos >= start)
						newElements.add(elt);
				} else {
					pos = addCigarElements(newElements, pos, start, end, elt);
				}
				break;
			default:
				throw new IllegalStateException("Cannot handle " + elt);
			}
		}

		return consolidateCigar(new Cigar(newElements));
	}

	protected static int addCigarElements(final List<CigarElement> dest, int pos, final int start, final int end,
			final CigarElement elt) {
		final int length = Math.min(pos + elt.getLength() - 1, end) - Math.max(pos, start) + 1;
		if (length > 0)
			dest.add(new CigarElement(length, elt.getOperator()));
		return pos + elt.getLength();
	}

	public static Cigar applyCigarToCigar(final Cigar firstToSecond, final Cigar secondToThird) {
        final boolean DEBUG = false;

        final List<CigarElement> newElements = new LinkedList<CigarElement>();
        final int nElements12 = firstToSecond.getCigarElements().size();
        final int nElements23 = secondToThird.getCigarElements().size();

        int cigar12I = 0, cigar23I = 0;
        int elt12I = 0, elt23I = 0;

        while ( cigar12I < nElements12 && cigar23I < nElements23 ) {
            final CigarElement elt12 = firstToSecond.getCigarElement(cigar12I);
            final CigarElement elt23 = secondToThird.getCigarElement(cigar23I);

            final CigarPairTransform transform = getTransformer(elt12.getOperator(), elt23.getOperator());

            if ( DEBUG )
                System.out.printf("Transform %s => %s with elt1 = %d %s @ %d elt2 = %d %s @ %d with transform %s%n",
                        firstToSecond, secondToThird, cigar12I, elt12.getOperator(), elt12I, cigar23I, elt23.getOperator(), elt23I, transform);

            if ( transform.op13 != null ) // skip no ops
                newElements.add(new CigarElement(1, transform.op13));

            elt12I += transform.advance12;
            elt23I += transform.advance23;

            // if have exhausted our current element, advance to the next one
            if ( elt12I == elt12.getLength() ) { cigar12I++; elt12I = 0; }
            if ( elt23I == elt23.getLength() ) { cigar23I++; elt23I = 0; }
        }

        return consolidateCigar(new Cigar(newElements));
    }
	
	private static CigarPairTransform getTransformer(final CigarOperator op12, final CigarOperator op23) {
        for ( final CigarPairTransform transform : cigarPairTransformers) {
            if ( transform.op12.contains(op12) && transform.op23.contains(op23) )
                return transform;
        }

        throw new IllegalStateException("No transformer for operators " + op12 + " and " + op23);
    }
	
	public static Cigar cleanUpCigar(final Cigar c) {
        final List<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements() - 1);

        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() != 0 && (! elements.isEmpty() || ce.getOperator() != CigarOperator.D)) {
                elements.add(ce);
            }
        }

        return new Cigar(elements);
    }
	
    public static Cigar leftAlignIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final boolean doNotThrowExceptionForMultipleIndels) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        final int numIndels = countIndelElements(cigar);
        if ( numIndels == 0 )
            return cigar;
        if ( numIndels == 1 )
            return leftAlignSingleIndel(cigar, refSeq, readSeq, refIndex, readIndex, true);

        // if we got here then there is more than 1 indel in the alignment
        if ( doNotThrowExceptionForMultipleIndels )
            return cigar;

        throw new UnsupportedOperationException("attempting to left align a CIGAR that has more than 1 indel in its alignment but this functionality has not been implemented yet");
    }
	
	private static void ensureLeftAlignmentHasGoodArguments(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {
        if ( cigar == null ) throw new IllegalArgumentException("attempting to left align a CIGAR that is null");
        if ( refSeq == null ) throw new IllegalArgumentException("attempting to left align a reference sequence that is null");
        if ( readSeq == null ) throw new IllegalArgumentException("attempting to left align a read sequence that is null");
        if ( refIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a reference index less than 0");
        if ( readIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a read index less than 0");
    }

    /**
     * Counts the number of I/D operators
     *
     * @param cigar   cigar to check -- cannot be null
     * @return  non-negative count of indel operators
     */
    private static int countIndelElements(final Cigar cigar) {
        int indelCount = 0;
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I )
                indelCount++;
        }
        return indelCount;
    }
    
    public static Cigar leftAlignSingleIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final boolean cleanupCigar) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        int indexOfIndel = -1;
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                // if there is more than 1 indel, exception out
                if (indexOfIndel != -1)
                    throw new IllegalArgumentException("attempting to left align a CIGAR that has more than 1 indel in its alignment");
                indexOfIndel = i;
            }
        }

        // if there is no indel, exception out
        if ( indexOfIndel == -1 )
            throw new IllegalArgumentException("attempting to left align a CIGAR that has no indels in its alignment");
        // if the alignment starts with an insertion (so that there is no place on the read to move that insertion further left), we are done
        if ( indexOfIndel == 0 )
            return cigar;

        final int indelLength = cigar.getCigarElement(indexOfIndel).getLength();

        byte[] altString = createIndelString(cigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);
        if (altString == null)
            return cigar;

        Cigar newCigar = cigar;
        for (int i = 0; i < indelLength; i++) {
            newCigar = moveCigarLeft(newCigar, indexOfIndel);
            byte[] newAltString = createIndelString(newCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);

            // check to make sure we haven't run off the end of the read
            boolean reachedEndOfRead = cigarHasZeroSizeElement(newCigar);

            if (Arrays.equals(altString, newAltString)) {
                cigar = newCigar;
                i = -1;
                if (reachedEndOfRead)
                    cigar = cleanupCigar ? cleanUpCigar(cigar) : cigar;
            }

            if (reachedEndOfRead)
                break;
        }

        return cigar;
    }
    
    protected static boolean cigarHasZeroSizeElement(final Cigar c) {
        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() == 0)
                return true;
        }
        return false;
    }
    
    private static Cigar moveCigarLeft(Cigar cigar, int indexOfIndel) {
        // get the first few elements
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for (int i = 0; i < indexOfIndel - 1; i++)
            elements.add(cigar.getCigarElement(i));

        // get the indel element and move it left one base
        CigarElement ce = cigar.getCigarElement(indexOfIndel - 1);
        elements.add(new CigarElement(Math.max(ce.getLength() - 1, 0), ce.getOperator()));
        elements.add(cigar.getCigarElement(indexOfIndel));
        if (indexOfIndel + 1 < cigar.numCigarElements()) {
            ce = cigar.getCigarElement(indexOfIndel + 1);
            elements.add(new CigarElement(ce.getLength() + 1, ce.getOperator()));
        } else {
            elements.add(new CigarElement(1, CigarOperator.M));
        }

        // get the last few elements
        for (int i = indexOfIndel + 2; i < cigar.numCigarElements(); i++)
            elements.add(cigar.getCigarElement(i));
        return new Cigar(elements);
    }
    
    private static byte[] createIndelString(final Cigar cigar, final int indexOfIndel, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        CigarElement indel = cigar.getCigarElement(indexOfIndel);
        int indelLength = indel.getLength();

        int totalRefBases = 0;
        for (int i = 0; i < indexOfIndel; i++) {
            CigarElement ce = cigar.getCigarElement(i);
            int length = ce.getLength();

            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    readIndex += length;
                    refIndex += length;
                    totalRefBases += length;
                    break;
                case S:
                    readIndex += length;
                    break;
                case N:
                    refIndex += length;
                    totalRefBases += length;
                    break;
                default:
                    break;
            }
        }

        // sometimes, when there are very large known indels, we won't have enough reference sequence to cover them
        if (totalRefBases + indelLength > refSeq.length)
            indelLength -= (totalRefBases + indelLength - refSeq.length);

        // the indel-based reference string
        byte[] alt = new byte[refSeq.length + (indelLength * (indel.getOperator() == CigarOperator.D ? -1 : 1))];

        // add the bases before the indel, making sure it's not aligned off the end of the reference
        if (refIndex > alt.length || refIndex > refSeq.length)
            return null;
        System.arraycopy(refSeq, 0, alt, 0, refIndex);
        int currentPos = refIndex;

        // take care of the indel
        if (indel.getOperator() == CigarOperator.D) {
            refIndex += indelLength;
        } else {
            System.arraycopy(readSeq, readIndex, alt, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel, making sure it's not aligned off the end of the reference
        if (refSeq.length - refIndex > alt.length - currentPos)
            return null;
        System.arraycopy(refSeq, refIndex, alt, currentPos, refSeq.length - refIndex);

        return alt;
    }
}
