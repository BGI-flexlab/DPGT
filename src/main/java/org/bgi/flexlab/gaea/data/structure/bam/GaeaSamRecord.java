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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.bam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.sequenceplatform.NGSPlatform;
import org.bgi.flexlab.gaea.util.EventType;
import org.bgi.flexlab.gaea.util.ReadUtils;

import java.nio.charset.Charset;
import java.util.*;

public class GaeaSamRecord extends SAMRecord implements Locatable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 8070430647521989737L;
	
	private static final Charset DEFAULT_CHARSET = Charset.forName("UTF-8");
	
	// ReduceReads specific attribute tags
	public static final String REDUCED_READ_CONSENSUS_TAG = "RR";
	// Base Quality Score Recalibrator specific attribute tags
	public static final String BQSR_BASE_INSERTION_QUALITIES = "BI";
	// base qualities for deletions
	public static final String BQSR_BASE_DELETION_QUALITIES = "BD";
	private static final int PROPER_PAIR_FLAG = 0x2;
	private static final int MATE_UNMAPPED_FLAG = 0x8;
	private static final int MATE_STRAND_FLAG = 0x20;
	private static final int FIRST_OF_PAIR_FLAG = 0x40;
	private static final int SECOND_OF_PAIR_FLAG = 0x80;
	private static final int DUPLICATE_READ_FLAG = 0x400;
	public static final String NO_ALIGNMENT_REFERENCE_NAME = "*";
	public static final int NO_ALIGNMENT_REFERENCE_INDEX = -1;
	public static final int UNSET_POSITION = 0;

	private int softStart = -1;
	private int softEnd = -1;
	/***
	 * Realignment need attribute
	 */
	private int originalAlignmentStart = -1;

	private boolean mustBeOut = false;

	private Map<Object, Object> temporaryAttributes;

	public GaeaSamRecord(SAMFileHeader header) {
		super(header);
	}

	public GaeaSamRecord(SAMFileHeader header, SAMRecord sam) {
		super(header);
		set(sam);
	}

	public GaeaSamRecord(SAMFileHeader header, SAMRecord sam, boolean setOrignalPos) {
		this(header, sam);
		mustBeOut = setOrignalPos;
	}

	public boolean needToOutput() {
		return mustBeOut;
	}

	public GaeaSamRecord(SAMFileHeader header, int referenceSequenceIndex, int alignmentStart, short mappingQuality,
			int flags, int mateReferenceSequenceIndex, int mateAlignmentStart, int insertSize) {
		super(header);
		this.setHeader(header);
		this.setReferenceIndex(referenceSequenceIndex);
		this.setAlignmentStart(alignmentStart);
		this.setMappingQuality(mappingQuality);
		this.setFlags(flags);
		this.setMateReferenceIndex(mateReferenceSequenceIndex);
		this.setMateAlignmentStart(mateAlignmentStart);
		this.setInferredInsertSize(insertSize);
	}

	private static boolean hasReferenceName(final Integer referenceIndex, final String referenceName) {
		return (referenceIndex != null && referenceIndex != NO_ALIGNMENT_REFERENCE_INDEX)
				|| !NO_ALIGNMENT_REFERENCE_NAME.equals(referenceName);
	}

	/**
	 * @return true if this SAMRecord has a reference, either as a String or
	 *         index (or both).
	 */
	public boolean hasReferenceName() {
		return hasReferenceName(getMateReferenceIndex(), getReferenceName());
	}

	/**
	 * @return true if this SAMRecord has a mate reference, either as a String
	 *         or index (or both).
	 */
	public boolean hasMateReferenceName() {
		return hasReferenceName(getMateReferenceIndex(), getMateReferenceName());
	}

	public void requireReadPaired() {
		if (!getReadPairedFlag()) {
			throw new IllegalStateException("Inappropriate call if not paired read");
		}
	}

	public boolean getProperPairFlagUnchecked() {
		return (getFlags() & PROPER_PAIR_FLAG) != 0;
	}

	public boolean getMateUnmappedFlagUnchecked() {
		return (getFlags() & MATE_UNMAPPED_FLAG) != 0;
	}

	public boolean getMateNegativeStrandFlagUnchecked() {
		return (getFlags() & MATE_STRAND_FLAG) != 0;
	}

	public boolean getFirstOfPairFlagUnchecked() {
		return (getFlags() & FIRST_OF_PAIR_FLAG) != 0;
	}

	public boolean getSecondOfPairFlagUnchecked() {
		return (getFlags() & SECOND_OF_PAIR_FLAG) != 0;
	}

	public void requireSigned(final String tag) {
		if (isUnsignedArrayAttribute(tag))
			throw new SAMException("Value for tag " + tag + " is not signed");
	}

	public void requireUnsigned(final String tag) {
		if (!isUnsignedArrayAttribute(tag))
			throw new SAMException("Value for tag " + tag + " is not unsigned");
	}

	public boolean safeEquals(final Object o1, final Object o2) {
		if (o1 == o2) {
			return true;
		} else if (o1 == null || o2 == null) {
			return false;
		} else {
			return o1.equals(o2);
		}
	}

	public List<SAMValidationError> isValidReferenceIndexAndPosition(final Integer referenceIndex,
			final String referenceName, final int alignmentStart, final boolean isMate) {
		final boolean hasReference = hasReferenceName(referenceIndex, referenceName);

		ArrayList<SAMValidationError> ret = null;
		if (!hasReference) {
			if (alignmentStart != 0) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START,
						buildMessage("Alignment start should be 0 because reference name = *.", isMate),
						getReadName()));
			}
		} else {
			if (alignmentStart == 0) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START,
						buildMessage("Alignment start should != 0 because reference name != *.", isMate),
						getReadName()));
			}

			if (getHeader().getSequenceDictionary().size() > 0) {
				final SAMSequenceRecord sequence = (referenceIndex != null ? getHeader().getSequence(referenceIndex)
						: getHeader().getSequence(referenceName));
				if (sequence == null) {
					if (ret == null)
						ret = new ArrayList<SAMValidationError>();
					ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_REFERENCE_INDEX,
							buildMessage("Reference sequence not found in sequence dictionary.", isMate),
							getReadName()));
				} else {
					if (alignmentStart > sequence.getSequenceLength()) {
						if (ret == null)
							ret = new ArrayList<SAMValidationError>();
						ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_ALIGNMENT_START, buildMessage(
								"Alignment start (" + alignmentStart + ") must be <= reference sequence length ("
										+ sequence.getSequenceLength() + ") on reference " + sequence.getSequenceName(),
								isMate), getReadName()));
					}
				}
			}
		}
		return ret;
	}

	public String buildMessage(final String baseMessage, final boolean isMate) {
		return isMate ? "Mate " + baseMessage : baseMessage;
	}

	public boolean isDuplicateRead() {
		if ((this.getFlags() & DUPLICATE_READ_FLAG) > 0) {
			return true;
		}
		return false;
	}

	public void removeDuplicateFlag() {
		if ((this.getFlags() & DUPLICATE_READ_FLAG) > 0) {
			this.setFlags(this.getFlags() - DUPLICATE_READ_FLAG);
		}
	}

	public int getClipReadLength() {
		int len = 0;
		final List<CigarElement> cigs = this.getCigar().getCigarElements();
		for (int i = cigs.size() - 1; i >= 0; --i) {
			final CigarElement cig = cigs.get(i);
			final CigarOperator op = cig.getOperator();

			if (!(op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP || op == CigarOperator.DELETION)) {
				len += cig.getLength();
			}
		}

		return len;
	}

	public void set(SAMRecord sam) {
		this.setReadName(sam.getReadName());
		this.setFlags(sam.getFlags());
		this.setReferenceIndex(sam.getReferenceIndex());
		this.setAlignmentStart(sam.getAlignmentStart());
		this.setMappingQuality(sam.getMappingQuality());
		this.setCigar(sam.getCigar());
		this.setMateReferenceIndex(sam.getMateReferenceIndex());
		this.setMateAlignmentStart(sam.getMateAlignmentStart());
		this.setInferredInsertSize(sam.getInferredInsertSize());
		this.setReadString(sam.getReadString());
		this.setBaseQualityString(sam.getBaseQualityString());
		for (SAMTagAndValue tag : sam.getAttributes()) {
			this.setAttribute(tag.tag, tag.value);
		}
	}

	public void removeAtrributes(String tag) {
		List<SAMTagAndValue> samTag = this.getAttributes();
		this.clearAttributes();
		for (SAMTagAndValue samTagAndValue : samTag) {
			if (!samTagAndValue.tag.equals(tag)) {
				this.setAttribute(samTagAndValue.tag, samTagAndValue.value);
			}
		}
	}
	
	public GaeaSamRecord copy(){
		GaeaSamRecord newRecord = null;
		try {
			newRecord = (GaeaSamRecord) super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e.toString());
		}
		return newRecord;
	}

	public Object clone() throws CloneNotSupportedException {
		final GaeaSamRecord newRecord = (GaeaSamRecord) super.clone();
		return newRecord;
	}

	public boolean getSupplementaryAlignmentFlag() {
		return super.getSupplementaryAlignmentFlag();
	}

	public int getOriginalAlignmentStart() {
		return originalAlignmentStart;
	}

	public void setOriginalAlignmentStart(int originalAlignmentStart) {
		this.originalAlignmentStart = originalAlignmentStart;
	}

	public List<SAMValidationError> isValid() {
		// ret is only instantiate if there are errors to report, in order to
		// reduce GC in the typical case
		// in which everything is valid. It's ugly, but more efficient.
		ArrayList<SAMValidationError> ret = null;
		if (!getReadPairedFlag()) {
			if (getProperPairFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_PROPER_PAIR,
						"Proper pair flag should not be set for unpaired read.", getReadName()));
			}
			if (getMateUnmappedFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED,
						"Mate unmapped flag should not be set for unpaired read.", getReadName()));
			}
			if (getMateNegativeStrandFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_NEG_STRAND,
						"Mate negative strand flag should not be set for unpaired read.", getReadName()));
			}
			if (getFirstOfPairFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_FIRST_OF_PAIR,
						"First of pair flag should not be set for unpaired read.", getReadName()));
			}
			if (getSecondOfPairFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_SECOND_OF_PAIR,
						"Second of pair flag should not be set for unpaired read.", getReadName()));
			}
			if (getMateReferenceIndex() != NO_ALIGNMENT_REFERENCE_INDEX) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_MATE_REF_INDEX,
						"MRNM should not be set for unpaired read.", getReadName()));
			}
		} else {
			final List<SAMValidationError> errors = isValidReferenceIndexAndPosition(mMateReferenceIndex,
					getMateReferenceName(), getMateAlignmentStart(), true);
			if (errors != null) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.addAll(errors);
			}
			if (!hasMateReferenceName() && !getMateUnmappedFlag()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_MATE_UNMAPPED,
						"Mapped mate should have mate reference name", getReadName()));
			}
			if (!getFirstOfPairFlagUnchecked() && !getSecondOfPairFlagUnchecked()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND,
						"Paired read should be marked as first of pair or second of pair.", getReadName()));
			}
		}
		if (getInferredInsertSize() > MAX_INSERT_SIZE || getInferredInsertSize() < -MAX_INSERT_SIZE) {
			if (ret == null)
				ret = new ArrayList<SAMValidationError>();
			ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_INSERT_SIZE, "Insert size out of range",
					getReadName()));
		}
		if (getReadUnmappedFlag()) {
			if (getNotPrimaryAlignmentFlag()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_NOT_PRIM_ALIGNMENT,
						"Not primary alignment flag should not be set for unmapped read.", getReadName()));
			}

		} else {
			if (getMappingQuality() >= 256) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_MAPPING_QUALITY, "MAPQ should be < 256.",
						getReadName()));
			}
			if (getCigarLength() == 0) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_CIGAR,
						"CIGAR should have > zero elements for mapped read.", getReadName()));
			}
			if (getHeader().getSequenceDictionary().size() == 0) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.MISSING_SEQUENCE_DICTIONARY,
						"Empty sequence dictionary.", getReadName()));
			}
			if (!hasReferenceName()) {
				if (ret == null)
					ret = new ArrayList<SAMValidationError>();
				ret.add(new SAMValidationError(SAMValidationError.Type.INVALID_FLAG_READ_UNMAPPED,
						"Mapped read should have valid reference name", getReadName()));
			}

		}
		final String rgId = (String) getAttribute(SAMTagUtil.getSingleton().RG);
		if (rgId != null && getHeader().getReadGroup(rgId) == null) {
			if (ret == null)
				ret = new ArrayList<SAMValidationError>();
			ret.add(new SAMValidationError(SAMValidationError.Type.READ_GROUP_NOT_FOUND,
					"RG ID on SAMRecord not found in header: " + rgId, getReadName()));
		}
		final List<SAMValidationError> errors = isValidReferenceIndexAndPosition(mReferenceIndex,
				getMateReferenceName(), getAlignmentStart(), false);
		if (errors != null) {
			if (ret == null)
				ret = new ArrayList<SAMValidationError>();
			ret.addAll(errors);
		}
		if (this.getReadLength() == 0 && !this.getNotPrimaryAlignmentFlag()) {
			Object fz = getAttribute(SAMTagUtil.getSingleton().FZ);
			if (fz == null) {
				String cq = (String) getAttribute(SAMTagUtil.getSingleton().CQ);
				String cs = (String) getAttribute(SAMTagUtil.getSingleton().CS);
				if (cq == null || cq.length() == 0 || cs == null || cs.length() == 0) {
					if (ret == null)
						ret = new ArrayList<SAMValidationError>();
					ret.add(new SAMValidationError(SAMValidationError.Type.EMPTY_READ,
							"Zero-length read without FZ, CS or CQ tag", getReadName()));
				} else if (!getReadUnmappedFlag()) {
					boolean hasIndel = false;
					for (CigarElement cigarElement : getCigar().getCigarElements()) {
						if (cigarElement.getOperator() == CigarOperator.DELETION
								|| cigarElement.getOperator() == CigarOperator.INSERTION) {
							hasIndel = true;
							break;
						}
					}
					if (!hasIndel) {
						if (ret == null)
							ret = new ArrayList<SAMValidationError>();
						ret.add(new SAMValidationError(SAMValidationError.Type.EMPTY_READ,
								"Colorspace read with zero-length bases but no indel", getReadName()));
					}
				}
			}
		}
		if (this.getReadLength() != getBaseQualities().length && !Arrays.equals(getBaseQualities(), NULL_QUALS)) {
			if (ret == null)
				ret = new ArrayList<SAMValidationError>();
			ret.add(new SAMValidationError(SAMValidationError.Type.MISMATCH_READ_LENGTH_AND_QUALS_LENGTH,
					"Read length does not match quals length", getReadName()));
		}
		if (ret == null || ret.size() == 0) {
			return null;
		}
		return ret;
	}

	// /////////////////////////////////////////////////////////////////////////////
	// *** ReduceReads functions ***//
	// /////////////////////////////////////////////////////////////////////////////

	public byte[] getReducedReadCounts() {

		byte[] reducedReadCounts = getByteArrayAttribute(REDUCED_READ_CONSENSUS_TAG);

		return reducedReadCounts;
	}

	public boolean isReducedRead() {
		return getReducedReadCounts() != null;
	}

	/**
	 * Calculates the reference coordinate for the beginning of the read taking
	 * into account soft clips but not hard clips.
	 *
	 * Note: getUnclippedStart() adds soft and hard clips, this function only
	 * adds soft clips.
	 */
	public int getSoftStart() {
		if (softStart < 0) {
			int start = this.getUnclippedStart();
			for (CigarElement cigarElement : this.getCigar().getCigarElements()) {
				if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
					start += cigarElement.getLength();
				else
					break;
			}
			softStart = start;
		}
		return softStart;
	}

	/**
	 * Calculates the reference coordinate for the end of the read taking into
	 * account soft clips but not hard clips.
	 *
	 * Note: getUnclippedEnd() adds soft and hard clips, this function only adds
	 * soft clips.
	 */
	public int getSoftEnd() {
		if (softEnd < 0) {
			int stop = this.getUnclippedStart();

			if (ReadUtils.readIsEntirelyInsertion(this))
				return stop;

			int shift = 0;
			CigarOperator lastOperator = null;
			for (CigarElement cigarElement : this.getCigar().getCigarElements()) {
				stop += shift;
				lastOperator = cigarElement.getOperator();
				if (cigarElement.getOperator().consumesReferenceBases()
						|| cigarElement.getOperator() == CigarOperator.SOFT_CLIP
						|| cigarElement.getOperator() == CigarOperator.HARD_CLIP)
					shift = cigarElement.getLength();
				else
					shift = 0;
			}
			softEnd = (lastOperator == CigarOperator.HARD_CLIP) ? stop - 1 : stop + shift - 1;
		}
		return softEnd;
	}

	/**
	 * The number of bases corresponding the i'th base of the reduced read.
	 */
	public final byte getReducedCount(final int i) {
		byte firstCount = getReducedReadCounts()[0];
		byte offsetCount = getReducedReadCounts()[i];
		return (i == 0) ? firstCount : (byte) Math.min(firstCount + offsetCount, Byte.MAX_VALUE);
	}

	public byte[] getBaseQualities(final EventType errorModel) {
		switch (errorModel) {
		case SNP:
			return getBaseQualities();
		case Insertion:
			return getBaseInsertionQualities(true);
		case Deletion:
			return getBaseDeletionQualities(true);
		default:
			throw new RuntimeException("Unrecognized Base Recalibration type: " + errorModel);
		}
	}

	public byte[] getBaseDeletionQualities() {
		return getBaseDeletionQualities(false);
	}

	public byte[] getBaseDeletionQualities(boolean disableAttribute) {
		byte[] quals = getExistingBaseDeletionQualities();
		if (quals == null) {
			quals = new byte[getBaseQualities().length];
			Arrays.fill(quals, (byte) 45);
			if (!disableAttribute)
				setBaseQualities(quals, EventType.Deletion);
		}
		return quals;
	}

	public byte[] getExistingBaseDeletionQualities() {
		return SAMUtils.fastqToPhred(getStringAttribute(BQSR_BASE_DELETION_QUALITIES));
	}

	public byte[] getBaseInsertionQualities() {
		return getBaseInsertionQualities(false);
	}

	public byte[] getBaseInsertionQualities(boolean disableAttribute) {
		byte[] quals = getExistingBaseInsertionQualities();
		if (quals == null) {
			quals = new byte[getBaseQualities().length];
			Arrays.fill(quals, (byte) 45);
			if (!disableAttribute)
				setBaseQualities(quals, EventType.Insertion);
		}
		return quals;
	}

	public byte[] getExistingBaseInsertionQualities() {
		return SAMUtils.fastqToPhred(getStringAttribute(BQSR_BASE_INSERTION_QUALITIES));
	}

	public void setBaseQualities(final byte[] quals, final EventType errorModel) {
		switch (errorModel) {
		case SNP:
			setBaseQualities(quals);
			break;
		case Insertion:
			setAttribute(GaeaSamRecord.BQSR_BASE_INSERTION_QUALITIES,
					quals == null ? null : SAMUtils.phredToFastq(quals));
			break;
		case Deletion:
			setAttribute(GaeaSamRecord.BQSR_BASE_DELETION_QUALITIES,
					quals == null ? null : SAMUtils.phredToFastq(quals));
			break;
		default:
			throw new RuntimeException("Unrecognized Base Recalibration type: " + errorModel);
		}
	}

	public NGSPlatform getNGSPlatform() {
		NGSPlatform mNGSPlatform = NGSPlatform.fromReadGroupPL(getReadGroup().getPlatform());

		return mNGSPlatform;
	}

	public boolean containsTemporaryAttribute(Object key) {
		if (temporaryAttributes != null) {
			return temporaryAttributes.containsKey(key);
		}
		return false;
	}

	public Object setTemporaryAttribute(Object key, Object value) {
		if (temporaryAttributes == null) {
			temporaryAttributes = new HashMap<Object, Object>();
		}
		return temporaryAttributes.put(key, value);
	}

	public Object getTemporaryAttribute(Object key) {
		if (temporaryAttributes != null) {
			return temporaryAttributes.get(key);
		}
		return null;
	}

	public boolean isEmpty() {
		return super.getReadBases() == null || super.getReadLength() == 0;
	}

	public static GaeaSamRecord emptyRead(GaeaSamRecord read) {
		GaeaSamRecord emptyRead = new GaeaSamRecord(read.getHeader(), read.getReferenceIndex(), 0, (short) 0,
				read.getFlags(), read.getMateReferenceIndex(), read.getMateAlignmentStart(),
				read.getInferredInsertSize());
		emptyRead.setReferenceIndex(read.getReferenceIndex());

		emptyRead.setCigarString("");
		emptyRead.setReadBases(new byte[0]);
		emptyRead.setBaseQualities(new byte[0]);

		SAMReadGroupRecord samRG = read.getReadGroup();
		emptyRead.clearAttributes();
		if (samRG != null) {
			SAMReadGroupRecord rg = new SAMReadGroupRecord(samRG.getReadGroupId(), samRG);
			emptyRead.setReadGroup(rg);
		}

		return emptyRead;
	}

	public void setReadGroup(final SAMReadGroupRecord readGroup) {
		setAttribute("RG", readGroup.getId());
	}

	public void resetSoftStartAndEnd() {
		softStart = -1;
		softEnd = -1;
	}

	public boolean hasBaseIndelQualities() {
		return getAttribute(BQSR_BASE_INSERTION_QUALITIES) != null
				|| getAttribute(BQSR_BASE_DELETION_QUALITIES) != null;
	}

	public void setPosition(final String contig, final int start) {
		if (contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1) {
			throw new IllegalArgumentException(
					"contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME
							+ ", and start must be >= 1 \ncontig = " + contig + "\nstart = " + start);
		}

		setReferenceName(contig);
		setAlignmentStart(start);
		setReadUnmappedFlag(false);
	}

	public boolean isUnmapped() {
		return getReadUnmappedFlag() || getReferenceName() == null
				|| getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
				|| getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
	}

	public byte getBaseQuality(int index) {
		if (index < 0 || index >= getReadLength())
			throw new IllegalArgumentException(
					"index " + index + " must more than zero and less than " + getReadLength());
		return this.getBaseQualities()[index];
	}

	public boolean mateIsUnmapped() {
		if (!getReadPairedFlag())
			throw new IllegalStateException("Cannot get mate information for an unpaired read");

		return getMateUnmappedFlag() || getMateReferenceName() == null
				|| getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
				|| getMateAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
	}

	public int getMateStart() {
		if (mateIsUnmapped()) {
			return 0;
		}

		return getMateAlignmentStart();
	}
	
	public String getMateContig() {
        if ( mateIsUnmapped() ) {
            return null;
        }

        return getMateReferenceName();
    }
	
	public boolean readHasNoAssignedPosition() {
	    // Check actual assigned positions rather than unmapped status, so that unmapped reads with
	    // assigned positions will be considered to have a position
	    return this.getReferenceName() == null ||
	           this.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
	           this.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
	}
	
	@Override
	public int getEnd() {
		if ( isUnmapped() ) {
            return UNSET_POSITION;
        }
		
        return getAlignmentEnd();
    }
	
	public boolean isFirstOfPair() {
		return getFirstOfPairFlag();
	}
	
	public Integer getAttributeAsInteger( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof Integer ) {
            return (Integer)attributeValue;
        }
        else {
            try {
                return Integer.parseInt(attributeValue.toString());
            }
            catch ( NumberFormatException e ) {
                throw new UserException.ReadAttributeTypeMismatch(attributeName, "integer", e);
            }
        }
    }

    public String getAttributeAsString( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = getAttribute(attributeName);
        if ( attributeValue instanceof byte[]) {
            // in case that the attribute is a byte[] array, the toString method will format it as name@hashCode
            // for a good representation of the byte[] as String, it encodes the bytes with the default charset (UTF-8)
            final byte[] val = (byte[]) attributeValue;
            return (val.length == 0) ? "" : new String(val, DEFAULT_CHARSET);
        }
        // otherwise, just use the toString() method unless it is null
        return attributeValue != null ? attributeValue.toString() : null;
    }
    
    public byte[] getBases() {
        final byte[] bases = getReadBases();

        // Make a defensive copy to protect against direct modification of the returned array
        return bases != null ? Arrays.copyOf(bases, bases.length) : new byte[0];
    }
    
    public String commonToString() {
        //Note: SAMRecord blows up on getAlignmentEnd when cigar is null.
        // That would result in a blow up here so we work around this bug
        // by checking for empty cigar (nulls get converted to empty cigars in SAMRecordToGATKReadAdapter)
        if (isUnmapped() || getCigar().isEmpty()){
            return String.format("%s UNMAPPED", this.getReadName());
        } else {
            return String.format("%s %s:%d-%d", getReadName(), getContig(), getStart(), getEnd());
        }
    }
    
    public boolean hasAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        return getAttribute(attributeName) != null;
    }
}
