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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.effect;

import org.bgi.flexlab.gaea.tools.annotator.interval.Exon;
import org.bgi.flexlab.gaea.tools.annotator.interval.Intron;
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;
import org.bgi.flexlab.gaea.tools.annotator.util.GprSeq;

/**
 * Coding DNA reference sequence
 *
 * References http://www.hgvs.org/mutnomen/recs.html
 *
 * Nucleotide numbering:
 * 	- there is no nucleotide 0
 * 	- nucleotide 1 is the A of the ATG-translation initiation codon
 * 	- the nucleotide 5' of the ATG-translation initiation codon is -1, the previous -2, etc.
 * 	- the nucleotide 3' of the translation stop codon is *1, the next *2, etc.
 * 	- intronic nucleotides (coding DNA reference sequence only)
 * 		- beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron, like c.77+1G, c.77+2T, ....
 * 		- end of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron, like ..., c.78-2A, c.78-1G.
 * 		- in the middle of the intron, numbering changes from "c.77+.." to "c.78-.."; for introns with an uneven number of nucleotides the central nucleotide is the last described with a "+" (see Discussion)
 *
 * Genomic reference sequence
 * 		- nucleotide numbering starts with 1 at the first nucleotide of the sequence
 * 		  NOTE: the sequence should include all nucleotides covering the sequence (gene) of interest and should start well 5' of the promoter of a gene
 * 		- no +, - or other signs are used
 * 		- when the complete genomic sequence is not known, a coding DNA reference sequence should be used
 * 		- for all descriptions the most 3' position possible is arbitrarily assigned to have been changed (see Exception)
 */

public class HgvsDna extends Hgvs {

	public static boolean debug = false;

	public HgvsDna(VariantEffect variantEffect) {
		super(variantEffect);
	}

	/**
	 * DNA level base changes
	 */
	protected String dnaBaseChange() {

		switch (variant.getVariantType()) {
		case SNP:
			if (strandPlus) return variant.getReference() + ">" + variant.getAlt();
			return GprSeq.reverseWc(variant.getReference()) + ">" + GprSeq.reverseWc(variant.getAlt());

		case MNP:
			String ref, alt;
			if (strandPlus) {
				ref = variant.getReference();
				alt = variant.getAlt();
			} else {
				ref = GprSeq.reverseWc(variant.getReference());
				alt = GprSeq.reverseWc(variant.getAlt());
			}
			return "del" + ref + "ins" + alt;

		case DEL:
		case DUP:
		case INS:
			if (variant.size() > MAX_SEQUENCE_LEN_HGVS) return "";
			String netChange = variant.netChange(false);
			if (strandPlus) return netChange;
			return GprSeq.reverseWc(netChange);

		case MIXED:
			if (strandPlus) return "del" + variant.getReference() + "ins" + variant.getAlt();
			return "del" + GprSeq.reverseWc(variant.getReference()) + "ins" + GprSeq.reverseWc(variant.getAlt());

		case INV:
			// Inversions are designated by "inv" after an indication of the
			// first and last nucleotides affected by the inversion.
			// Reference: http://www.hgvs.org/mutnomen/recs-DNA.html#inv
			// => No base changes are used
			return "";

		case INTERVAL:
			return "";

		case BND:
			return "";

		default:
			throw new RuntimeException("Unimplemented method for variant type " + variant.getVariantType());
		}
	}

	/**
	 * Is this position downstream?
	 */
	boolean isDownstream(int pos) {
		if (tr.isStrandPlus()) return tr.getEnd() < pos;
		return pos < tr.getStart();
	}

	/**
	 * Is this a duplication?
	 */
	protected boolean isDuplication() {
		// Only insertion cause duplications
		// Reference: Here is a discussion of a possible new term ('los') as an analogous
		//            to 'dup' for deletions:
		//                http://www.hgvs.org/mutnomen/disc.html#loss
		//            So, it's still not decided if there is an analogous 'dup' term
		//            for deletions.
		if (!variant.isIns()) return false;

		// Extract sequence from genomic coordinates before variant
		String seq = null;

		// Get sequence at the 3'end of the variant
		int sstart, send;
		int len = variant.getAlt().length();

		if (strandPlus) {
			sstart = Math.max(0, variant.getStart() - len);
			send = variant.getStart() - 1;
		} else {
			sstart = variant.getStart();
			send = sstart + (len - 1);
		}

		// Maybe we can just use exonic sequences (it's faster)
		// Create a marker and check that is comprised within exon boundaries
		Marker m = new Marker(variant.getParent(), sstart, send, false, "");

//		Exon ex = variantEffect.getExon();
//		if (ex != null && ex.includes(m)) {
//			if (debug) Gpr.debug("Variant: " + variant + "\n\tmarker: " + m.toStr() + "\tsstart:" + sstart + "\tsend: " + send + "\n\texon: " + ex + "\n\tstrand: " + (strandPlus ? "+" : "-"));
//			seq = ex.getSequence(m);
//			if (debug) Gpr.debug("Sequence (Exon)  [ " + sstart + " , " + send + " ]: '" + seq + "'\talt: '" + variant.getAlt() + "'\tsequence (+ strand): " + (ex.isStrandPlus() ? ex.getSequence() : GprSeq.reverseWc(ex.getSequence())));
//		}

		// May be it is not completely in the exon. Use genomic sequences
		if (seq == null) {
			seq = genome.querySequence(m);
			if (debug) Gpr.debug("Sequence (Genome) [ " + sstart + " , " + send + " ]: '" + seq + "'\talt: '" + variant.getAlt() + "'\tsequence (+ strand): " + seq);
		}

		// Compare to ALT sequence
		if (seq == null) return false; // Cannot compare

		return seq.equalsIgnoreCase(variant.getAlt());
	}

	/**
	 * Is this position upstream?
	 */
	boolean isUpstream(int pos) {
		if (tr.isStrandPlus()) return pos < tr.getStart();
		return tr.getEnd() < pos;
	}

	/**
	 * Genomic position for exonic variants
	 */
	protected String pos() {
		int posStart = -1, posEnd = -1;

		int variantPosStart = strandPlus ? variant.getStart() : variant.getEnd();

		switch (variant.getVariantType()) {
		case SNP:
			posStart = posEnd = variantPosStart;
			break;

		case MNP:
			posStart = variantPosStart;
			posEnd = posStart + (strandPlus ? 1 : -1) * (variant.size() - 1);
			break;

		case INS:
			posStart = variantPosStart;
			if (duplication) {
				// Duplication coordinates
				int lenAlt = variant.getAlt().length();
				if (lenAlt == 1) {
					// One base duplications do not require end positions:
					// Reference: http://www.hgvs.org/mutnomen/disc.html#dupins
					// Example: c.7dupT (or c.7dup) denotes the duplication (insertion) of a T at position 7 in the sequence ACTTACTGCC to ACTTACTTGCC
					posStart += strandPlus ? -1 : 0; // The 'previous base' is duplicated, we have to decrement the position
					posEnd = posStart;
				} else {
					// Duplication coordinates
					if (strandPlus) {
						posEnd = posStart - 1;
						posStart -= lenAlt;
					} else {
						// Insert is 'before' variant position, so we must shift one base (compared to plus strand)
						posEnd = posStart;
						posStart += lenAlt - 1;
					}
				}
			} else {
				// Other insertions must list both positions:
				// Reference: http://www.hgvs.org/mutnomen/disc.html#ins
				//            ...to prevent confusion, both flanking residues have to be listed.
				// Example: c.6_7dup (or c.6_7dupTG) denotes a TG duplication (TG insertion) in the sequence ACATGTGCC to ACATGTGTGCC
				if (strandPlus) {
					// Insert before current posStart
					posEnd = posStart;
					posStart -= 1;
				} else {
					// Insert before current posStart (negative strand)
					posEnd = posStart - 1;
				}
			}
			break;

		case BND:
		case DEL:
		case DUP:
		case INV:
		case MIXED:
			if (strandPlus) {
				posStart = variant.getStart();
				posEnd = variant.getEnd();
			} else {
				posStart = variant.getEnd();
				posEnd = variant.getStart();
			}
			break;

		case INTERVAL:
			return "";

		default:
			throw new RuntimeException("Unimplemented method for variant type " + variant.getVariantType());
		}

		// Single base
		if (posStart == posEnd) return pos(posStart);

		// Base range
		String ps = pos(posStart);
		String pe = pos(posEnd);
		if (ps == null || pe == null) return null;
		return ps + "_" + pe;
	}

	/**
	 * HGVS position base on genomic coordinates (chr is assumed to be the same as in transcript/marker).
	 */
	protected String pos(int pos) {
		// Cannot do much if there is no transcript
		if (tr == null) return Integer.toString(pos + 1);

		// Are we in an exon?
		// Note: This may come from an intron-exon boundary variant (intron side, walked in a duplication).
		//       In that case, the exon marker won't be available from 'variantEffect.marker'.
		Exon ex = tr.findExon(pos);
		if (ex != null) return posExon(pos);

		Intron intron = tr.findIntron(pos);
		if (intron != null) return posIntron(pos, intron);

		if (isDownstream(pos)) return posDownstream(pos);
		if (isUpstream(pos)) return posUpstream(pos);

		if (debug) Gpr.debug("Unknown HGVS position " + pos + ", transcript " + tr);
		return null;
	}

	/**
	 * Position downstream of the transcript
	 */
	protected String posDownstream(int pos) {
		int baseNumCdsEnd = tr.getCdsEnd();
		int idx = Math.abs(pos - baseNumCdsEnd);

		return "*" + idx; // We are after stop codon, coordinates must be '*1', '*2', etc.
	}

	/**
	 * Convert genomic position to HGVS compatible (DNA) position
	 */
	protected String posExon(int pos) {
		if (tr.isUtr3(pos)) return posUtr3(pos);
		if (tr.isUtr5(pos)) return posUtr5(pos);

		int idx = tr.baseNumberCds(pos, false) + 1; // Coding Exon: just use CDS position

		// Could not find dna position in transcript?
		if (idx <= 0) return null;
		return "" + idx;
	}

	/**
	 * Intronic position
	 */
	protected String posIntron(int pos, Intron intron) {
		// Jump to closest exon position
		// Ref:
		//		beginning of the intron; the number of the last nucleotide of the preceding exon, a plus sign and the position in the intron, like c.77+1G, c.77+2T, etc.
		// 		end of the intron; the number of the first nucleotide of the following exon, a minus sign and the position upstream in the intron, like c.78-1G.
		int posExon = -1;
		String posExonStr = "";
		int distanceLeft = Math.max(0, pos - intron.getStart()) + 1;
		int distanceRight = Math.max(0, intron.getEnd() - pos) + 1;
		if (distanceLeft < distanceRight) {
			posExon = intron.getStart() - 1;
			posExonStr = (intron.isStrandPlus() ? "+" : "-");
		} else if (distanceRight < distanceLeft) {
			posExon = intron.getEnd() + 1;
			posExonStr = (intron.isStrandPlus() ? "-" : "+");
		} else {
			// Reference: in the middle of the intron, numbering changes from "c.77+.." to "c.78-.."; for introns with an uneven number of nucleotides the central nucleotide is the last described with a "+"
			posExonStr = "+";

			if (strandPlus) posExon = intron.getStart() - 1;
			else posExon = intron.getEnd() + 1;
		}

		// Distance to closest exonic base
		int exonDistance = Math.abs(posExon - pos);

		// Closest exonic base within coding region?
		int cdsLeft = Math.min(tr.getCdsStart(), tr.getCdsEnd());
		int cdsRight = Math.max(tr.getCdsStart(), tr.getCdsEnd());
		if ((posExon >= cdsLeft) && (posExon <= cdsRight)) {
			int distExonBase = tr.baseNumberCds(posExon, false) + 1;
			return distExonBase + (exonDistance > 0 ? posExonStr + exonDistance : "");
		}

		// Left side of coding part
		int cdnaPos = tr.baseNumber2MRnaPos(posExon);
		if (posExon < cdsLeft) {
			int cdnaStart = tr.baseNumber2MRnaPos(cdsLeft); // tr.getCdsStart());
			int utrDistance = Math.abs(cdnaStart - cdnaPos);
			String utrStr = strandPlus ? "-" : "*";
			return utrStr + utrDistance + (exonDistance > 0 ? posExonStr + exonDistance : "");
		}

		// Right side of coding part
		int cdnaEnd = tr.baseNumber2MRnaPos(cdsRight); // tr.getCdsEnd());
		int utrDistance = Math.abs(cdnaEnd - cdnaPos);
		String utrStr = strandPlus ? "*" : "-";
		return utrStr + utrDistance + (exonDistance > 0 ? posExonStr + exonDistance : "");
	}

	/**
	 * Position upstream of the transcript
	 */
	protected String posUpstream(int pos) {
		int tss = tr.getCdsStart();
		int idx = Math.abs(pos - tss);

		if (idx <= 0) return null;
		return "-" + idx; // 5'UTR: We are before TSS, coordinates must be '-1', '-2', etc.
	}

	/**
	 * Position within 3'UTR
	 */
	protected String posUtr3(int pos) {
		int baseNum = tr.baseNumber2MRnaPos(pos);
		int baseNumCdsEnd = tr.baseNumber2MRnaPos(tr.getCdsEnd());
		int idx = Math.abs(baseNum - baseNumCdsEnd);

		if (idx <= 0) return null;
		return "*" + idx; // 3'UTR: We are after stop codon, coordinates must be '*1', '*2', etc.
	}

	/**
	 * Position within 5'UTR
	 */
	protected String posUtr5(int pos) {
		int baseNum = tr.baseNumber2MRnaPos(pos);
		int baseNumTss = tr.baseNumber2MRnaPos(tr.getCdsStart());
		int idx = Math.abs(baseNum - baseNumTss);

		if (idx <= 0) return null;
		return "-" + idx; // 5'UTR: We are before TSS, coordinates must be '-1', '-2', etc.
	}

//	/**
//	 * Translocation nomenclature.
//	 * From HGVS:
//	 * 		Translocations are described at the molecular level using the
//	 * 		format "t(X;4)(p21.2;q34)", followed by the usual numbering, indicating
//	 * 		the position translocation breakpoint. The sequences of the translocation
//	 * 		breakpoints need to be submitted to a sequence database (Genbank, EMBL,
//	 * 		DDJB) and the accession.version numbers should be given (see Discussion).
//	 * 		E.g.:
//	 * 			t(X;4)(p21.2;q35)(c.857+101_857+102) denotes a translocation breakpoint
//	 * 			in the intron between coding DNA nucleotides 857+101 and 857+102, joining
//	 * 			chromosome bands Xp21.2 and 4q34
//	 */
//	protected String prefixTranslocation() {
//		VariantTranslocation vtr = (VariantTranslocation) variant;
//
//		// Chromosome part
//		String chrCoords = "(" //
//				+ vtr.getChromosomeName() //
//				+ ";" //
//				+ vtr.getEndPoint().getChromosomeName() //
//				+ ")" //
//				;
//
//		// Get cytobands
//		String band1 = "";
//		CytoBands cytoBands = genome.getCytoBands();
//		Markers bands1 = cytoBands.query(vtr);
//		if (!bands1.isEmpty()) band1 = bands1.get(0).getId(); // Get first match
//
//		String band2 = "";
//		Markers bands2 = cytoBands.query(vtr.getEndPoint());
//		if (!bands2.isEmpty()) band2 = bands2.get(0).getId(); // Get first match
//
//		String bands = "(" + band1 + ";" + band2 + ")";
//
//		return "t" + chrCoords + bands + "(";
//	}

	@Override
	public String toString() {
		if (variant == null || genome == null) return null;

		// Is this a duplication?
		if (variant.isIns()) duplication = isDuplication();

		String type = "", prefix = "", suffix = "";
		switch (variant.getVariantType()) {
		case INS:
			type = duplication ? "dup" : "ins";
			break;

		case DEL:
			type = "del";
			break;

		case MNP:
			type = "";
			break;

		case SNP:
		case MIXED:
		case INTERVAL:
			type = "";
			break;

		case INV:
			type = "inv";
			break;

		case DUP:
			type = "dup";
			break;

		case BND:
//			prefix = prefixTranslocation();
			type = "";
			suffix = ")";
			break;

		default:
			throw new RuntimeException("Unimplemented method for variant type " + variant.getVariantType());
		}

		// HGVS formatted Position
		String pos = pos();
		if (pos == null) return null;

		return prefix + typeOfReference() + pos + type + dnaBaseChange() + suffix;
	}

	/**
	 * Prefix for coding or non-coding sequences
	 */
	protected String typeOfReference() {
		if (tr == null) return "n.";

		// Is the transcript protein coding?
		String prefix = tr.isProteinCoding() ? "c." : "n.";

		// Not using transcript ID?
		if (!hgvsTrId) return prefix;

		StringBuilder sb = new StringBuilder();
		sb.append(tr.getId());

		String ver = tr.getVersion();
		if (!ver.isEmpty()) sb.append("." + ver);

		sb.append(':');
		sb.append(prefix);

		return sb.toString();
	}

}
