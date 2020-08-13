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

public class HgvsOld extends HgvsDna {

	public HgvsOld(VariantEffect variantEffect) {
		super(variantEffect);
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

//		if (isDownstream(pos)) return posDownstream(pos);
//		if (isUpstream(pos)) return posUpstream(pos);

		if (debug) Gpr.debug("Unknown HGVS position " + pos + ", transcript " + tr);
		return null;
	}

	/**
	 * Convert genomic position to HGVS compatible (DNA) position
	 */
	protected String posExon(int pos) {
		int idx = tr.baseNumberTr(pos, false) + 1; // Coding Exon: just use CDS position

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
			int distExonBase = tr.baseNumberTr(posExon, false) + 1;
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

		String prefix = "n.";

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
