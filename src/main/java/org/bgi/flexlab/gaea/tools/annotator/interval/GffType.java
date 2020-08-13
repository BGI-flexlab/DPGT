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
package org.bgi.flexlab.gaea.tools.annotator.interval;

public enum GffType {

	GENE //
	, TRANSCRIPT //
	, EXON //
	, CDS //
	, START_CODON //
	, STOP_CODON //
	, UTR5 //
	, UTR3 //
	, INTRON_CONSERVED //
	, INTERGENIC_CONSERVED //
	, UNKNOWN //
	;

	public static GffType parse(String str) {
		switch (str.toLowerCase()) {
		case "gene":
		case "protein":
			return GENE;

		case "pseudogene":
		case "transcript":
		case "mrna":
		case "trna":
		case "snorna":
		case "rrna":
		case "ncrna":
		case "mirna":
		case "snrna":
		case "pseudogenic_transcript":
			return TRANSCRIPT;

		case "exon":
		case "pseudogenic_exon":
			return EXON;

		case "cds":
			return CDS;

		case "start_codon":
			return START_CODON;

		case "stop_codon":
			return STOP_CODON;

		case "five_prime_utr":
		case "5'-utr":
		case "5'utr":
		case "5utr":
			return UTR5;

		case "three_prime_utr":
		case "3'-utr":
		case "3'utr":
		case "3utr":
			return UTR3;

		case "intron_CNS":
		case "intron_cns":
			return INTRON_CONSERVED;

		case "inter_cns":
			return INTERGENIC_CONSERVED;

		default:
			return UNKNOWN;
		}
	}
}
