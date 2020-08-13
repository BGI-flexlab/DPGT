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

/**
 * Transcript level support
 * Reference: http://useast.ensembl.org/Help/Glossary?id=492;redirect=no
 *
 * @author pcingola
 */
public enum TranscriptSupportLevel {

	TSL_1 // All splice junctions of the transcript are supported by at least one non-suspect mRNA
	, TSL_2 // The best supporting mRNA is flagged as suspect or the support is from multiple ESTs
	, TSL_3 // The only support is from a single EST
	, TSL_4 // The best supporting EST is flagged as suspect
	, TSL_5 // No single transcript supports the model structure
	, /**
		The transcript was not analyzed for one of the following reasons:
		pseudo-gene annotation, including transcribed pseudo-genes
		human leukocyte antigen (HLA) transcript
		immunoglobin gene transcript
		T-cell receptor transcript
		single-exon transcript (will be included in a future version)
		**/
	TSL_NA;

	public static TranscriptSupportLevel parse(String str) {
		if (str.startsWith("TSL_")) return TranscriptSupportLevel.valueOf(str);

		// Sometimes GTF files have strings like this one:
		//      transcript_support_level "4 (assigned to previous version 3)"
		// So we have to remove the part in parenthesis
		if (str.length() > 2) str = str.substring(0, 2).trim();

		// Safely parse TSL
		try {
			return TranscriptSupportLevel.valueOf("TSL_" + str);
		} catch (Exception e) {
			return null;
		}
	}

}
