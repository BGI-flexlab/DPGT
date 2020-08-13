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

import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;

/**
 * HGSV notation
 *
 * References: http://www.hgvs.org/
 *
 * @author pcingola
 */
public class Hgvs {

	// Don't show sequences that are too long
	public static final int MAX_SEQUENCE_LEN_HGVS = 100;

	protected VariantEffect variantEffect;
	protected Variant variant;
	protected Marker marker;
	protected Transcript tr;
	protected Genome genome;

	protected boolean duplication;
	protected boolean strandPlus, strandMinus;
	protected boolean hgvsTrId;

	public static String parseTranscript(String hgvs) {
		int idxTr = hgvs.indexOf(':');
		if (idxTr < 0) return null;
		return hgvs.substring(0, idxTr);
	}

	public static String removeTranscript(String hgvs) {
		int idxTr = hgvs.indexOf(':');
		if (idxTr < 0) return hgvs;
		return hgvs.substring(idxTr + 1);
	}

	public Hgvs(VariantEffect variantEffect) {
		this.variantEffect = variantEffect;
		variant = variantEffect.getVariant();
		marker = variantEffect.getMarker();
		tr = variantEffect.getTranscript();
		genome = marker != null ? marker.getGenome() : null;
		hgvsTrId = Config.get().isHgvsTrId();
		initStrand();
	}

	protected void initStrand() {
		// Strand information
		if (tr != null) strandMinus = tr.isStrandMinus();
		else if (marker != null) strandMinus = marker.isStrandMinus();

		strandPlus = !strandMinus;
	}

}
