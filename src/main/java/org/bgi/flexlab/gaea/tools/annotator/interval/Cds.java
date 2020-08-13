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

import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;

/**
 * CDS: The coding region of a gene, also known as the coding sequence or CDS (from Coding DNA Sequence), is
 * that portion of a gene's DNA or RNA, composed of exons, that codes for protein.
 *
 * @author pcingola
 *
 */
public class Cds extends Marker implements MarkerWithFrame {

	private static final long serialVersionUID = 1636197649250882952L;

	byte frame = -1; // Frame can be {-1, 0, 1, 2}, where '-1' means unknown

	public Cds() {
		super();
		type = EffectType.CDS;
	}

	public Cds(Transcript parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.CDS;
	}

	@Override
	public Cds cloneShallow() {
		Cds clone = (Cds) super.cloneShallow();
		clone.frame = frame;
		return clone;
	}

	/**
	 * Correct coordinates according to frame differences
	 * @param frameCorrection
	 */
	public boolean frameCorrection(int frameCorrection) {
		if (frameCorrection <= 0) return true; // Nothing to do

		// Can correct?
		if (size() <= frameCorrection) {
			Config.get().warning("CDS too short, cannot correct frame:", " frame size " + size() + ", frame correction " + frameCorrection + ", CDS: " + this);
			return false;
		}

		// Correct start or end coordinates
		if (isStrandPlus()) start += frameCorrection;
		else end -= frameCorrection;

		// Correct frame
		frame = (byte) ((frame - frameCorrection) % 3);
		while (frame < 0)
			frame += 3;

		return true;
	}

	@Override
	public int getFrame() {
		return frame;
	}


	/**
	 * Frame can be {-1, 0, 1, 2}, where '-1' means unknown
	 * @param frame
	 */
	@Override
	public void setFrame(int frame) {
		if ((frame > 2) || (frame < -1)) throw new RuntimeException("Invalid frame value: " + frame);
		this.frame = (byte) frame;
	}

	@Override
	public String toString() {
		return getChromosomeName() + "\t" + start + "-" + end //
				+ " " //
				+ type //
				+ ((id != null) && (id.length() > 0) ? " '" + id + "'" : "") //
				+ ", frame: " + frame //
				;
	}

}
