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

import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffects;
import org.bgi.flexlab.gaea.tools.annotator.util.KeyValue;

import java.util.Collections;
import java.util.Iterator;

/**
 * This is a custom interval (i.e. intervals provided by the user)
 *
 * @author pcingola
 */
public class Custom extends Marker implements Iterable<KeyValue<String, String>> {

	private static final long serialVersionUID = -6843535415295857726L;

	String label;
	double score = Double.NaN;

	public Custom() {
		super();
		type = EffectType.CUSTOM;
		label = "";
	}

	public Custom(Marker parent, int start, int end, boolean strandMinus, String id, String label) {
		super(parent, start, end, strandMinus, id);
		type = EffectType.CUSTOM;

		this.label = label;
		if (label == null || label.isEmpty()) label = id;
	}

	@Override
	public Custom cloneShallow() {
		Custom clone = (Custom) super.cloneShallow();
		clone.label = label;
		clone.score = score;
		return clone;
	}

	public String getLabel() {
		return label;
	}

	public double getScore() {
		return score;
	}

	/**
	 * Do we have additional annotations?
	 */
	public boolean hasAnnotations() {
		return false;
	}

	@Override
	public Iterator<KeyValue<String, String>> iterator() {
		// Nothing to iterate on
		return Collections.<KeyValue<String, String>> emptySet().iterator();
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public void setScore(double score) {
		this.score = score;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getChromosomeName());
		sb.append("\t");
		sb.append(start);
		sb.append("-");
		sb.append(end);
		sb.append(" ");
		sb.append(type);
		sb.append(((id != null) && (id.length() > 0) ? " '" + id + "'" : ""));

		if (hasAnnotations()) {
			for (KeyValue<String, String> kv : this)
				sb.append(kv.key + "=" + kv.value + ";");
		}

		return sb.toString();
	}

	@Override
	public boolean variantEffect(Variant variant, VariantEffects changeEffecs) {
		if (!intersects(variant)) return false; // Sanity check
		changeEffecs.add(variant, this, EffectType.CUSTOM, label);
		return true;
	}
}
