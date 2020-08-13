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

import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.EffectImpact;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.ErrorWarningType;
import org.bgi.flexlab.gaea.tools.annotator.interval.Marker;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * A sorted collection of variant effects
 *
 * @author pcingola
 */
public class VariantEffects implements Iterable<VariantEffect> {

	public static boolean debug = false;
	List<VariantEffect> effects;

	public VariantEffects() {
		effects = new ArrayList<VariantEffect>();
	}

	/**
	 * Add an effect
	 */
	public void add(Variant variant, Marker marker, EffectType effectType, EffectImpact effectImpact, String message) {
		VariantEffect effNew = new VariantEffect(variant);
		effNew.set(marker, effectType, effectImpact, message);
		add(effNew);
	}

	/**
	 * Add an effect
	 */
	public void add(Variant variant, Marker marker, EffectType effectType, String message) {
		add(variant, marker, effectType, effectType.effectImpact(), message);
	}

	/**
	 * Add an effect
	 */
	public void add(VariantEffect variantEffect) {
		effects.add(variantEffect);
	}

	/**
	 * Add: If possible, only add an effect type (otherwise add the full effect)
	 */
	public void addEffectType(Variant variant, Marker marker, EffectType effectType) {
		if (canAddType(variant, marker)) {
			get().addEffect(effectType);
		} else{
			add(variant, marker, effectType, effectType.effectImpact(), "");
		}
	}

	public void addErrorWarning(Variant variant, ErrorWarningType errwarn) {
		VariantEffect veff = get();
		if (veff != null) veff.addErrorWarningInfo(errwarn);
		else {
			if (debug) Gpr.debug("Could not get latest " + VariantEffect.class.getSimpleName());
			veff = new VariantEffect(variant);
			veff.addErrorMessage(errwarn);
			add(veff);
		}
	}

	/**
	 * Can we add an effectType to the previous variatnEffect?
	 * @return true if transcript IDs and variant's genotypes match (i.e. we can add effectType)
	 */
	boolean canAddType(Variant variant, Marker marker) {
		VariantEffect veff = get();
		if (veff == null || veff.getVariant() == null) return false;

		// Do genotypes match?
		String gt = veff.getVariant().getGenotype();
		String vgt = variant.getGenotype();
		if (((vgt != null) ^ (gt != null)) // One null and one non-null?
				|| ((vgt != null) && (gt != null) && !variant.getGenotype().equals(variant.getGenotype())) // Both non-null, but different?
		) return false;

		// Do transcripts match?
		Transcript trMarker = (Transcript) marker.findParent(Transcript.class);
		Transcript tr = veff.getTranscript();
		if (tr == null || trMarker == null) return false;

		return tr.getId().equals(trMarker.getId());
	}

	/**
	 * Get (or create) the latest ChangeEffect
	 */
	public VariantEffect get() {
		if (effects.isEmpty()) return null;
		return effects.get(effects.size() - 1);
	}

	public VariantEffect get(int index) {
		return effects.get(index);
	}

	public boolean hasMarker() {
		VariantEffect veff = get();
		if (veff == null) return false;
		return veff.getMarker() != null;
	}

	public boolean isEmpty() {
		return effects.isEmpty();
	}

	@Override
	public Iterator<VariantEffect> iterator() {
		return effects.iterator();
	}

	public void setMarker(Marker marker) {
		VariantEffect veff = get();
		if (veff != null) veff.setMarker(marker);
		else Gpr.debug("Could not get latest VariantEffect");
	}

	public int size() {
		return effects.size();
	}

	public void sort() {
		Collections.sort(effects);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("Effects; " + size() + "\n");
		for (VariantEffect eff : this)
			sb.append(Gpr.prependEachLine("\t", eff.toStr()) + "\n");
		return sb.toString();
	}
}
