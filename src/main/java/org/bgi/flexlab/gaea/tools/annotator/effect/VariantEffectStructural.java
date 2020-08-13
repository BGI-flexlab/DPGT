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

import org.bgi.flexlab.gaea.tools.annotator.interval.*;

import java.util.*;

/**
 * Effect of a structural variant affecting multiple genes
 *
 * @author pcingola
 */
public class VariantEffectStructural extends VariantEffect {

	List<Marker> featuresLeft;
	List<Marker> featuresRight;
	Set<Gene> genes;
	int countWholeGenes = 0; // How many genes does the variant fully include?
	int countPartialGenes = 0; // How many genes does the variant partially overlap?

	public VariantEffectStructural(Variant variant) {
		this(variant, null);
	}

	public VariantEffectStructural(Variant variant, Markers intersects) {
		super(variant);
		featuresLeft = new LinkedList<>();
		featuresRight = new LinkedList<>();
		genes = new HashSet<Gene>();

		if (intersects != null) {
			setGenes(intersects);
			setEffect(effect());
		}
	}

	EffectType effect() {
		switch (variant.getVariantType()) {
		case INV:
			return countWholeGenes > 0 ? EffectType.GENE_INVERSION : EffectType.NONE;

		case DEL:
			// Note: For one gene, we annotate all transcripts
			return countWholeGenes > 1 ? EffectType.GENE_DELETED : EffectType.NONE;

		case DUP:
			return countWholeGenes > 0 ? EffectType.GENE_DUPLICATION : EffectType.NONE;

		case BND:
			// Translocation can only intersect "partial genes" since they 
			// have two (unconnected) break-points. 
			switch (countPartialGenes) {
			case 0:
				// We use 'FEATURE_FUSION'  when neither end of the translocation 
				// lands on a gene (i.e.a fusion of intergenic regions)
				return EffectType.FEATURE_FUSION;

			case 1:
				// Genes on only one side of the translocation (the other side is intergenic) 
				return EffectType.GENE_FUSION_HALF;

			default:
				// Genes on both sides of the translocation
				return EffectType.GENE_FUSION;
			}

		default:
			throw new RuntimeException("Unknown effect for variant type " + variant.getVariantType());
		}

	}

	/**
	 * Is there another 'fusion' effect?
	 */
	public List<VariantEffect> fusion() {
		// Only if both genes are different
		if (variant.isBnd()) {
			if (featuresLeft.isEmpty() || featuresRight.isEmpty()) return null;

			if (featuresLeft.isEmpty()) { return null; }

			if (featuresRight.isEmpty()) return null;
		} else {
			if (featuresLeft.isEmpty() || featuresRight.isEmpty()) return null;
		}

		// Add all gene pairs
		List<VariantEffect> fusions = new LinkedList<VariantEffect>();
		for (Marker gLeft : featuresLeft)
			for (Marker gRight : featuresRight) {
				if (variant.isBnd()) {
					// Is this a translocation? OK
				} else if (!isGene(gLeft) || !isGene(gRight)) {
					// For non-translocations, both sides must be genes in order to create a fusion
					continue;
				} else {
					if (gLeft.getId().equals(gRight.getId())) {
						// Otherwise, make sure the variant is not acting within
						// the same gene (e.g. a deletion)
						continue;
					} else {
						// If both genes overlap and the variant is within that
						// region, then it's not a fusion, it's just a variant
						// acting on both genes.
						Marker gIntersect = gLeft.intersect(gRight);
						if (gIntersect != null && gIntersect.includes(variant)) continue;
					}
				}

				// Add all possible transcript fussions
				fusions.addAll(fusion(variant, gLeft, gRight));
			}

		return fusions;
	}

	/**
	 * Create all possible fusions for these two genes
	 */
	List<VariantEffect> fusion(Variant variant, Marker mLeft, Marker mRight) {
		List<VariantEffect> fusions = new LinkedList<VariantEffect>();
		// One for fusion effect for each transcript
		// This can be a long list...
		Markers msLeft = new Markers();
		Markers msRight = new Markers();

		// If it's a gene, add all transcripts
		Gene gLeft = null;
		if (isGene(mLeft)) {
			gLeft = (Gene) mLeft;
			msLeft.addAll(gLeft.subIntervals());
		} else msLeft.add(mLeft);

		// If it's a gene, add all transcripts
		Gene gRight = null;
		if (isGene(mRight)) {
			gRight = (Gene) mRight;
			msRight.addAll(gRight.subIntervals());
		} else msRight.add(mRight);

		for (Marker ml : msLeft)
			for (Marker mr : msRight) {
				// Transcript from the same gene are added only once
				if (isTranscript(ml) //
						&& isTranscript(mr) //
						&& gLeft.getId().equals(gRight.getId()) // Genes have the same ID
						&& ml.getId().compareTo(mr.getId()) > 0 // Compare transcript IDs alphabetically
				) continue;

				VariantEffectFusion fusion = new VariantEffectFusion(variant, ml, mr);
				fusions.add(fusion);
			}

		return fusions;
	}

	@Override
	public Gene getGene() {
		for (Marker m : featuresLeft)
			if (isGene(m)) return (Gene) m;

		for (Marker m : featuresRight)
			if (isGene(m)) return (Gene) m;

		return null;
	}

	@Override
	public List<Gene> getGenes() {
		ArrayList<Gene> list = new ArrayList<>();
		list.addAll(genes);
		return list;
	}

	@Override
	public Marker getMarker() {
		return getGene();
	}

	/**
	 * We say that intersects the "Left" side if the 'start' of the variant
	 * intersects the marker.
	 *
	 * Note: For translocations this is the 'main' variants (as opposed to
	 * the endPoint). This is just nomenclature and we could have defined
	 * it the other way around (or called it intersect1 and intersect2
	 * instead of intersectLeft intersercRight)
	 */
	boolean intersectsLeft(Marker m) {
		return m.getChromosomeName().equals(variant.getChromosomeName()) //
				&& m.intersects(variant.getStart());
	}

	/**
	 * We say that intersects the "Right" side if the 'end' of the variant
	 * intersects the marker.
	 *
	 * Note: For translocations this is the 'endPoint' (as opposed to
	 * the 'main' variant). This is just nomenclature and we could have defined
	 * it the other way around (or called it intersect1 and intersect2
	 * instead of intersectLeft intersercRight)
	 */
	boolean intersectsRight(Marker m) {
		if (variant.isBnd()) {
			Marker endPoint = ((VariantTranslocation) variant).getEndPoint();
			return m.getChromosomeName().equals(endPoint.getChromosomeName()) //
					&& m.intersects(endPoint.getStart());
		}

		return m.getChromosomeName().equals(variant.getChromosomeName()) //
				&& m.intersects(variant.getEnd());
	}

	protected boolean isGene(Marker m) {
		return m instanceof Gene;
	}

	@Override
	public boolean isMultipleGenes() {
		return true;
	}

	protected boolean isTranscript(Marker m) {
		return m instanceof Transcript;
	}

	/**
	 * Set genes from all intersecting intervals
	 */
	void setGenes(Markers intersects) {
		for (Marker m : intersects)
			if (m instanceof Gene) {
				if (intersectsLeft(m)) featuresLeft.add(m);
				if (intersectsRight(m)) featuresRight.add(m);

				if (variant.includes(m)) countWholeGenes++;
				else countPartialGenes++;

				genes.add((Gene) m);
			} else if (!(m instanceof Chromosome)) {
				if (intersectsLeft(m)) featuresLeft.add(m);
				if (intersectsRight(m)) featuresRight.add(m);
			}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toStr());

		sb.append("\n\tGene left  : [");
		for (Marker m : featuresLeft)
			sb.append(" " + m.getId());
		sb.append("]");

		sb.append("\n\tGene right  : [");
		for (Marker m : featuresRight)
			sb.append(" " + m.getId());
		sb.append("]");

		sb.append("\n\tGenes all: [");
		for (Gene g : genes)
			sb.append(g.getGeneName() + " ");
		sb.append(" ]");

		return sb.toString();
	}
}
