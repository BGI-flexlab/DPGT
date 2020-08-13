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
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.ErrorWarningType;
import org.bgi.flexlab.gaea.tools.annotator.interval.*;
import org.bgi.flexlab.gaea.tools.annotator.interval.tree.IntervalForest;
import org.bgi.flexlab.gaea.tools.annotator.interval.tree.Itree;

import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Predicts effects of Variants
 *
 * @author pcingola
 */
public class SnpEffectPredictor implements Serializable {
	private static final long serialVersionUID = 4519418862303325081L;

//	public static final int DEFAULT_UP_DOWN_LENGTH = 5000;
	public static final int DEFAULT_UP_DOWN_LENGTH = 0;
	public static final int SMALL_VARIANT_SIZE_THRESHOLD = 10; // Number of bases for a variant to be considered 'small'

	boolean useChromosomes = true;
	boolean debug;
	int upDownStreamLength = DEFAULT_UP_DOWN_LENGTH;
	int spliceSiteSize = SpliceSite.CORE_SPLICE_SITE_SIZE;
	int spliceRegionExonSize = SpliceSite.SPLICE_REGION_EXON_SIZE;
	int spliceRegionIntronMin = SpliceSite.SPLICE_REGION_INTRON_MIN;
	int spliceRegionIntronMax = SpliceSite.SPLICE_REGION_INTRON_MAX;
	Genome genome;
	Markers markers; // All other markers are stored here (e.g. custom markers, intergenic, etc.)
	IntervalForest intervalForest; // Interval forest by chromosome name
	IntervalForest intervalForestGene; // Interval forest by gene


	public SnpEffectPredictor(Genome genome) {
		this.genome = genome;
		markers = new Markers();
	}

	/**
	 * Add a gene interval
	 */
	public void add(Gene gene) {
		genome.getGenes().add(gene);
	}

	/**
	 * Add a marker
	 *
	 * Note: Markers have to be added BEFORE building the interval trees.
	 *       Interval trees are built the first time you call snpEffect(snp) method.
	 */
	public void add(Marker marker) {
		markers.add(marker);
	}

	/**
	 * Add a set of markers
	 */
	public void addAll(Markers markersToAdd) {
		for (Marker marker : markersToAdd)
			markers.add(marker);
	}

	/**
	 * Add a gene dependent marker
	 */
	public void addPerGene(String geneId, Marker marker) {
		if (intervalForestGene == null) {
			intervalForestGene = new IntervalForest();
			intervalForestGene.setDebug(debug);
		}
		intervalForestGene.getOrCreateTreeChromo(geneId).add(marker);
	}

	/**
	 * Create interval trees (forest)
	 */
	public void buildForest() {
		intervalForest = new IntervalForest();
		intervalForest.setDebug(debug);

		// Add all chromosomes to forest
		if (useChromosomes) {
			for (Chromosome chr : genome)
				intervalForest.add(chr);
		}

		//		 TODO 
		// In a circular genome, a gene can have negative coordinates or crosses
		// over chromosome end. These genes are mirrored to the opposite end of
		// the chromosome so that they can be referenced by both circular coordinates.
		//		genome.getGenes().createCircularGenes();

		// Add all genes to forest
		for (Gene gene : genome.getGenes())
			intervalForest.add(gene);

		//---
		// Create (and add) up-down stream, splice sites, intergenic, etc
		//---
		markers.add(createGenomicRegions());

		// Mark canonical transcripts
		canonical();

		// Add all 'markers' to forest (includes custom intervals)
		intervalForest.add(markers);

		// Build interval forest
		intervalForest.build();

		buildPerGene();
	}

	/**
	 * Build 'per gene' information
	 */
	void buildPerGene() {
		if (intervalForestGene != null) intervalForestGene.build();
	}

	/**
	 * Make sure all genes have canonical transcripts
	 */
	void canonical() {
		for (Gene g : genome.getGenes())
			g.canonical();
	}

	/**
	 * Create (and add) up-down stream, splice sites, intergenic, etc
	 */
	public Markers createGenomicRegions() {
		Markers markers = new Markers();

		// Add up-down stream intervals
		for (Marker upDownStream : genome.getGenes().createUpDownStream(upDownStreamLength))
			markers.add(upDownStream);

		// Add splice site intervals
		genome.getGenes().createSpliceSites(spliceSiteSize, spliceRegionExonSize, spliceRegionIntronMin, spliceRegionIntronMax);

		// Intergenic markers
		for (Intergenic intergenic : genome.getGenes().createIntergenic())
			markers.add(intergenic);

		return markers;
	}

	/**
	 * Filter transcripts by TSL
	 */
	public void filterTranscriptSupportLevel(TranscriptSupportLevel maxTsl) {
		for (Gene g : genome.getGenes())
			g.filterTranscriptSupportLevel(maxTsl);
	}

	/**
	 * Obtain a gene interval
	 */
	public Gene getGene(String geneIntervalId) {
		return genome.getGenes().get(geneIntervalId);
	}

	public Genome getGenome() {
		return genome;
	}

	public IntervalForest getIntervalForest() {
		return intervalForest;
	}

	public Markers getMarkers() {
		return markers;
	}

	public int getSpliceRegionExonSize() {
		return spliceRegionExonSize;
	}

	public int getSpliceRegionIntronMax() {
		return spliceRegionIntronMax;
	}

	public int getSpliceRegionIntronMin() {
		return spliceRegionIntronMin;
	}

	public Transcript getTranscript(String trId) {
		for (Gene g : genome.getGenes())
			for (Transcript tr : g)
				if (tr.getId().equals(trId)) return tr;

		return null;
	}

	public int getUpDownStreamLength() {
		return upDownStreamLength;
	}

	/**
	 * Is the chromosome missing in this marker?
	 */
	boolean isChromosomeMissing(Marker marker) {
		// Missing chromosome in marker?
		if (marker.getChromosome() == null) return true;

		// Missing chromosome in genome?
		String chrName = marker.getChromosomeName();
		Chromosome chr = genome.getChromosome(chrName);
		if (chr == null) return true;

		// Chromosome length is 1 or less?
		if (chr.size() < 1) return true;

		// Tree not found in interval forest?
		if (!intervalForest.hasTree(chrName)) return true;

		// OK, we have the chromosome
		return false;
	}

	/**
	 * Dump to sdtout
	 */
	public void print() {
		System.out.println(genome);

		// Show genes
		for (Gene gene : genome.getGenes().sorted())
			System.out.println(gene);

		// Show other inervals
		for (Marker marker : markers)
			System.out.println(marker);
	}

	/**
	 * Return a collection of intervals that intersect 'marker'
	 */
	public Markers query(Marker marker) {
		return marker.query(intervalForest);
	}

	/**
	 * Find closest gene to this marker
	 *
	 * In case more than one 'closest' gene is
	 * found (e.g. two or more genes at the
	 * same distance). The following rules
	 * apply:
	 *
	 * 		i) If many genes have the same 'closest
	 * 		   distance', coding genes are preferred.
	 *
	 * 		ii) If more than one coding gene has the
	 * 		    same 'closet distance', a random gene
	 *			is returned.
	 *
	 * @param inputInterval
	 */
	public Gene queryClosestGene(Marker inputInterval) {
		int initialExtension = 1000;

		String chrName = inputInterval.getChromosomeName();
		Chromosome chr = genome.getChromosome(chrName);
		if (chr == null) return null;

		if (chr.size() > 0) {
			// Extend interval to capture 'close' genes
			for (int extend = initialExtension; extend < chr.size(); extend *= 2) {
				int start = Math.max(inputInterval.getStart() - extend, 0);
				int end = inputInterval.getEnd() + extend;
				Marker extended = new Marker(chr, start, end, false, "");

				// Find all genes that intersect with the interval
				Markers markers = query(extended);
				Markers genes = new Markers();
				int minDist = Integer.MAX_VALUE;
				for (Marker m : markers) {
					if (m instanceof Gene) {
						int dist = m.distance(inputInterval);
						if (dist < minDist) {
							genes.add(m);
							minDist = dist;
						}
					}
				}

				// Found something?
				if (genes.size() > 0) {
					// Find a gene having distance 'minDist'. Prefer coding genes
					Gene minDistGene = null;

					for (Marker m : genes) {
						int dist = m.distance(inputInterval);
						if (dist == minDist) {
							Gene gene = (Gene) m;
							if (minDistGene == null) minDistGene = gene;
							else if (!minDistGene.isProteinCoding() && gene.isProteinCoding()) minDistGene = gene;
						}
					}

					return minDistGene;
				}

			}
		}

		// Nothing found
		return null;
	}

	/**
	 * Return a collection of intervals that intersect 'marker'
	 * Query resulting genes, transcripts and exons to get ALL types of intervals possible
	 */
	public Markers queryDeep(Marker marker) {
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(marker)) throw new RuntimeException("Chromosome missing for marker: " + marker);

		boolean hitChromo = false;
		Markers hits = new Markers();
		Markers intersects = query(marker);

		if (intersects.size() > 0) {
			for (Marker m : intersects) {
				hits.add(m);

				if (m instanceof Chromosome) {
					hitChromo = true; // OK (we have to hit a chromosome, otherwise it's an error
				} else if (m instanceof Gene) {
					// Analyze Genes
					Gene gene = (Gene) m;
					hits.addAll(gene.query(marker));
				}
			}
		}

		if (!hitChromo && Config.get().isErrorChromoHit()) throw new RuntimeException("ERROR: Out of chromosome range. " + marker);
		return hits;
	}

	/**
	 * Name of the regions hit by a marker
	 * @return A set of region names
	 */
	public Set<String> regions(Marker marker, boolean showGeneDetails, boolean compareTemplate) {
		return regions(marker, showGeneDetails, compareTemplate, null);
	}

	/**
	 * Name of the regions hit by a marker
	 * @param id : Only use genes or transcripts matching this ID (null for any)
	 */
	public Set<String> regions(Marker marker, boolean showGeneDetails, boolean compareTemplate, String id) {
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(marker)) throw new RuntimeException("Chromosome missing for marker: " + marker);

		boolean hitChromo = false;
		HashSet<String> hits = new HashSet<String>();

		Markers intersects = query(marker);
		if (intersects.size() > 0) {
			for (Marker markerInt : intersects) {

				if (markerInt instanceof Chromosome) {
					hitChromo = true; // OK (we have to hit a chromosome, otherwise it's an error
					hits.add(markerInt.getClass().getSimpleName()); // Add marker name to the list
				} else if (markerInt instanceof Gene) {
					// Analyze Genes
					Gene gene = (Gene) markerInt;
					regionsAddHit(hits, gene, marker, showGeneDetails, compareTemplate);

					// For all transcripts...
					for (Transcript tr : gene) {
						if ((id == null) || gene.getId().equals(id) || tr.getId().equals(id)) { // Mathes ID? (...or no ID to match)

							// Does it intersect this transcript?
							if (tr.intersects(marker)) {
								regionsAddHit(hits, tr, marker, showGeneDetails, compareTemplate);

								// Does it intersect a UTR?
								for (Utr utr : tr.getUtrs())
									if (utr.intersects(marker)) regionsAddHit(hits, utr, marker, showGeneDetails, compareTemplate);

								// Does it intersect an exon?
								for (Exon ex : tr)
									if (ex.intersects(marker)) regionsAddHit(hits, ex, marker, showGeneDetails, compareTemplate);

								// Does it intersect an intron?
								for (Intron intron : tr.introns())
									if (intron.intersects(marker)) regionsAddHit(hits, intron, marker, showGeneDetails, compareTemplate);
							}
						}
					}
				} else {
					// No ID to match?
					if (id == null) regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate);
					else {
						// Is ID from transcript?
						Transcript tr = (Transcript) markerInt.findParent(Transcript.class);
						if ((tr != null) && (tr.getId().equals(id))) {
							regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate); // Transcript ID matches => count
						} else {
							// Is ID from gene?
							Gene gene = (Gene) markerInt.findParent(Gene.class);
							if ((gene != null) && (gene.getId().equals(id))) regionsAddHit(hits, markerInt, marker, showGeneDetails, compareTemplate); // Gene ID matches => count
						}
					}
				}
			}
		}

		if (!hitChromo) throw new RuntimeException("ERROR: Out of chromosome range. " + marker);
		return hits;
	}

	/**
	 * Add into to a hash
	 */
	void regionsAddHit(HashSet<String> hits, Marker hit2add, Marker marker, boolean showGeneDetails, boolean compareTemplate) {
		String hitStr = hit2add.getClass().getSimpleName();

		if (compareTemplate) {
			Gene gene = (Gene) hit2add.findParent(Gene.class);
			if (gene != null) hitStr += (hit2add.isStrandPlus() == marker.isStrandPlus()) ? "_TEMPLATE_STRAND" : "_NON_TEMPLATE_STRAND";
		}

		if (showGeneDetails && (hit2add instanceof Gene)) {
			Gene gene = (Gene) hit2add;
			hitStr += "[" + gene.getBioType() + ", " + gene.getGeneName() + ", " + (gene.isProteinCoding() ? "protein" : "not-protein") + "]";
		}

		hits.add(hitStr); // Add marker name to the list
	}

	/**
	 * Remove all non-canonical transcripts
	 */
	public void removeNonCanonical() {
		for (Gene g : genome.getGenes())
			g.removeNonCanonical();
	}

	/**
	 * Remove all unverified transcripts
	 *
	 * @return true if ALL genes had ALL transcripts removed (i.e. something
	 * went wrong, like in cases where no transcript was checked during the
	 * building process)
	 */
	public boolean removeUnverified() {
		boolean allRemoved = true;
		for (Gene g : genome.getGenes())
			allRemoved &= g.removeUnverified();

		return allRemoved;
	}

	/**
	 * Remove all transcripts that are NOT in the list
	 * @return : Number of transcripts removed
	 */
	public int retainAllTranscripts(Set<String> trIds) {
		int total = 0;
		for (Gene g : genome.getGenes())
			total += g.keepTranscripts(trIds);
		return total;
	}

	/**
	 * Remove all transcripts that are NOT in the list
	 * @return : Number of transcripts removed
	 */
	public int retainTranscriptsProtein() {
		int total = 0;
		for (Gene g : genome.getGenes())
			total += g.keepTranscriptsProtein();
		return total;
	}


	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setSpliceRegionExonSize(int spliceRegionExonSize) {
		this.spliceRegionExonSize = spliceRegionExonSize;
	}

	public void setSpliceRegionIntronMax(int spliceRegionIntronMax) {
		this.spliceRegionIntronMax = spliceRegionIntronMax;
	}

	public void setSpliceRegionIntronMin(int spliceRegionIntronMin) {
		this.spliceRegionIntronMin = spliceRegionIntronMin;
	}

	public void setSpliceSiteSize(int spliceSiteSize) {
		this.spliceSiteSize = spliceSiteSize;
	}

	public void setUpDownStreamLength(int upDownStreamLength) {
		this.upDownStreamLength = upDownStreamLength;
	}

	public void setUseChromosomes(boolean useChromosomes) {
		this.useChromosomes = useChromosomes;
	}

	public int size() {
		if (intervalForest == null) return 0;
		return intervalForest.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(genome.getVersion() + "\n");
		for (Chromosome chr : genome)
			sb.append(chr + "\n");
		sb.append(genome.getGenes());
		return sb.toString();
	}

	/**
	 * Predict the effect of a variant
	 */
	public VariantEffects variantEffect(Variant variant) {
		VariantEffects variantEffects = new VariantEffects();

		// Chromosome missing?
		if (Config.get().isErrorOnMissingChromo() && isChromosomeMissing(variant)) {
			variantEffects.addErrorWarning(variant, ErrorWarningType.ERROR_CHROMOSOME_NOT_FOUND);
			return variantEffects;
		}

		// Is this a structural variant? Large structural variants (e.g. involving more than
		// one gene) may require to calculate effects by using all involved genes
		// For some variants we require to be 'N' bases apart (translocations
		// are assumed always involve large genomic regions)
		boolean structuralVariant = variant.isStructural() && (variant.isBnd() || variant.size() > SMALL_VARIANT_SIZE_THRESHOLD);

		// Structural variants?
		if (structuralVariant && !variant.isBnd()) {
			// Large variants could make the query results huge and slow down
			// the algorithm, so we stop here
			// Note: Translocations (BND) only intercept two loci, so this 
			//       issue does not apply.
			if (variant.isStructuralHuge()) {
				variantEffectStructuralLarge(variant, variantEffects);
				return variantEffects;
			}
		}

		//---
		// Query interval tree: Which intervals does variant intersect?
		//---
		Markers intersects = query(variant);

		// In case of large structural variants, we need to check the number of genes
		// involved. If more than one, then we need a different approach (e.g. taking
		// into account all genes involved to calculate fusions)");
		if (structuralVariant) {
			if (variant.isBnd()) {
				System.err.println("Break-ends (rearrangement) variant can not annotated now!");
//				variantEffectBnd(variant, variantEffects, intersects);
				return variantEffects;
			}

			// Calculated effect based on multiple genes
			intersects = variantEffectStructural(variant, variantEffects, intersects);

			// Are we done?
			if (intersects == null) return variantEffects;
		}

		// Calculate variant effect for each query result
		variantEffect(variant, variantEffects, intersects);

		return variantEffects;
	}

	/**
	 * Calculate variant effect for each marker in 'intersect'
	 */
	protected void variantEffect(Variant variant, VariantEffects variantEffects, Markers intersects) {
		// Show all results
		boolean hitChromo = false, hitSomething = false;
		for (Marker marker : intersects) {
			if (marker instanceof Chromosome) hitChromo = true; // Do we hit any chromosome?
			else {
				// Analyze all markers
				if (variant.isNonRef()) {
					marker.variantEffectNonRef(variant, variantEffects);
				}
				else{
					marker.variantEffect(variant, variantEffects);
				}

				// Do we have 'per gene' information?
				if (intervalForestGene != null && marker instanceof Gene) //
					variantEffectGene(marker.getId(), variant, variantEffects);

				hitSomething = true;
			}
		}
		
		// Any errors or intergenic (i.e. did not hit any gene)
		if (!hitChromo) {
			// Special case: Insertion right after chromosome's last base
			Chromosome chr = genome.getChromosome(variant.getChromosomeName());
			if (variant.isIns() && variant.getStart() == (chr.getEnd() + 1)) {
				// This is a chromosome extension
				variantEffects.add(variant, null, EffectType.CHROMOSOME_ELONGATION, "");
			} else if (Config.get().isErrorChromoHit()) {
				variantEffects.addErrorWarning(variant, ErrorWarningType.ERROR_OUT_OF_CHROMOSOME_RANGE);
			}
		} else if (!hitSomething) {
			if (Config.get().isOnlyRegulation()) {
				variantEffects.add(variant, null, EffectType.NONE, "");
			} else {
				variantEffects.add(variant, null, EffectType.INTERGENIC, "");
			}
		}
		
	}

	/**
	 * Add gene specific annotations
	 */
	protected void variantEffectGene(String geneId, Variant variant, VariantEffects variantEffects) {
		Itree itree = intervalForestGene.getTree(geneId);
		if (itree == null) return;

		Markers res = itree.query(variant);
		for (Marker m : res)
			m.variantEffect(variant, variantEffects);
	}

	/**
	 * Calculate structural variant effects taking into account all involved genes
	 *
	 * @return A list of intervals that need to be further analyzed
	 *         or 'null' if no further gene-by-gene analysis is required
	 */
	Markers variantEffectStructural(Variant variant, VariantEffects variantEffects, Markers intersects) {
		// Any variant effects added?
		boolean added = false;

		// Create a new variant effect for structural variants, add effect (if any)
		VariantEffectStructural veff = new VariantEffectStructural(variant, intersects);

		EffectType et = veff.getEffectType();
		boolean considerFussion = (et == EffectType.NONE //
				|| et == EffectType.GENE_FUSION //
				|| et == EffectType.GENE_FUSION_REVERESE //
				|| et == EffectType.GENE_FUSION_HALF //
				|| et == EffectType.FEATURE_FUSION //
		);

		// Note that fusions are added in the next step, when we invoke veff.fusion(), so we skip them here
		if (!considerFussion) {
			variantEffects.add(veff);
			added = true;
		}

		// Do we have a fusion event?
		List<VariantEffect> veffFusions = veff.fusion();
		if (veffFusions != null && !veffFusions.isEmpty()) {
			for (VariantEffect veffFusion : veffFusions) {
				added = true;
				variantEffects.add(veffFusion);
			}
		}

		// In some cases we want to annotate the varaint's partially overlapping genes
		if (variant.isDup()) {
			Markers markers = new Markers();
			for (Marker m : intersects)
				if (!variant.includes(m)) {
					// Note that all these markers overlap the variant so we
					// just filter out the ones fully included in the variant
					markers.add(m);
				}

			// Use these markers for further analysis
			return markers;
		}

		// If variant effects were added, there is no need for further analysis
		return added ? null : intersects;
	}

//	/**
//	 * Calculate translocations variant effects 
//	 */
//	void variantEffectBnd(Variant variant, VariantEffects variantEffects, Markers intersects) {
//		// Create a new variant effect for structural variants, add effect (if any)
//		VariantEffectStructural veff = new VariantEffectStructural(variant, intersects);

		//		EffectType et = veff.getEffectType();
		//		boolean considerFussion = (et == EffectType.NONE //
		//				|| et == EffectType.GENE_FUSION //
		//				|| et == EffectType.GENE_FUSION_REVERESE //
		//				|| et == EffectType.GENE_FUSION_HALF //
		//				|| et == EffectType.FEATURE_FUSION //
		//		);
		//
		//		// Note that fusions are added in the next step, when we invoke veff.fusion(), so we skip them here
		//		if (!considerFussion) {
		//			variantEffects.add(veff);
		//			added = true;
		//		}
		//
		//		// Do we have a fusion event?
		//		List<VariantEffect> veffFusions = veff.fusion();
		//		if (veffFusions != null && !veffFusions.isEmpty()) {
		//			for (VariantEffect veffFusion : veffFusions) {
		//				added = true;
		//				variantEffects.add(veffFusion);
		//			}
		//		}
		//
		//		// In some cases we want to annotate the varaint's partially overlapping genes
		//		if (variant.isDup()) {
		//			Markers markers = new Markers();
		//			for (Marker m : intersects)
		//				if (!variant.includes(m)) {
		//					// Note that all these markers overlap the variant so we
		//					// just filter out the ones fully included in the variant
		//					markers.add(m);
		//				}
		//
		//			// Use these markers for further analysis
		//			return markers;
		//		}
		//
		//		// If variant effects were added, there is no need for further analysis
		//		return added ? null : intersects;
//	}

	/**
	 * Add large structural variant effects
	 */
	void variantEffectStructuralLarge(Variant variant, VariantEffects variantEffects) {
		EffectType eff;

		switch (variant.getVariantType()) {
		case DEL:
			eff = EffectType.CHROMOSOME_LARGE_DELETION;
			break;

		case DUP:
			eff = EffectType.CHROMOSOME_LARGE_DUPLICATION;
			break;

		case INV:
			eff = EffectType.CHROMOSOME_LARGE_INVERSION;
			break;

		default:
			throw new RuntimeException("Unimplemented option for variant type " + variant.getVariantType());
		}

		// Add effect
		variantEffects.add(variant, variant.getChromosome(), eff, "");
	}

}
