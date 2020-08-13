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
import org.bgi.flexlab.gaea.tools.annotator.util.Timer;

import java.io.Serializable;
import java.util.*;

/**
 * A collection of genes (marker intervals)
 * Note: It is assumed that all genes belong to the same genome
 *
 * @author pcingola
 */
public class Genes implements Iterable<Gene>, Serializable {

	private static final long serialVersionUID = 9022385501946879197L;
	public static final String CIRCULAR_GENE_ID = "_circ";

	public boolean debug = false;
	Genome genome;
	HashMap<String, Gene> genesById;

	public Genes(Genome genome) {
		genesById = new HashMap<String, Gene>();
		this.genome = genome;
	}

	/**
	 * Add a gene interval to this collection
	 */
	public void add(Gene gene) {
		genesById.put(gene.getId(), gene);
	}

	/** In a circular genome, a gene can have negative coordinates or crosses
		over chromosome end. These genes are mirrored to the opposite end of
		the chromosome so that they can be referenced by both circular coordinates.
	 */
	public void createCircularGenes() {
		List<Gene> newGenes = new LinkedList<Gene>();

		// Check if any gene spans across chromosome limits
		for (Gene g : genome.getGenes()) {
			Chromosome chr = g.getChromosome();

			Gene newGene = null;

			if ((g.getStart() < 0) || (g.getEnd() > chr.getEnd())) {
				newGene = (Gene) g.clone();

				// Change IDs
				newGene.setId(g.getId() + CIRCULAR_GENE_ID);
				for (Transcript tr : newGene) {
					tr.setId(tr.getId() + CIRCULAR_GENE_ID);
					for (Exon ex : tr)
						ex.setId(ex.getId() + CIRCULAR_GENE_ID);
				}

				// Shift coordinates
				int shift = 0;
				if (g.getStart() < 0) {
					shift = chr.size();
				} else if (g.getEnd() > chr.getEnd()) {
					shift = -chr.size();
				}

				newGene.shiftCoordinates(shift);
				if (Config.get().isVerbose()) Timer.showStdErr("Gene '" + g.getId() + "' spans across coordinate zero: Assuming circular chromosome, creating mirror gene at the end." //
						+ "\n\tGene        :" + g.toStr() //
						+ "\n\tNew gene    :" + newGene.toStr() //
						+ "\n\tChrsomosome :" + chr.toStr() //
				);

				// Add them to genes
				newGenes.add(newGene);
			}
		}

		// Add all newly created genes
		if (!newGenes.isEmpty()) {
			for (Gene g : newGenes)
				genome.getGenes().add(g);
			if (Config.get().isVerbose()) Timer.showStdErr("Total: " + newGenes.size() + " added as circular mirrored genes (appended '" + CIRCULAR_GENE_ID + "' to IDs).");
		}

	}

	/**
	 * Creates a list of Intergenic regions
	 */
	public List<Intergenic> createIntergenic() {
		ArrayList<Intergenic> intergenics = new ArrayList<Intergenic>(genesById.size());

		// Create a list of genes sorted by position
		ArrayList<Gene> genesSorted = new ArrayList<Gene>(genesById.size());
		genesSorted.addAll(genesById.values());
		Collections.sort(genesSorted);

		// For each gene, transcript
		Gene genePrev = null;
		Chromosome chrPrev = null;
		for (Gene gene : genesSorted) {
			// Chromosome change? Invaludate genePrev
			if (chrPrev != gene.getChromosome()) {

				// Add last intergenic region in the chromosome
				if (chrPrev != null && genePrev != null) {
					int start = genePrev.getEnd() + 1;
					int end = chrPrev.getEnd();
					if (start < end) {
						String id = genePrev.getId() + "-END";
						String name = genePrev.getGeneName() + "-END";
						Intergenic intergenic = new Intergenic(genePrev.getChromosome(), start, end, false, id, name);
						intergenics.add(intergenic);
					}
				}

				genePrev = null;
			}

			// Intergenic region's [start, end] interval
			int start = (genePrev != null ? genePrev.getEnd() + 1 : 0);
			int end = gene.getStart() - 1;

			// Valid intergenic region?
			if (start < end) {
				String name = (genePrev != null ? genePrev.getGeneName() + "-" : "") + gene.getGeneName();
				String id = (genePrev != null ? genePrev.getId() + "-" : "") + gene.getId();
				Intergenic intergenic = new Intergenic(gene.getChromosome(), start, end, false, id, name);
				intergenics.add(intergenic);
			}

			// Is it null or ends before this one? update 'genePrev'
			if ((genePrev == null) || (gene.getEnd() > genePrev.getEnd())) genePrev = gene;

			// Update chrPrev
			chrPrev = gene.getChromosome();
		}

		return intergenics;
	}

	/**
	 * Create splice sites.
	 *
	 * @param createIfMissing : If true, create canonical splice sites if they are missing.
	 *
	 * For a definition of splice site, see comments at the beginning of SpliceSite.java
	 */
	public void createSpliceSites(int spliceSiteSize, int spliceRegionExonSize, int spliceRegionIntronMin, int spliceRegionIntronMax) {
		// For each gene, transcript
		for (Gene gene : this)
			for (Transcript tr : gene)
				tr.createSpliceSites(spliceSiteSize, spliceRegionExonSize, spliceRegionIntronMin, spliceRegionIntronMax);
	}

	/**
	 * Creates a list of UP/DOWN stream regions (for each transcript)
	 * Upstream (downstream) stream is defined as upDownLength before (after) transcript
	 *
	 * Note: If upDownLength <=0 no interval is created
	 */
	public List<Marker> createUpDownStream(int upDownLength) {
		ArrayList<Marker> list = new ArrayList<Marker>();
		if (upDownLength <= 0) return list;

		// For each gene, transcript
		for (Gene gene : this) {
			for (Transcript tr : gene) {
				tr.createUpDownStream(upDownLength);
				list.add(tr.getUpstream());
				list.add(tr.getDownstream());
			}
		}
		return list;
	}

	/**
	 * Find a transcript by ID
	 */
	public Transcript findTranscript(String trId) {
		for (Gene g : this)
			for (Transcript tr : g)
				if (tr.getId().equals(trId)) return tr;
		return null;
	}

	/**
	 * Obtain a gene interval
	 */
	public Gene get(String geneId) {
		return genesById.get(geneId);
	}

	@Override
	public Iterator<Gene> iterator() {
		return genesById.values().iterator();
	}

	public int size() {
		return genesById.size();
	}

	public Collection<Gene> sorted() {
		ArrayList<Gene> genes = new ArrayList<Gene>();
		genes.addAll(genesById.values());
		Collections.sort(genes, new IntervalComparatorByStart());
		return genes;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Gene gint : this)
			sb.append(gint + "\n");
		return sb.toString();
	}

	public Collection<Gene> values() {
		return genesById.values();
	}

}
