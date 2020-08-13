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
import org.bgi.flexlab.gaea.tools.annotator.effect.SnpEffectPredictor;
import org.bgi.flexlab.gaea.tools.annotator.interval.Exon.ExonSpliceType;
import org.bgi.flexlab.gaea.tools.annotator.util.CountByType;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;
import org.bgi.flexlab.gaea.tools.annotator.util.Timer;

import java.util.HashMap;

/**
 * Characterize exons based on alternative splicing
 *
 * References: "Alternative splicing and evolution - diversification, exon definition and function"  (see Box 1)
 *
 * @author pablocingolani
 */
public class ExonSpliceCharacterizer {

	public static final int MAX_EXONS = 1000; // Do not characterize transcripts having more than this number of exons
	public static final int SHOW_EVERY = 1000;

	boolean verbose = false;
	Genome genome;
	HashMap<Exon, Exon.ExonSpliceType> typeByExon;
	CountByType countByType = new CountByType();

	public ExonSpliceCharacterizer(Genome genome) {
		this.genome = genome;
		typeByExon = new HashMap<Exon, Exon.ExonSpliceType>();
	} 

	public ExonSpliceCharacterizer(String genomeVer) {
//		Config config = new Config(genomeVer); 
//		TODO 处理
		Config config = Config.getConfigInstance();
		SnpEffectPredictor snpEffectPredictor = config.getSnpEffectPredictor();
		genome = snpEffectPredictor.getGenome();
		typeByExon = new HashMap<Exon, Exon.ExonSpliceType>();
	}

	/**
	 * Characterize all exons
	 */
	public CountByType characterize() {
		type();
		return countByType;
	}

	/**
	 * Count number of exons
	 */
	int countExons() {
		int count = 0;
		for (Gene g : genome.getGenes())
			for (Transcript tr : g)
				count += tr.numChilds();
		return count;
	}

	/**
	 * Does the marker intersect any exon in 'tr'?
	 * @param m
	 * @param tr
	 * @return
	 */
	boolean intersectsAnyExon(Marker m, Transcript tr) {
		for (Exon e : tr)
			if (m.intersects(e)) return true;
		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_3SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt3ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon end, different exon start?
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					} else {
						// Same exon end, different exon start? (negative strand)
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is thins an ALTERNATIVE_5SS exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isAlt5ss(Exon exon, Gene gene) {
		for (Transcript tr : gene)
			for (Exon e : tr) {
				if (exon.intersects(e)) {
					if (exon.isStrandPlus()) {
						// Same exon start, different exon end?
						if ((exon.getStart() == e.getStart()) && (exon.getEnd() != e.getEnd())) return true;
					} else {
						// Same exon start, different exon end? (negative strand)
						if ((exon.getStart() != e.getStart()) && (exon.getEnd() == e.getEnd())) return true;
					}
				}
			}

		return false;
	}

	/**
	 * Is this exon mutually exclusive with another exon?
	 * @param exon
	 * @param gene
	 * @return
	 */
	boolean isMutEx(Exon exon, Gene gene) {
		if (gene.numChilds() <= 1) return false;

		//---
		// Make a list of all unique exons
		//---
		String exonKey = key(exon);
		HashMap<String, Exon> uniqEx = new HashMap<String, Exon>();
		for (Transcript tr : gene)
			for (Exon e : tr) {
				String ekey = key(e);
				if (!exonKey.equals(ekey)) uniqEx.put(ekey, e);
			}

		//---
		// For each unique exon, compare if it is mutually exclusive with 'exon'
		//---
		Transcript exonTr = (Transcript) exon.getParent();
		for (Exon e : uniqEx.values()) {
			ExonSpliceType type = typeByExon.get(e);;
			if (type == Exon.ExonSpliceType.SKIPPED) { // Only check for these exons (avoid all ALT*)
				boolean xor = true;
				for (Transcript tr : gene) {
					if (exonTr.intersects(tr) && exon.intersects(tr) && e.intersects(exonTr)) { // Make sure both exons intersect both transcripts (otherwise they cannot be mutually exclusive)
						xor &= intersectsAnyExon(e, tr) ^ intersectsAnyExon(exon, tr);
					} else xor = false;
				}

				// XOR is true? => Mutually exclusive
				if (xor) return true;
			}
		}

		return false;
	}

	/**
	 * Create a simple hash key based on choromosomal position
	 * @param m
	 * @return
	 */
	String key(Marker m) {
		return m.getChromosomeName() + ":" + m.getStart() + "-" + m.getEnd();
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Mark exons types
	 */
	void type() {
		if (verbose) {
			Timer.showStdErr("Caracterizing exons by splicing (stage 1) : ");
			System.out.print("\t");
		}

		// Find retained exons
		int numExon = 1;
		for (Gene g : genome.getGenes()) {

			// Count exons
			CountByType count = new CountByType();
			for (Transcript tr : g)
				for (Exon e : tr)
					count.inc(key(e));

			// Label exons
			int countTr = g.numChilds();
			for (Transcript tr : g) {
				for (Exon e : tr) {
					if (verbose) Gpr.showMark(numExon++, SHOW_EVERY, "\t");

					String eKey = key(e);
					int countEx = (int) count.get(eKey);

					// Is this exon maintained in all transcripts?
					if (countEx == countTr) type(e, Exon.ExonSpliceType.RETAINED);
					else {
						if (isAlt3ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_3SS);
						else if (isAlt5ss(e, g)) type(e, Exon.ExonSpliceType.ALTTENATIVE_5SS);
						else if (tr.numChilds() > 1) {
							if (e.getRank() == 1) type(e, Exon.ExonSpliceType.ALTTENATIVE_PROMOMOTER);
							else if (e.getRank() == tr.numChilds()) type(e, Exon.ExonSpliceType.ALTTENATIVE_POLY_A);
							else type(e, Exon.ExonSpliceType.SKIPPED);
						}
					}
				}
			}
		}

		if (verbose) {
			System.err.println("");
			Timer.showStdErr("Caracterizing exons by splicing (stage 2) : ");
			System.out.print("\t");
		}

		// Now analyze if there are mutually exclusive exons
		numExon = 1;
		for (Gene g : genome.getGenes()) {
			for (Transcript tr : g) {
				if (tr.numChilds() < MAX_EXONS) {
					for (Exon e : tr) {
						if (verbose) Gpr.showMark(numExon++, SHOW_EVERY, "\t");
						ExonSpliceType type = typeByExon.get(e);
						if (type == ExonSpliceType.SKIPPED) { // Try to re-annotate only these
							if (isMutEx(e, g)) type(e, Exon.ExonSpliceType.MUTUALLY_EXCLUSIVE);
						}
					}
				} else {
					System.err.println("");
					Gpr.debug("WARNING: Gene '" + g.getId() + "', transcript '" + tr.getId() + "' has too many exons (" + tr.numChilds() + " exons). Skipped");
				}
			}
		}

		if (verbose) Timer.showStdErr("done.");
	}

	/**
	 * Mark this exons as 'type'
	 * @param e
	 * @param type
	 */
	void type(Exon e, Exon.ExonSpliceType type) {
		e.spliceType = type;
		countByType.inc(type.toString());
		typeByExon.put(e, type);
	}
}
