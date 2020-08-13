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
package org.bgi.flexlab.gaea.tools.annotator.codons;

import org.bgi.flexlab.gaea.tools.annotator.interval.Chromosome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;

import java.util.HashMap;
import java.util.Iterator;

/**
 * All codon tables are stored here. Mapping for genome/chromosome to codon table are also stored here
 *
 * Note: This object is a singleton
 *
 * @author pcingola
 */
public class CodonTables implements Iterable<CodonTable> {

	// Standard codon table
	public static final String STANDARD_TABLE = "TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G";
	public static final String STANDARD_TABLE_NAME = "Standard";
	
	// Vertebrate_Mitochondrial codon table
	public static final String VERTEBRATE_MITOCHONDRIAL_TABLE = "TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/W, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I+, ATC/I+, ATA/M+, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/*, AGG/*, GTT/V, GTC/V, GTA/V, GTG/V+, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G";
	public static final String VERTEBRATE_MITOCHONDRIAL_TABLE_NAME = "Vertebrate_Mitochondrial";

	private static final String KEY_SEPARATOR = "_";
	private static final CodonTables codonTables = new CodonTables();

	HashMap<String, CodonTable> codonTableByName;
	HashMap<String, CodonTable> genChr2codonTable;

	public static CodonTables getInstance() {
		return codonTables;
	}

	private CodonTables() {
		codonTableByName = new HashMap<String, CodonTable>();
		genChr2codonTable = new HashMap<String, CodonTable>();
		CodonTable codonTable = new CodonTable(STANDARD_TABLE_NAME, STANDARD_TABLE);
		add(codonTable);
	}

	/**
	 * Translate a codon into an amino acid for a given genome+chromosome
	 */
	public String aa(String codon, Genome genome, String chromosome) {
		return getTable(genome, chromosome).aa(codon);
	}

	/**
	 * Add a codon table
	 */
	public void add(CodonTable codonTable) {
		codonTableByName.put(codonTable.getName(), codonTable);
	}

	/**
	 * Translate an amino acid into a codon for a given genome+chromosome
	 */
	public String codon(String aa, Genome genome, String chromosome) {
		return getTable(genome, chromosome).codon(aa);
	}

	/**
	 * Get default genome-wide codon table
	 */
	public CodonTable getTable(Genome genome) {
		String key = genome.getId();
		CodonTable codonTable = genChr2codonTable.get(key);
		if (codonTable == null) return codonTables.getTable(STANDARD_TABLE_NAME); // Not found? Use default
		return codonTable;
	}

	/**
	 * Get a codon table
	 * WARNING: It will return the standard codon table if nothing if found
	 */
	public CodonTable getTable(Genome genome, String chromosome) {
		String key = genome.getId() + KEY_SEPARATOR + chromosome;
		CodonTable codonTable = genChr2codonTable.get(key);
		if (codonTable != null) return codonTable;
		return getTable(genome); // Not found? Use genome wide codon table
	}

	/**
	 * Get a codon table by name
	 */
	public CodonTable getTable(String codonTableName) {
		return getInstance().codonTableByName.get(codonTableName);
	}

	@Override
	public Iterator<CodonTable> iterator() {
		return codonTableByName.values().iterator();
	}

	/**
	 * Set a codon table for a given genome & chromosome
	 */
	public void set(Genome genome, Chromosome chr, CodonTable codonTable) {
		add(codonTable); // Just in case it's not already added
		String key = genome.getId() + KEY_SEPARATOR + chr.getId();
		genChr2codonTable.put(key, codonTable);
	}

	/**
	 * Set a codon table for a all chromosomes in a genome
	 * I.e.: Default genome-wide chromosome table
	 */
	public void set(Genome genome, CodonTable codonTable) {
		add(codonTable); // Just in case it's not already added
		String key = genome.getId();
		genChr2codonTable.put(key, codonTable);
	}

}
