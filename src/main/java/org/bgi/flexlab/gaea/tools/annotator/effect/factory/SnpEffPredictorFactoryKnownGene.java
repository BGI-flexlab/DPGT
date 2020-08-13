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
package org.bgi.flexlab.gaea.tools.annotator.effect.factory;

import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.effect.SnpEffectPredictor;
import org.bgi.flexlab.gaea.tools.annotator.interval.*;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;
import org.bgi.flexlab.gaea.tools.annotator.util.MultivalueHashMap;

import java.io.BufferedReader;
import java.util.List;

/**
 * This class creates a SnpEffectPredictor from a TXT file dumped using UCSC table browser
 *
 * Fields in this table
 *
 * Field        Example               SQL type            Info     Description
 * -----        -------               --------            ----     -----------
 * name         uc001aaa.3            varchar(255)        values   Name of gene
 * chrom        chr1                  varchar(255)        values   Reference sequence chromosome or scaffold
 * strand       +                     char(1)             values   + or - for strand
 * txStart      11873                 int(10) unsigned    range    Transcription start position
 * txEnd        14409                 int(10) unsigned    range    Transcription end position
 * cdsStart     11873                 int(10) unsigned    range    Coding region start
 * cdsEnd       11873                 int(10) unsigned    range    Coding region end
 * exonCount    3                     int(10) unsigned    range    Number of exons
 * exonStarts   11873,12612,13220,    longblob                     Exon start positions
 * exonEnds     12227,12721,14409,    longblob                     Exon end positions
 * proteinID                          varchar(40)         values   UniProt display ID for Known Genes, UniProt accession or RefSeq protein ID for UCSC Genes
 * alignID      uc001aaa.3            varchar(255)        values   Unique identifier for each (known gene, alignment position) pair
 *
 * @author pcingola
 */
public class SnpEffPredictorFactoryKnownGene extends SnpEffPredictorFactory {

	public static final String CDS_STAT_COMPLETE = "cmpl";

	int ignoredTr = 0;
	MultivalueHashMap<String, Gene> genesByName;

	public SnpEffPredictorFactoryKnownGene(Config config) {
		super(config, 0); // Zero values coordinates

		genesByName = new MultivalueHashMap<String, Gene>();

		frameType = FrameType.UCSC;
		frameCorrection = true;
	}

	@Override
	public SnpEffectPredictor create() {
		try {
			// Read gene intervals from a file
			if (fileName == null) fileName = Config.get().getGeneInfo();

			System.out.println("Reading gene intervals file : '" + fileName + "'");
			readRefSeqFile(); // Read gene info

			beforeExonSequences(); // Some clean-up before readng exon sequences

			// Read chromosome sequences and set exon sequences
			if (readSequences) readExonSequences();
			else if (createRandSequences) createRandSequences();

			finishUp(); // Perform adjustments

			if (verbose) {
				System.out.println(config.getGenome());
				System.out.println("# Ignored transcripts        : " + ignoredTr);
			}
		} catch (Exception e) {
			if (verbose) e.printStackTrace();
			throw new RuntimeException("Error reading file '" + fileName + "'\n" + e);
		}
		return snpEffectPredictor;
	}

	/**
	 * Find (or create) a gene for this transcript
	 */
	Gene findOrCreateGene(String geneName, String trId, Chromosome chromo, int start, int end, boolean strandMinus, boolean isCoding) {
		Marker tr = new Marker(chromo, start, end, strandMinus, trId);
		List<Gene> genes = genesByName.get(geneName);
		int geneIndex = 0;
		if (genes != null) {
			for (Gene gene : genes) {
				if (gene.intersects(tr)) {
					// Do we need to update gene length?
					if (start < gene.getStart()) gene.setStart(start);
					if (gene.getEnd() < end) gene.setEnd(end);

					return gene;
				}
			}

			geneIndex = genes.size() + 1;
		}

		// Need to create a new gene
		String geneId = geneName + (geneIndex > 0 ? "." + geneIndex : "");
		Gene gene = new Gene(chromo, start, end, strandMinus, geneId, geneName, BioType.coding(isCoding));
		genesByName.add(geneName, gene);
		add(gene);

		return gene;
	}

	/**
	 * Read and parse genes file
	 */
	protected void readRefSeqFile() {
		try {
			int count = 0;
			BufferedReader reader = Gpr.reader(fileName);
			if (reader == null) return; // Error

			for (lineNum = 1; reader.ready(); lineNum++) {
				line = reader.readLine();

				// Skip headers
				if (!line.startsWith("#")) {
					String fields[] = line.split("\t");

					if (fields.length >= 9) {
						// Parse fields
						int fieldNum = 0;
						String id = fields[fieldNum++];
						String chromoName = fields[fieldNum++];
						boolean strandMinus = fields[fieldNum++].equals("-");

						int txstart = parsePosition(fields[fieldNum++]);
						int txend = parsePosition(fields[fieldNum++]) - 1; // Our internal database representations of coordinates always have a zero-based start and a one-based end (Reference: http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1 )

						int cdsStart = parsePosition(fields[fieldNum++]);
						int cdsEnd = parsePosition(fields[fieldNum++]) - 1; // Our internal database representations of coordinates always have a zero-based start and a one-based end (Reference: http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1 )

						int exonCount = Gpr.parseIntSafe(fields[fieldNum++]);
						String exonStarts = fields[fieldNum++];
						String exonEnds = fields[fieldNum++]; // Our internal database representations of coordinates always have a zero-based start and a one-based end (Reference: http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1 )

						String proteinId = fields[fieldNum++];
						// String alignId = fields[fieldNum++]; // Not used

						//---
						// Create
						//----
						Chromosome chromo = getOrCreateChromosome(chromoName);

						// Is it protein coding?
						boolean isCoding = !proteinId.isEmpty(); // Protein ID assigned?

						// Create IDs
						String trId = uniqueTrId(id);

						// Get or create gene
						Gene gene = findOrCreateGene(proteinId, trId, chromo, txstart, txend, strandMinus, isCoding);

						// Create transcript
						Transcript tr = new Transcript(gene, txstart, txend, strandMinus, trId);
						tr.setProteinCoding(isCoding);
						add(tr);

						// Add Exons and CDS
						String exStartStr[] = exonStarts.split(",");
						String exEndStr[] = exonEnds.split(",");
						for (int i = 0; i < exonCount; i++) {
							// Exons
							int exStart = parsePosition(exStartStr[i]);
							int exEnd = parsePosition(exEndStr[i]) - 1; // Our internal database representations of coordinates always have a zero-based start and a one-based end (Reference: http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1 )
							String exId = trId + ".ex." + (i + 1);
							Exon ex = new Exon(tr, exStart, exEnd, strandMinus, exId, i);
							add(ex);

							// CDS (ony if intersects)
							if ((exStart <= cdsEnd) && (exEnd >= cdsStart)) {
								Cds cds = new Cds(tr, Math.max(cdsStart, exStart), Math.min(cdsEnd, exEnd), strandMinus, exId);
								add(cds);
							}
						}

						count++;
						if (count % MARK == 0) System.out.print('.');
						if (count % (100 * MARK) == 0) System.out.print("\n\t");
					}
				}
			}

			reader.close();
		} catch (Exception e) {
			Gpr.debug("Offending line (lineNum: " + lineNum + "): '" + line + "'");
			throw new RuntimeException(e);
		}
	}

	/**
	 * Create a new (unique) transcript ID
	 * @param id
	 * @return
	 */
	String uniqueTrId(String id) {
		if (!transcriptsById.containsKey(id)) return id;
		for (int i = 2; true; i++) {
			String trId = id + "." + i;
			if (!transcriptsById.containsKey(trId)) return trId;
		}
	}
}
