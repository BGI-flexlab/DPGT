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
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;

import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.effect.SnpEffectPredictor;
import org.bgi.flexlab.gaea.tools.annotator.effect.factory.SnpEffPredictorFactory;
import org.bgi.flexlab.gaea.tools.annotator.effect.factory.SnpEffPredictorFactoryRefSeq;
import org.bgi.flexlab.gaea.tools.annotator.interval.Gene;
import org.bgi.flexlab.gaea.tools.annotator.interval.SpliceSite;
import org.bgi.flexlab.gaea.tools.annotator.interval.Transcript;
import org.bgi.flexlab.gaea.tools.annotator.interval.TranscriptSupportLevel;
import org.bgi.flexlab.gaea.tools.annotator.util.Timer;

import java.io.Serializable;

class AnnotatorBuild implements Serializable{
	
	private static final long serialVersionUID = 8558515853505312687L;
	
	private boolean debug; // Debug mode
	private boolean verbose; // Be verbose
	private boolean canonical = false; // Use only canonical transcripts
	private boolean strict = false; // Only use transcript that have been validated
	
	private int spliceSiteSize = SpliceSite.CORE_SPLICE_SITE_SIZE; // Splice site size default: 2 bases (canonical splice site)
	private int spliceRegionExonSize = SpliceSite.SPLICE_REGION_EXON_SIZE;
	private int spliceRegionIntronMin = SpliceSite.SPLICE_REGION_INTRON_MIN;
	private int spliceRegionIntronMax = SpliceSite.SPLICE_REGION_INTRON_MAX;
	private int upDownStreamLength = SnpEffectPredictor.DEFAULT_UP_DOWN_LENGTH; // Upstream & downstream interval length
	
	private TranscriptSupportLevel maxTranscriptSupportLevel = null; // Filter by maximum Transcript Support Level (TSL)
	private boolean onlyProtein = false; // Only use protein coding transcripts
	
	private Config config;
	boolean storeSequences = false; // Store full sequences
	
	AnnotatorBuild(Config config) {
		this.config = config;
		this.debug = config.isDebug();
		this.verbose = config.isVerbose();
	}
	
	/**
	 * Create SnpEffectPredictor
	 */
	SnpEffectPredictor createSnpEffPredictor() {

		// Create factory
		SnpEffPredictorFactory factory;
		
//		TODO 支持多种基因信息格式:  RefSeq, EMBL, UCSC KnownGenes ...
		factory = new SnpEffPredictorFactoryRefSeq(config);

		// Create SnpEffPredictors
		factory.setVerbose(verbose);
		factory.setDebug(debug);
		factory.setStoreSequences(storeSequences);
		
		return factory.create();
	}
	
	void buildForest(){
		// Set upstream-downstream interval length
		config.getSnpEffectPredictor().setUpDownStreamLength(upDownStreamLength);

		// Set splice site/region sizes
		config.getSnpEffectPredictor().setSpliceSiteSize(spliceSiteSize);
		config.getSnpEffectPredictor().setSpliceRegionExonSize(spliceRegionExonSize);
		config.getSnpEffectPredictor().setSpliceRegionIntronMin(spliceRegionIntronMin);
		config.getSnpEffectPredictor().setSpliceRegionIntronMax(spliceRegionIntronMax);

		// Filter canonical transcripts
		if (canonical) {
			if (verbose) Timer.showStdErr("Filtering out non-canonical transcripts.");
			config.getSnpEffectPredictor().removeNonCanonical();

			if (verbose) {
				// Show genes and transcript (which ones are considered 'canonical')
				Timer.showStdErr("Canonical transcripts:\n\t\tgeneName\tgeneId\ttranscriptId\tcdsLength");
				for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes()) {
					for (Transcript t : g) {
						String cds = t.cds();
						int cdsLen = (cds != null ? cds.length() : 0);
						System.err.println("\t\t" + g.getGeneName() + "\t" + g.getId() + "\t" + t.getId() + "\t" + cdsLen);
					}
				}
			}
			if (verbose) Timer.showStdErr("done.");
		}

		// Filter transcripts by TSL
		if (maxTranscriptSupportLevel != null) {
			if (verbose) Timer.showStdErr("Filtering transcripts by Transcript Support Level (TSL): " + maxTranscriptSupportLevel);
			config.getSnpEffectPredictor().filterTranscriptSupportLevel(maxTranscriptSupportLevel);

			if (verbose) {
				// Show genes and transcript (which ones are considered 'canonical')
				Timer.showStdErr("Transcript:\n\t\tgeneName\tgeneId\ttranscriptId\tTSL");
				for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes()) {
					for (Transcript t : g)
						System.err.println("\t\t" + g.getGeneName() + "\t" + g.getId() + "\t" + t.getId() + "\t" + t.getTranscriptSupportLevel());
				}
			}
			if (verbose) Timer.showStdErr("done.");
		}

		// Filter verified transcripts
		if (strict) {
			if (verbose) Timer.showStdErr("Filtering out non-verified transcripts.");
			if (config.getSnpEffectPredictor().removeUnverified()) {
				fatalError("All transcripts have been removed form every single gene!\nUsing strickt on this database leaves no information.");
			}
			if (verbose) Timer.showStdErr("done.");
		}

		// Use transcripts set form input file
//		if (onlyTranscriptsFile != null) {
//			// Load file
//			String onlyTr = Gpr.readFile(onlyTranscriptsFile);
//			HashSet<String> trIds = new HashSet<String>();
//			for (String trId : onlyTr.split("\n"))
//				trIds.add(trId.trim());
//
//			// Remove transcripts
//			if (verbose) Timer.showStdErr("Filtering out transcripts in file '" + onlyTranscriptsFile + "'. Total " + trIds.size() + " transcript IDs.");
//			int removed = config.getSnpEffectPredictor().retainAllTranscripts(trIds);
//			if (verbose) Timer.showStdErr("Done: " + removed + " transcripts removed.");
//		}

		// Use protein coding transcripts
		if (onlyProtein) {
			// Remove transcripts
			if (verbose) Timer.showStdErr("Filtering out non-protein coding transcripts.");
			int removed = config.getSnpEffectPredictor().retainTranscriptsProtein();
			if (verbose) Timer.showStdErr("Done: " + removed + " transcripts removed.");
		}

		// Build tree
		if (verbose) Timer.showStdErr("Building interval forest");
		config.getSnpEffectPredictor().buildForest();
		if (verbose) Timer.showStdErr("done.");
	}
	
	/**
	 * Show an error message and exit
	 */
	public void fatalError(String message) {
		System.err.println("Fatal error: " + message);
		System.exit(-1);
	}
	
	
	
}
