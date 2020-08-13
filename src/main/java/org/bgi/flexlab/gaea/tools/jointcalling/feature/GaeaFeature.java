package org.bgi.flexlab.gaea.tools.jointcalling.feature;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;

import htsjdk.tribble.Feature;

public abstract class GaeaFeature implements Feature {
	public GaeaFeature(String name) {
		this.name = name;
	}

	String name;

	protected void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}

	public abstract GenomeLocation getLocation();

	// TODO: this should be a Feature
	public abstract Object getUnderlyingObject();

	/**
	 * wrapping a Tribble feature in a GATK friendly interface
	 */
	public static class TribbleFeature extends GaeaFeature {
		private final GenomeLocationParser genomeLocParser;
		private final Feature feature;
		private GenomeLocation position = null;

		public TribbleFeature(GenomeLocationParser genomeLocParser, Feature f, String name) {
			super(name);
			this.genomeLocParser = genomeLocParser;
			feature = f;
		}

		public GenomeLocation getLocation() {
			if (position == null) {
				position = genomeLocParser.createGenomeLocation(feature.getContig(), feature.getStart(),
						feature.getEnd());
			}
			return position;
		}

		/**
		 * Return the features reference sequence name, e.g chromosome or contig
		 */
		@Override
		public String getChr() {
			return getContig();
		}

		/**
		 * Return the features reference sequence name, e.g chromosome or contig
		 */
		@Override
		public String getContig() {
			return feature.getContig();
		}

		/**
		 * Return the start position in 1-based coordinates (first base is 1)
		 */
		@Override
		public int getStart() {
			return feature.getStart();
		}

		/**
		 * Return the end position following 1-based fully closed conventions.
		 * The length of a feature is end - start + 1;
		 */
		@Override
		public int getEnd() {
			return feature.getEnd();
		}

		// TODO: this should be a Feature, actually
		public Object getUnderlyingObject() {
			return feature;
		}
	}
}
