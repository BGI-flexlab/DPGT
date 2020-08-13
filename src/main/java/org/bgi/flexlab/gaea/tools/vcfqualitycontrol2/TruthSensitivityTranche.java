package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.exception.MalformedFile;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.XReadLines;

final class TruthSensitivityTranche extends Tranche {
	private static final int CURRENT_VERSION = 5;

	static final Comparator<TruthSensitivityTranche> TRUTH_SENSITIVITY_ORDER = (tranche1, tranche2) -> Double
			.compare(tranche1.targetTruthSensitivity, tranche2.targetTruthSensitivity);

	// Note: visibility is set to package-local for testing
	final double targetTruthSensitivity;

	public TruthSensitivityTranche(final double targetTruthSensitivity, final double minVQSLod, final int numKnown,
			final double knownTiTv, final int numNovel, final double novelTiTv, final int accessibleTruthSites,
			final int callsAtTruthSites, final VCFQualityControlArgumentCollection.Mode model) {
		this(targetTruthSensitivity, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites,
				callsAtTruthSites, model, "anonymous");
	}

	public TruthSensitivityTranche(final double targetTruthSensitivity, final double minVQSLod, final int numKnown,
			final double knownTiTv, final int numNovel, final double novelTiTv, final int accessibleTruthSites,
			final int callsAtTruthSites, final VCFQualityControlArgumentCollection.Mode model, final String name) {
		super(name, knownTiTv, numNovel, minVQSLod, model, novelTiTv, accessibleTruthSites, numKnown,
				callsAtTruthSites);
		if (targetTruthSensitivity < 0.0 || targetTruthSensitivity > 100.0) {
			throw new UserException("Target FDR is unreasonable " + targetTruthSensitivity);
		}

		if (numKnown < 0 || numNovel < 0) {
			throw new UserException("Invalid tranche - no. variants is < 0 : known " + numKnown + " novel " + numNovel);
		}

		if (name == null) {
			throw new UserException("BUG -- name cannot be null");
		}

		this.targetTruthSensitivity = targetTruthSensitivity;

	}

	public Double getTrancheIndex() {
		return targetTruthSensitivity;
	}

	public static String printHeader() {
		try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
				final PrintStream stream = new PrintStream(bytes)) {

			stream.println("# Variant quality score tranches file");
			stream.println("# Version number " + CURRENT_VERSION);
			stream.println(
					"targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

			return bytes.toString();
		} catch (IOException e) {
			throw new UserException("IOException while converting tranche to a string");
		}
	}

	@Override
	public String toString() {
		return String.format(
				"TruthSensitivityTranche targetTruthSensitivity=%.2f minVQSLod=%.4f known=(%d @ %.4f) novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
				targetTruthSensitivity, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites,
				callsAtTruthSites, name);
	}

	/**
	 * Returns a list of tranches, sorted from most to least specific, read in
	 * from file f.
	 * 
	 * @throws IOException
	 *             if there are problems reading the file.
	 */
	public static List<TruthSensitivityTranche> readTranches(final File f) throws IOException {
		String[] header = null;
		List<TruthSensitivityTranche> tranches = new ArrayList<>();

		try (XReadLines xrl = new XReadLines(f)) {
			for (final String line : xrl) {
				if (line.startsWith(COMMENT_STRING)) {
					continue;
				}

				final String[] vals = line.split(VALUE_SEPARATOR);
				if (header == null) { // reading the header
					header = vals;
					if (header.length != EXPECTED_COLUMN_COUNT) {
						throw new MalformedFile(f, "Expected 11 elements in header line " + line);
					}
				} else {
					if (header.length != vals.length) {
						throw new MalformedFile(f, "Line had too few/many fields.  Header = " + header.length + " vals "
								+ vals.length + ". The line was: " + line);
					}

					Map<String, String> bindings = new LinkedHashMap<>();
					for (int i = 0; i < vals.length; i++) {
						bindings.put(header[i], vals[i]);
					}
					tranches.add(new TruthSensitivityTranche(getRequiredDouble(bindings, "targetTruthSensitivity"),
							getRequiredDouble(bindings, "minVQSLod"), getOptionalInteger(bindings, "numKnown", -1),
							getOptionalDouble(bindings, "knownTiTv", -1.0), getRequiredInteger(bindings, "numNovel"),
							getRequiredDouble(bindings, "novelTiTv"),
							getOptionalInteger(bindings, "accessibleTruthSites", -1),
							getOptionalInteger(bindings, "callsAtTruthSites", -1),
							VCFQualityControlArgumentCollection.Mode.valueOf(bindings.get("model")),
							bindings.get("filterName")));
				}
			}
		}

		tranches.sort(TRUTH_SENSITIVITY_ORDER);
		return tranches;
	}

	static final class TruthSensitivityMetric {
		private final String name;
		private double[] runningSensitivity;
		private final int nTrueSites;

		public TruthSensitivityMetric(final int nTrueSites) {
			this.name = "TruthSensitivity";
			this.nTrueSites = nTrueSites;
		}

		public String getName() {
			return name;
		}

		public double getThreshold(final double tranche) {
			return 1.0 - tranche / 100.0; // tranche of 1 => 99% sensitivity
											// target
		}

		public void calculateRunningMetric(final List<VariantDatum> data) {
			int nCalledAtTruth = 0;
			runningSensitivity = new double[data.size()];

			for (int i = data.size() - 1; i >= 0; i--) {
				VariantDatum datum = data.get(i);
				nCalledAtTruth += datum.atTruthSite ? 1 : 0;
				runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
			}
		}

		public double getRunningMetric(final int i) {
			return runningSensitivity[i];
		}
	}

	public static List<TruthSensitivityTranche> findTranches(final List<VariantDatum> data,
			final double[] trancheThresholds, final TruthSensitivityMetric metric,
			final VCFQualityControlArgumentCollection.Mode model) {

		Collections.sort(data, VariantDatum.VariantDatumLODComparator);
		metric.calculateRunningMetric(data);

		List<TruthSensitivityTranche> tranches = new ArrayList<>();
		for (double trancheThreshold : trancheThresholds) {
			TruthSensitivityTranche t = findTranche(data, metric, trancheThreshold, model);

			if (t == null) {
				if (tranches.isEmpty()) {
					throw new UserException(String.format(
							"Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?",
							metric.getName(), metric.getThreshold(trancheThreshold)));
				}
				break;
			}

			tranches.add(t);
		}

		tranches.sort(TruthSensitivityTranche.TRUTH_SENSITIVITY_ORDER);
		return tranches;
	}

	private static TruthSensitivityTranche findTranche(final List<VariantDatum> data,
			final TruthSensitivityMetric metric, final double trancheThreshold,
			final VCFQualityControlArgumentCollection.Mode model) {

		double metricThreshold = metric.getThreshold(trancheThreshold);
		int n = data.size();
		for (int i = 0; i < n; i++) {
			if (metric.getRunningMetric(i) >= metricThreshold) {
				// we've found the largest group of variants with sensitivity >=
				// our target truth sensitivity
				TruthSensitivityTranche t = trancheOfVariants(data, i, trancheThreshold, model);
				return t;
			}
		}

		return null;
	}

	protected static TruthSensitivityTranche trancheOfVariants(final List<VariantDatum> data, final int minI,
			final double ts, final VCFQualityControlArgumentCollection.Mode model) {
		Tranche basicTranche = Tranche.trancheOfVariants(data, minI, ts, model);
		return new TruthSensitivityTranche(ts, basicTranche.minVQSLod, basicTranche.numKnown, basicTranche.knownTiTv,
				basicTranche.numNovel, basicTranche.novelTiTv, basicTranche.accessibleTruthSites,
				basicTranche.callsAtTruthSites, model, DEFAULT_TRANCHE_NAME);
	}

}
