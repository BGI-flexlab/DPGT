package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.TrainData;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.mode.GaussianMixtureModel;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.mode.MultivariateGaussian;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.report.Report;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.report.ReportTable;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.RScriptUtils;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.VCFUtils;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.QualityUtils;
import org.bgi.flexlab.gaea.util.Utils;

import Jama.Matrix;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class VariantRecalibrator {
	private VariantDataManager dataManager;
	private VCFQualityControlOptions options;
	private GenomeLocationParser parser = null;

	private final List<VariantContext> variantsAtLocus = new ArrayList<>();

	private int counter = 0;

	final private ExpandingArrayList<VariantDatum> reduceSum = new ExpandingArrayList<>(2000);

	private final VariantRecalibratorEngine engine;

	private PrintStream tranchesStream = null;

	private VariantContextWriter recalWriter = null;

	private VCFHeader header = null;

	public static  boolean lenientVCFProcessing = false;

	public static  boolean createOutputVariantIndex = true;

	public static boolean createOutputVariantMD5 = false;

	private List<TruthSensitivityTranche> tranches = null;
	
	public VariantRecalibrator(VCFQualityControlOptions options,VCFHeader header){
		this.options = options;
		this.header = header;
		
		if(this.parser == null){
			this.parser = new GenomeLocationParser(header.getSequenceDictionary());
		}
		dataManager = new VariantDataManager(options.getUseAnnotations(), options.getArgument());

		engine = new VariantRecalibratorEngine(options.getArgument());

		try {
			tranchesStream = new PrintStream(options.getTrancheFile());
		} catch (FileNotFoundException e) {
			throw new UserException.CouldNotCreateOutputFile(options.getTrancheFile(), e);
		}

		Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
		VCFUtils.addVQSRStandardHeaderLines(hInfo);
		if (null != header) {
			hInfo = VCFUtils.updateHeaderContigLines(hInfo, null, header.getSequenceDictionary(), true);
		}

		VCFHeader outHeader = new VCFHeader(hInfo);
		recalWriter = createVCFWriter(new File(options.getRecalFile()),outHeader);
		recalWriter.writeHeader(outHeader);
	}

	public static VariantContextWriter createVCFWriter(final File outFile,VCFHeader header) {
		Utils.nonNull(outFile);
		VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder().clearOptions()
				.setOutputFile(outFile);

		if (lenientVCFProcessing) {
			vcWriterBuilder = vcWriterBuilder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
		}

		if (createOutputVariantIndex && null != header) {
			vcWriterBuilder = vcWriterBuilder.setOption(Options.INDEX_ON_THE_FLY);
		}

		if (createOutputVariantMD5) {
			vcWriterBuilder.setCreateMD5();
		}

		if (null != header) {
			vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(header.getSequenceDictionary());
		}

		return vcWriterBuilder.build();
	}

	public static VariantContextWriter createVCFWriter(final File outFile,
			final SAMSequenceDictionary referenceDictionary, final boolean createMD5, final Options... options) {
		Utils.nonNull(outFile);

		VariantContextWriterBuilder vcWriterBuilder = new VariantContextWriterBuilder().clearOptions()
				.setOutputFile(outFile);

		if (createMD5) {
			vcWriterBuilder.setCreateMD5();
		}

		if (null != referenceDictionary) {
			vcWriterBuilder = vcWriterBuilder.setReferenceDictionary(referenceDictionary);
		}

		for (Options opt : options) {
			vcWriterBuilder = vcWriterBuilder.setOption(opt);
		}

		return vcWriterBuilder.build();
	}

	public void addTrainData(TrainData trainingData) {
		dataManager.addTrainingSet(trainingData);
	}
	
	public void setHeader(VCFHeader header){
		this.header = header;
	}

	public VariantDataManager getManager() {
		return this.dataManager;
	}

	public void apply(VariantContext vc) {
		if (variantLocusChanged(vc)) {
			consumeQueuedVariants();
		}
		if (counter % options.getSampleMode() == 0)
			variantsAtLocus.add(vc);
		counter++;
	}

	private boolean variantLocusChanged(final VariantContext nextVC) {
		if (variantsAtLocus.isEmpty()) {
			return false;
		} else {
			final VariantContext previous = variantsAtLocus.get(0);
			return (nextVC.getStart() != previous.getStart() || !nextVC.getContig().equals(previous.getContig()));
		}
	}

	private void consumeQueuedVariants() {
		variantsAtLocus.forEach(v -> addVariantDatum(v, true));
		variantsAtLocus.clear();
	}

	private void addVariantDatum(final VariantContext vc, final boolean isInput) {
		if (vc != null && (options.ignorAllFilter() || vc.isNotFiltered()
				|| options.getIgnoreInputFilters().containsAll(vc.getFilters()))) {
			if (VariantDataManager.checkVariationClass(vc, options.getArgument().MODE)
					&& !options.getArgument().useASannotations) {
				addDatum(reduceSum, isInput, vc, null, null);
			} else if (options.getArgument().useASannotations) {
				for (final Allele allele : vc.getAlternateAlleles()) {
					if (!GaeaVCFConstants.isSpanningDeletion(allele)
							&& VariantDataManager.checkVariationClass(vc, allele, options.getArgument().MODE)) {
						addDatum(reduceSum, isInput, vc, vc.getReference(), allele);
					}
				}
			}
		}
	}

	private void addDatum(final ExpandingArrayList<VariantDatum> variants, final boolean isInput,
			final VariantContext vc, final Allele refAllele, final Allele altAllele) {
		final VariantDatum datum = new VariantDatum();

		// Populate the datum with lots of fields from the VariantContext,
		// unfortunately the VC is too big so we just
		// pull in only the things we absolutely need.
		datum.referenceAllele = refAllele;
		datum.alternateAllele = altAllele;
		dataManager.decodeAnnotations(datum, vc, true);

		// non-deterministic because order of calls depends on load of machine
		datum.loc = (isInput ? parser.createGenomeLocation(vc) : null);

		datum.originalQual = vc.getPhredScaledQual();
		datum.isSNP = vc.isSNP() && vc.isBiallelic();
		datum.isTransition = datum.isSNP && VariantContextUtils.isTransition(vc);
		datum.isAggregate = !isInput;

		// Loop through the training data sets and if they overlap this locus
		// (and allele, if applicable) then update
		// the prior and training status appropriately. The locus used to find
		// training set variants is retrieved
		// by parseTrainingSets from the FeatureContext argument.
		dataManager.parseTrainingSets(parser,vc, datum, options.isTrustAllPolymorphic());
		final double priorFactor = QualityUtils.qualityToProbability(datum.prior);
		datum.prior = Math.log10(priorFactor) - Math.log10(1.0 - priorFactor);

		variants.add(datum);
	}

	public void addDatum(VariantDatum datum) {
		reduceSum.add(datum);
	}

	public Object traversal() {
		consumeQueuedVariants(); // finish processing any queued variants

		//try {
			dataManager.setData(reduceSum);
			dataManager.normalizeData(true); // Each data point is now (x -
												// mean) / standard deviation

			final GaussianMixtureModel goodModel;
			final GaussianMixtureModel badModel;

			final List<VariantDatum> positiveTrainingData = dataManager.getTrainingData();
			final List<VariantDatum> negativeTrainingData;

			// Generate the positive model using the training data and evaluate
			// each variant
			goodModel = engine.generateModel(positiveTrainingData, options.getArgument().MAX_GAUSSIANS);
			engine.evaluateData(dataManager.getData(), goodModel, false);
			// Generate the negative model using the worst performing data and
			// evaluate each variant contrastively
			negativeTrainingData = dataManager.selectWorstVariants();
			badModel = engine.generateModel(negativeTrainingData, Math
					.min(options.getArgument().MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, options.getArgument().MAX_GAUSSIANS));

			if (badModel.failedToConverge || goodModel.failedToConverge) {
				throw new UserException(
						"NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider "
								+ (badModel.failedToConverge
										? "raising the number of variants used to train the negative model (via --minNumBadVariants 5000, for example)."
										: "lowering the maximum number of Gaussians allowed for use in the model (via --maxGaussians 4, for example)."));
			}

			dataManager.dropAggregateData(); // Don't need the aggregate data
												// anymore so let's free up the
												// memory
			engine.evaluateData(dataManager.getData(), badModel, true);

			if (options.getOutputModel() != null) {
				final Report report = writeModelReport(goodModel, badModel, options.getUseAnnotations());
				try (final PrintStream modelReportStream = new PrintStream(options.getOutputModel())) {
					report.print(modelReportStream);
				} catch (FileNotFoundException e) {
					throw new UserException.CouldNotCreateOutputFile("File: (" + options.getOutputModel() + ")", e);
				}
			}

			engine.calculateWorstPerformingAnnotation(dataManager.getData(), goodModel, badModel);

			// Find the VQSLOD cutoff values which correspond to the various
			// tranches of calls requested by the user
			final int nCallsAtTruth = TrancheManager.countCallsAtTruth(dataManager.getData(), Double.NEGATIVE_INFINITY);
			final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric(nCallsAtTruth);

			if (!options.getScatterTranches()) {
				tranches = TrancheManager.findTranches(dataManager.getData(), options.getTsTranches(), metric,
						options.getArgument().MODE);
				if (tranchesStream != null) {
					tranchesStream.print(TruthSensitivityTranche.printHeader());
					tranchesStream.print(Tranche.tranchesString(tranches));
				}
			} /*else {
				tranches = TrancheManager.findVQSLODTranches(dataManager.getData(), options.getVQSLODTranches(), metric,
						options.getArgument().MODE);

				if (tranchesStream != null) {
					tranchesStream.print(VQSLODTranche.printHeader());
					tranchesStream.print(Tranche.tranchesString(tranches));
				}
			}*/

			if (recalWriter != null)
				dataManager.writeOutRecalibrationTable(recalWriter, header.getSequenceDictionary());
			if (options.getRScript() != null) {
				RScriptUtils.createVisualizationScript(options.getRScript(), engine, dataManager,
						dataManager.getRandomDataForPlotting(1000, positiveTrainingData, negativeTrainingData,
								dataManager.getEvaluationData()),
						goodModel, badModel, 0.0,
						dataManager.getAnnotationKeys().toArray(new String[options.getUseAnnotations().size()]));
			}

			return true;
		/*} catch (final Exception e) {
			throw new UserException(e.toString());
		}*/
	}

	public Report writeModelReport(final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel,
			final List<String> annotationList) {
		final String formatString = "%.16E";
		final Report report = new Report();

		if (dataManager != null) { // for unit test
			final double[] meanVector = dataManager.getMeanVector();
			final ReportTable annotationMeans = makeVectorTable("AnnotationMeans",
					"Mean for each annotation, used to normalize data", dataManager.annotationKeys, meanVector, "Mean",
					formatString);
			report.addTable(annotationMeans);

			final double[] varianceVector = dataManager.getVarianceVector(); // "varianceVector"
																				// is
																				// actually
																				// stdev
			final ReportTable annotationVariances = makeVectorTable("AnnotationStdevs",
					"Standard deviation for each annotation, used to normalize data", dataManager.annotationKeys,
					varianceVector, "StandardDeviation", // column header must
															// be one word
					formatString);
			report.addTable(annotationVariances);
		}

		final List<String> gaussianStrings = new ArrayList<>();
		final double[] pMixtureLog10s = new double[goodModel.getModelGaussians().size()];
		int idx = 0;

		for (final MultivariateGaussian gaussian : goodModel.getModelGaussians()) {
			pMixtureLog10s[idx] = gaussian.pMixtureLog10;
			gaussianStrings.add(Integer.toString(idx++));
		}

		final ReportTable goodPMix = makeVectorTable("GoodGaussianPMix", "Pmixture log 10 used to evaluate model",
				gaussianStrings, pMixtureLog10s, "pMixLog10", formatString, "Gaussian");
		report.addTable(goodPMix);

		gaussianStrings.clear();
		final double[] pMixtureLog10sBad = new double[badModel.getModelGaussians().size()];
		idx = 0;

		for (final MultivariateGaussian gaussian : badModel.getModelGaussians()) {
			pMixtureLog10sBad[idx] = gaussian.pMixtureLog10;
			gaussianStrings.add(Integer.toString(idx++));
		}
		final ReportTable badPMix = makeVectorTable("BadGaussianPMix", "Pmixture log 10 used to evaluate model",
				gaussianStrings, pMixtureLog10sBad, "pMixLog10", formatString, "Gaussian");
		report.addTable(badPMix);

		// The model and Gaussians don't know what the annotations are, so get
		// them from this class
		// VariantDataManager keeps the annotation in the same order as the
		// argument list
		final ReportTable positiveMeans = makeMeansTable("PositiveModelMeans",
				"Vector of annotation values to describe the (normalized) mean for each Gaussian in the positive model",
				annotationList, goodModel, formatString);
		report.addTable(positiveMeans);

		final ReportTable positiveCovariance = makeCovariancesTable("PositiveModelCovariances",
				"Matrix to describe the (normalized) covariance for each Gaussian in the positive model; covariance matrices are joined by row",
				annotationList, goodModel, formatString);
		report.addTable(positiveCovariance);

		// do the same for the negative model means
		final ReportTable negativeMeans = makeMeansTable("NegativeModelMeans",
				"Vector of annotation values to describe the (normalized) mean for each Gaussian in the negative model",
				annotationList, badModel, formatString);
		report.addTable(negativeMeans);

		final ReportTable negativeCovariance = makeCovariancesTable("NegativeModelCovariances",
				"Matrix to describe the (normalized) covariance for each Gaussian in the negative model; covariance matrices are joined by row",
				annotationList, badModel, formatString);
		report.addTable(negativeCovariance);

		return report;
	}

	private ReportTable makeVectorTable(final String tableName, final String tableDescription,
			final List<String> annotationList, final double[] perAnnotationValues, final String columnName,
			final String formatString, final String firstColumn) {
		final ReportTable vectorTable = new ReportTable(tableName, tableDescription, annotationList.size(),
				ReportTable.Sorting.DO_NOT_SORT);
		vectorTable.addColumn(firstColumn, "");
		vectorTable.addColumn(columnName, formatString);
		for (int i = 0; i < perAnnotationValues.length; i++) {
			vectorTable.addRowIDMapping(annotationList.get(i), i, true);
			vectorTable.set(i, 1, perAnnotationValues[i]);
		}
		return vectorTable;
	}

	private ReportTable makeMeansTable(final String tableName, final String tableDescription,
			final List<String> annotationList, final GaussianMixtureModel model, final String formatString) {
		final ReportTable meansTable = new ReportTable(tableName, tableDescription, annotationList.size(),
				ReportTable.Sorting.DO_NOT_SORT);
		meansTable.addColumn("Gaussian", "");
		for (final String annotationName : annotationList) {
			meansTable.addColumn(annotationName, formatString);
		}
		final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
		for (int i = 0; i < modelGaussians.size(); i++) {
			final MultivariateGaussian gaussian = modelGaussians.get(i);
			final double[] meanVec = gaussian.mu;
			if (meanVec.length != annotationList.size())
				throw new IllegalStateException(
						"Gaussian mean vector does not have the same size as the list of annotations");
			meansTable.addRowIDMapping(i, i, true);
			for (int j = 0; j < annotationList.size(); j++)
				meansTable.set(i, annotationList.get(j), meanVec[j]);
		}
		return meansTable;
	}

	private ReportTable makeCovariancesTable(final String tableName, final String tableDescription,
			final List<String> annotationList, final GaussianMixtureModel model, final String formatString) {
		final ReportTable modelCovariances = new ReportTable(tableName, tableDescription, annotationList.size() + 2,
				ReportTable.Sorting.DO_NOT_SORT); // +2 is for Gaussian and
													// Annotation columns
		modelCovariances.addColumn("Gaussian", "");
		modelCovariances.addColumn("Annotation", "");
		for (final String annotationName : annotationList) {
			modelCovariances.addColumn(annotationName, formatString);
		}
		final List<MultivariateGaussian> modelGaussians = model.getModelGaussians();
		for (int i = 0; i < modelGaussians.size(); i++) {
			final MultivariateGaussian gaussian = modelGaussians.get(i);
			final Matrix covMat = gaussian.sigma;
			if (covMat.getRowDimension() != annotationList.size()
					|| covMat.getColumnDimension() != annotationList.size())
				throw new IllegalStateException(
						"Gaussian covariance matrix does not have the same size as the list of annotations");
			for (int j = 0; j < annotationList.size(); j++) {
				modelCovariances.set(j + i * annotationList.size(), "Gaussian", i);
				modelCovariances.set(j + i * annotationList.size(), "Annotation", annotationList.get(j));
				for (int k = 0; k < annotationList.size(); k++) {
					modelCovariances.set(j + i * annotationList.size(), annotationList.get(k), covMat.get(j, k));

				}
			}
		}
		return modelCovariances;
	}

	protected ReportTable makeVectorTable(final String tableName, final String tableDescription,
			final List<String> annotationList, final double[] perAnnotationValues, final String columnName,
			final String formatString) {
		return makeVectorTable(tableName, tableDescription, annotationList, perAnnotationValues, columnName,
				formatString, "Annotation");
	}

	public final SAMSequenceDictionary getBestAvailableSequenceDictionary() {
		return header.getSequenceDictionary();
	}

	public List<TruthSensitivityTranche> getTranches() {
		return this.tranches;
	}

	public List<VariantDatum> getData() {
		Collections.sort(dataManager.getData(), VariantDatum.getComparator(header.getSequenceDictionary()));
		return dataManager.getData();
	}
	
	public List<String> getAnnotationKeys(){
		return dataManager.getAnnotationKeys();
	}
}
