package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.variant.statistic.VariantBasicStatistic;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.ToolsRunner;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.DBResource;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.FileResource;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.TrainData;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.AnnotationUtils;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.GaeaVCFHeaderLines;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class ApplyVQSR extends ToolsRunner{
	protected static final String LOW_VQSLOD_FILTER_NAME = "LOW_VQSLOD";

	private VariantRecalibrator recalibrator = null;

	private VCFQualityControlOptions options = null;

	final static private String listPrintSeparator = ",";

	final static private String arrayParseRegex = "[\\[\\]\\s]";

	final static private String emptyStringValue = "NA";

	final static private String emptyFloatValue = "NaN";

	private boolean foundSNPTranches = false;
	private boolean foundINDELTranches = false;

	private boolean useASannotations = false;

	private boolean EXCLUDE_FILTERED = false;

	private Double TS_FILTER_LEVEL = null;

	private VCFQualityControlArgumentCollection.Mode MODE = VCFQualityControlArgumentCollection.Mode.SNP;

	private List<TruthSensitivityTranche> tranches = new ArrayList<TruthSensitivityTranche>();

	protected Double VQSLOD_CUTOFF = null;

	private final double DEFAULT_VQSLOD_CUTOFF = 0.0;

	final static private String trancheFilterString = "VQSRTranche";

	private Set<String> ignoreInputFilterSet = new HashSet<String>();

	private boolean IGNORE_ALL_FILTERS = false;
	
	private VariantContextWriter vcfWriter;
	
	private VariantBasicStatistic basicStatics = null;
    
    private VariantBasicStatistic filteredStatics = null;
	
	public ApplyVQSR() {
		this.toolsDescription = "Gaea vcf quality control!\n";
	}

	public ApplyVQSR(VCFQualityControlOptions options) {
		this();
		this.options = options;
		
		VCFHeader header = null;
		
		try {
			header = getHeader();
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e.toString());
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
		recalibrator = new VariantRecalibrator(options,header);
		
		if(options.isStatics()){
			basicStatics = new VariantBasicStatistic(options.isIndividuals());
			filteredStatics = new VariantBasicStatistic(options.isIndividuals());
		}
	}
	
	public VCFHeader getHeader() throws FileNotFoundException, IOException{
		File f = new File(options.getInputs());
		AsciiLineReaderIterator iterator = null;
		if (f.getName().endsWith(".gz"))
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new GZIPInputStream(new FileInputStream(f))));
		else
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new FileInputStream(f)));
		VCFCodec codec = new VCFCodec();
		VCFHeader header = (VCFHeader)codec.readActualHeader(iterator);
		
		iterator.close();
		return header;
	}

	private void applyRecal() throws IOException {
		File f = new File(options.getInputs());
		AsciiLineReaderIterator iterator = null;
		if (f.getName().endsWith(".gz"))
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new GZIPInputStream(new FileInputStream(f))));
		else
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new FileInputStream(f)));
		VCFCodec codec = new VCFCodec();
		VCFHeader header = (VCFHeader)codec.readActualHeader(iterator);
		recalibrator.setHeader(header);
		
		for(String resource : options.getResources()) {
			TrainData trainData = new TrainData(options.getReference(), resource);
			if(trainData.isDB()) {
				trainData.setType(new DBResource());
			} else {
				trainData.setType(new FileResource());
			}
			trainData.initialize();
			recalibrator.addTrainData(trainData);
		}
		
		while (iterator.hasNext()) {
			VariantContext vc = codec.decode(iterator.next());
			recalibrator.apply(vc);
		}
		recalibrator.traversal();
		iterator.close();
	}
	
	public void apply() throws IOException{
		try {
			applyRecal();
		} catch (IOException e) {
			throw new UserException(e.toString());
		}
		
		File f = new File(options.getInputs());
		AsciiLineReaderIterator iterator = null;
		if (f.getName().endsWith(".gz"))
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new GZIPInputStream(new FileInputStream(f))));
		else
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new FileInputStream(f)));
		VCFCodec codec = new VCFCodec();
		VCFHeader header = (VCFHeader)codec.readActualHeader(iterator);
		
		initialize(header);
		
		while ( iterator.hasNext() ){
			VariantContext vc = codec.decode(iterator.next());
			if(options.isStatics())
				basicStatics.variantStatic(vc);
			applyVqsr(vc,recalibrator.getData());
		}
		
		iterator.close();
	}

	public void initialize(VCFHeader inputHeader) {
		final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

		final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
		hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
		hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.VQS_LOD_KEY));
		hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.CULPRIT_KEY));
		hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.POSITIVE_LABEL_KEY));
		hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.NEGATIVE_LABEL_KEY));
		if (useASannotations) {
			hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.AS_FILTER_STATUS_KEY));
			hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.AS_CULPRIT_KEY));
			hInfo.add(GaeaVCFHeaderLines.getInfoLine(GaeaVCFConstants.AS_VQS_LOD_KEY));
		}

		checkForPreviousApplyRecalRun(Collections.unmodifiableSet(inputHeaders));

		final TreeSet<String> samples = new TreeSet<>();
		samples.addAll(inputHeader.getGenotypeSamples());

		if (TS_FILTER_LEVEL != null) {
			for (final TruthSensitivityTranche t : recalibrator.getTranches()) {
				if (t.targetTruthSensitivity >= TS_FILTER_LEVEL) {
					tranches.add(t);
				}
			}
			Collections.reverse(tranches); 
		}

		if (TS_FILTER_LEVEL != null) {
			// if the user specifies both ts_filter_level and lodCutoff then
			// throw a user error
			if (VQSLOD_CUTOFF != null) {
				throw new UserException(
						"Arguments --ts_filter_level and --lodCutoff are mutually exclusive. Please only specify one option.");
			}

			if (tranches.size() >= 2) {
				for (int iii = 0; iii < tranches.size() - 1; iii++) {
					final TruthSensitivityTranche t = tranches.get(iii);
					hInfo.add(
							new VCFFilterHeaderLine(t.name,
									String.format("Truth sensitivity tranche level for " + t.model.toString()
											+ " model at VQS Lod: " + t.minVQSLod + " <= x < "
											+ tranches.get(iii + 1).minVQSLod)));
				}
			}
			if (tranches.size() >= 1) {
				hInfo.add(new VCFFilterHeaderLine(tranches.get(0).name + "+",
						String.format("Truth sensitivity tranche level for " + tranches.get(0).model.toString()
								+ " model at VQS Lod < " + tranches.get(0).minVQSLod)));
			} else {
				throw new UserException(
						"No tranches were found in the file or were above the truth sensitivity filter level "
								+ TS_FILTER_LEVEL);
			}
		} else {
			if (VQSLOD_CUTOFF == null) {
				VQSLOD_CUTOFF = DEFAULT_VQSLOD_CUTOFF;
			}
			hInfo.add(new VCFFilterHeaderLine(LOW_VQSLOD_FILTER_NAME, "VQSLOD < " + VQSLOD_CUTOFF));
		}
		
		final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter = VariantRecalibrator.createVCFWriter(new File(options.getOutputPath()),vcfHeader);
        vcfWriter.writeHeader(vcfHeader);
	}

	private void checkForPreviousApplyRecalRun(final Set<VCFHeaderLine> inputHeaders) {
		for (final VCFHeaderLine header : inputHeaders) {
			if (header instanceof VCFFilterHeaderLine) {
				final String filterName = ((VCFFilterHeaderLine) header).getID();
				if (filterName.length() < 12 || !filterName.substring(0, 11).equalsIgnoreCase(trancheFilterString)) {
					continue;
				}
				if (filterName.charAt(11) == 'S') {
					// for SNP tranches, get sensitivity limit
					final String sensitivityLimits = filterName.substring(14);
					if (trancheIntervalIsValid(sensitivityLimits))
						foundSNPTranches = true;
				} else if (filterName.charAt(11) == 'I') {
					// for INDEL tranches, get sensitivity limit
					final String sensitivityLimits = filterName.substring(16);
					if (trancheIntervalIsValid(sensitivityLimits))
						foundINDELTranches = true;
				}
			}
		}
	}

	private boolean trancheIntervalIsValid(final String sensitivityLimits) {
		final String[] vals = sensitivityLimits.split("to");
		if (vals.length != 2)
			return false;
		try {
			Double.parseDouble(vals[0]);
			Double.parseDouble(vals[1].replace("+", "")); 
		} catch (NumberFormatException e) {
			throw new UserException(
					"Poorly formatted tranche filter name does not contain two sensitivity interval end points.");
		}
		return true;
	}

	public void applyVqsr(VariantContext vc, List<VariantDatum> data) {
		final List<VariantDatum> recals = getDatum(data, vc.getStart());
		final boolean evaluateThisVariant = useASannotations || VariantDataManager.checkVariationClass(vc, MODE);

		// vc.isNotFiltered is true for PASS; vc.filtersHaveBeenApplied covers
		// PASS and filters
		final boolean variantIsNotFiltered = IGNORE_ALL_FILTERS || vc.isNotFiltered()
				|| (!ignoreInputFilterSet.isEmpty() && ignoreInputFilterSet.containsAll(vc.getFilters()));

		if (evaluateThisVariant && variantIsNotFiltered) {
			String filterString;
			final VariantContextBuilder builder = new VariantContextBuilder(vc);
			if (!useASannotations) {
				filterString = doSiteSpecificFiltering(vc, recals, builder);
			} else { // allele-specific mode
				filterString = doAlleleSpecificFiltering(vc, recals, builder);
			}

			// for both non-AS and AS modes:
			if (filterString.equals(VCFConstants.PASSES_FILTERS_v4)) {
				builder.passFilters();
			} else if (filterString.equals(VCFConstants.UNFILTERED)) {
				builder.unfiltered();
			} else {
				builder.filters(filterString);
			}

			final VariantContext outputVC = builder.make();
			if (!EXCLUDE_FILTERED || outputVC.isNotFiltered()) {
				if(outputVC.isNotFiltered() && options.isStatics())
					filteredStatics.variantStatic(outputVC);
				vcfWriter.add(outputVC);
			}
		} else { // valid VC but not compatible with this mode, so just emit the
					// variant untouched
			if(vc.isNotFiltered() && options.isStatics())
				filteredStatics.variantStatic(vc);
			vcfWriter.add(vc);
		}
	}

	private String doSiteSpecificFiltering(final VariantContext vc, final List<VariantDatum> recals,
			final VariantContextBuilder builder) {
		VariantDatum recalDatum = getMatchingRecalVC(vc, recals, null);
		if (recalDatum == null) {
			throw new UserException("Encountered input variant which isn't found in the input recal file. "
					+ "Please make sure VariantRecalibrator and ApplyRecalibration were "
					+ "run on the same set of input variants. First seen at: " + vc);
		}

		final String lodString = String.valueOf(recalDatum.lod);
		if (lodString == null) {
			throw new UserException(
					"Encountered a malformed record in the input recal file. There is no lod for the record at: " + vc);
		}
		final double lod;
		try {
			lod = Double.valueOf(lodString);
		} catch (NumberFormatException e) {
			throw new UserException(
					"Encountered a malformed record in the input recal file. The lod is unreadable for the record at: "
							+ vc);
		}

		builder.attribute(GaeaVCFConstants.VQS_LOD_KEY, lod);
		builder.attribute(GaeaVCFConstants.CULPRIT_KEY, recalDatum.worstAnnotation != -1
				? recalibrator.getAnnotationKeys().get(recalDatum.worstAnnotation) : "NULL");
		if (recalDatum != null) {
			if (recalDatum.atTrainingSite)
				builder.attribute(GaeaVCFConstants.POSITIVE_LABEL_KEY, true);
			if (recalDatum.atAntiTrainingSite)
				builder.attribute(GaeaVCFConstants.NEGATIVE_LABEL_KEY, true);
		}

		return generateFilterString(lod);
	}

	private List<VariantDatum> getDatum(List<VariantDatum> data, int start) {
		return data.stream().filter(d -> d.loc.getStart() == start).collect(Collectors.toList());
	}

	private VariantDatum getMatchingRecalVC(final VariantContext target, final List<VariantDatum> recalVCs,
			final Allele allele) {
		for (final VariantDatum recalVC : recalVCs) {
			if (target.getEnd() == recalVC.loc.getStop()) {
				if (!useASannotations)
					return recalVC;
				else if (allele.equals(recalVC.referenceAllele))
					return recalVC;
			}
		}
		return null;
	}

	protected String generateFilterString(final double lod) {
		String filterString = null;
		if (TS_FILTER_LEVEL != null) {
			for (int i = tranches.size() - 1; i >= 0; i--) {
				final TruthSensitivityTranche tranche = tranches.get(i);
				if (lod >= tranche.minVQSLod) {
					if (i == tranches.size() - 1) {
						filterString = VCFConstants.PASSES_FILTERS_v4;
					} else {
						filterString = tranche.name;
					}
					break;
				}
			}

			if (filterString == null) {
				filterString = tranches.get(0).name + "+";
			}
		} else {
			filterString = (lod < VQSLOD_CUTOFF ? LOW_VQSLOD_FILTER_NAME : VCFConstants.PASSES_FILTERS_v4);
		}

		return filterString;
	}

	private String doAlleleSpecificFiltering(final VariantContext vc, final List<VariantDatum> recals,
			final VariantContextBuilder builder) {
		double bestLod = VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE;
		final List<String> culpritStrings = new ArrayList<>();
		final List<String> lodStrings = new ArrayList<>();
		final List<String> AS_filterStrings = new ArrayList<>();

		String[] prevCulpritList = null;
		String[] prevLodList = null;
		String[] prevASfiltersList = null;

		// get VQSR annotations from previous run of ApplyRecalibration, if
		// applicable
		if (foundINDELTranches || foundSNPTranches) {
			final String prevCulprits = vc.getAttributeAsString(GaeaVCFConstants.AS_CULPRIT_KEY, "");
			prevCulpritList = prevCulprits.isEmpty() ? new String[0] : prevCulprits.split(listPrintSeparator);
			final String prevLodString = vc.getAttributeAsString(GaeaVCFConstants.AS_VQS_LOD_KEY, "");
			prevLodList = prevLodString.isEmpty() ? new String[0] : prevLodString.split(listPrintSeparator);
			final String prevASfilters = vc.getAttributeAsString(GaeaVCFConstants.AS_FILTER_STATUS_KEY, "");
			prevASfiltersList = prevASfilters.isEmpty() ? new String[0] : prevASfilters.split(listPrintSeparator);
		}

		// for each allele in the current VariantContext
		for (int altIndex = 0; altIndex < vc.getNAlleles() - 1; altIndex++) {
			final Allele allele = vc.getAlternateAllele(altIndex);

			// if the current allele is not part of this recalibration mode, add
			// its annotations to the list and go to the next allele
			if (!VariantDataManager.checkVariationClass(vc, allele, MODE)) {
				updateAnnotationsWithoutRecalibrating(altIndex, prevCulpritList, prevLodList, prevASfiltersList,
						culpritStrings, lodStrings, AS_filterStrings);
				continue;
			}

			// if the current allele does need to have recalibration applied...

			// initialize allele-specific VQSR annotation data with values for
			// spanning deletion
			String alleleLodString = emptyFloatValue;
			String alleleFilterString = emptyStringValue;
			String alleleCulpritString = emptyStringValue;

			// if it's not a spanning deletion, replace those allele strings
			// with the real values
			if (!GaeaVCFConstants.isSpanningDeletion(allele)) {
				VariantDatum recalDatum = getMatchingRecalVC(vc, recals, allele);
				if (recalDatum == null) {
					throw new UserException("Encountered input allele which isn't found in the input recal file. "
							+ "Please make sure VariantRecalibrator and ApplyRecalibration were "
							+ "run on the same set of input variants with flag -AS. First seen at: " + vc);
				}

				// compare VQSLODs for all alleles in the current mode for
				// filtering later
				final double lod = recalDatum.lod;
				if (lod > bestLod)
					bestLod = lod;

				alleleLodString = String.format("%.4f", lod);
				alleleFilterString = generateFilterString(lod);
				alleleCulpritString = recalDatum.worstAnnotation != -1
						? recalibrator.getAnnotationKeys().get(recalDatum.worstAnnotation) : ".";

				if (recalDatum != null) {
					if (recalDatum.atTrainingSite)
						builder.attribute(GaeaVCFConstants.POSITIVE_LABEL_KEY, true);
					if (recalDatum.atAntiTrainingSite)
						builder.attribute(GaeaVCFConstants.NEGATIVE_LABEL_KEY, true);
				}
			}

			// append per-allele VQSR annotations
			lodStrings.add(alleleLodString);
			AS_filterStrings.add(alleleFilterString);
			culpritStrings.add(alleleCulpritString);
		}

		// Annotate the new record with its VQSLOD, AS_FilterStatus, and the
		// worst performing annotation
		if (!AS_filterStrings.isEmpty())
			builder.attribute(GaeaVCFConstants.AS_FILTER_STATUS_KEY,
					AnnotationUtils.encodeStringList(AS_filterStrings));
		if (!lodStrings.isEmpty())
			builder.attribute(GaeaVCFConstants.AS_VQS_LOD_KEY, AnnotationUtils.encodeStringList(lodStrings));
		if (!culpritStrings.isEmpty())
			builder.attribute(GaeaVCFConstants.AS_CULPRIT_KEY, AnnotationUtils.encodeStringList(culpritStrings));

		return generateFilterStringFromAlleles(vc, bestLod);
	}

	protected String generateFilterStringFromAlleles(final VariantContext vc, final double bestLod) {
		String filterString = ".";

		final boolean bothModesWereRun = (MODE == VCFQualityControlArgumentCollection.Mode.SNP && foundINDELTranches)
				|| (MODE == VCFQualityControlArgumentCollection.Mode.INDEL && foundSNPTranches);
		final boolean onlyOneModeNeeded = !vc.isMixed() && VariantDataManager.checkVariationClass(vc, MODE);

		// if both SNP and INDEL modes have not yet been run (and need to be),
		// leave this variant as unfiltered and add the filters for the alleles
		// in this mode to the INFO field
		if (!bothModesWereRun && !onlyOneModeNeeded) {
			return VCFConstants.UNFILTERED;
		}

		// if both SNP and INDEL modes have been run or the site is not mixed,
		// generate a filter string for this site based on both models
		// pull out the allele filter status from the info field (there may be
		// more than one entry in the list if there were multiple snp/indel
		// alleles assessed in the other mode)
		final String prevFilterStatus = vc.getAttributeAsString(GaeaVCFConstants.AS_FILTER_STATUS_KEY, null);

		// if this site hasn't had a filter applied yet
		if (prevFilterStatus != null && !prevFilterStatus.equals(VCFConstants.UNFILTERED)) {
			final String prevAllelesFilterStatusString = vc.getAttributeAsString(GaeaVCFConstants.AS_FILTER_STATUS_KEY,
					null);
			final String[] prevAllelesFilterStatusList = prevAllelesFilterStatusString.split(listPrintSeparator);
			// start with the current best allele filter as the most lenient
			// filter across all modes and all alleles
			String mostLenientFilterName = generateFilterString(bestLod);
			// if the current mode's best allele passes the tranche filter, then
			// let the whole site pass
			if (mostLenientFilterName.equals(VCFConstants.PASSES_FILTERS_v4)) {
				filterString = mostLenientFilterName;
			}
			// if the current mode's best allele does not pass the tranche
			// filter, compare the most lenient filter of this mode with those
			// from the previous mode
			else {
				double mostLenientSensitivityLowerLimit = parseFilterLowerLimit(mostLenientFilterName);
				for (int i = 0; i < prevAllelesFilterStatusList.length; i++) {
					final String alleleFilterString = prevAllelesFilterStatusList[i].replaceAll(arrayParseRegex, "")
							.trim();
					// if any allele from the previous mode passed the tranche
					// filter, then let the whole site pass
					if (alleleFilterString.equals(VCFConstants.PASSES_FILTERS_v4)) {
						mostLenientFilterName = alleleFilterString;
						break;
					}
					// if there's no PASS, then we need to parse the filters to
					// find out how lenient they are
					else {
						final double alleleLowerLimit = parseFilterLowerLimit(alleleFilterString);
						if (alleleLowerLimit == -1)
							continue;
						if (alleleLowerLimit < mostLenientSensitivityLowerLimit) {
							mostLenientSensitivityLowerLimit = alleleLowerLimit;
							mostLenientFilterName = alleleFilterString;
						}
					}

				}
				filterString = mostLenientFilterName;
			}
		}

		// if both modes have been run, but the previous mode didn't apply a
		// filter, use the current mode's best allele VQSLOD filter (shouldn't
		// get run, but just in case)
		else {
			filterString = generateFilterString(bestLod);
		}

		return filterString;
	}

	public double parseFilterLowerLimit(final String trancheFilter) {
		final Pattern pattern = Pattern.compile("VQSRTranche\\S+(\\d+\\.\\d+)to(\\d+\\.\\d+)");
		final Matcher m = pattern.matcher(trancheFilter);
		return m.find() ? Double.parseDouble(m.group(1)) : -1;
	}

	private void updateAnnotationsWithoutRecalibrating(final int altIndex, final String[] prevCulpritList,
			final String[] prevLodList, final String[] prevASfiltersList, final List<String> culpritString,
			final List<String> lodString, final List<String> AS_filterString) {
		if (foundINDELTranches || foundSNPTranches) {
			if (altIndex < prevCulpritList.length) {
				culpritString.add(prevCulpritList[altIndex].replaceAll(arrayParseRegex, "").trim());
				lodString.add(prevLodList[altIndex].replaceAll(arrayParseRegex, "").trim());
				AS_filterString.add(prevASfiltersList[altIndex].replaceAll(arrayParseRegex, "").trim());
			}
		} else { // if the other allele type hasn't been processed yet, make
					// sure there are enough entries
			culpritString.add(emptyStringValue);
			lodString.add(emptyFloatValue);
			AS_filterString.add(emptyStringValue);
		}
	}
	
	private void reportWriter() throws IOException{
		if(!options.isStatics())
			return;
		String sampleName = basicStatics.getSampleName();
		BufferedWriter os = new BufferedWriter(new FileWriter(options.getStaticsPath()+"/"+sampleName));
		os.write("before filter\n");
		os.write(basicStatics.toString());
		os.write("\n");
		os.write("after filter\n");
		os.write(filteredStatics.toString());
		os.close();
	}

	@Override
	public int run(String[] args) throws Exception {
		VCFQualityControlOptions options = new VCFQualityControlOptions();
		options.parse(args);
		ApplyVQSR applyVQSR = new ApplyVQSR(options);
		applyVQSR.apply();
		applyVQSR.reportWriter();
		return 0;
	}
}
