package org.bgi.flexlab.gaea.tools.jointcalling.annotator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ReducibleAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaVcfHeaderLines;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public abstract class RMSAnnotation extends InfoFieldAnnotation implements ReducibleAnnotation {

	@Override
	public String getRawKeyName() {
		return null;
	}

	@Override
	public Map<String, Object> annotateRawData(RefMetaDataTracker tracker, ChromosomeInformationShare ref,
			VariantContext vc) {
		return new HashMap<>();
	}

	@Override
	public Map<String, Object> combineRawData(List<Allele> allelesList,
			List<? extends ReducibleAnnotationData> listOfRawData) {
		ReducibleAnnotationData combinedData = new ReducibleAnnotationData(null);

		for (final ReducibleAnnotationData currentValue : listOfRawData) {
			parseRawDataString(currentValue);
			combineAttributeMap(currentValue, combinedData);

		}
		final Map<String, Object> annotations = new HashMap<>();
		String annotationString = makeRawAnnotationString(allelesList, combinedData.getAttributeMap());
		annotations.put(getRawKeyName(), annotationString);
		return annotations;
	}

	abstract protected String makeRawAnnotationString(List<Allele> vcAlleles, Map<Allele, Number> sumOfSquares);

	abstract protected String makeFinalizedAnnotationString(VariantContext vc, Map<Allele, Number> sumOfSquares);

	protected void parseRawDataString(ReducibleAnnotationData<Number> myData) {
		final String rawDataString = myData.getRawData();
		String[] rawMQdataAsStringVector;
		rawMQdataAsStringVector = rawDataString.split(",");
		double squareSum = Double.parseDouble(rawMQdataAsStringVector[0]);
		myData.putAttribute(Allele.NO_CALL, squareSum);
	}

	public void combineAttributeMap(ReducibleAnnotationData<Number> toAdd, ReducibleAnnotationData<Number> combined) {
		if (combined.getAttribute(Allele.NO_CALL) != null)
			combined.putAttribute(Allele.NO_CALL,
					(Double) combined.getAttribute(Allele.NO_CALL) + (Double) toAdd.getAttribute(Allele.NO_CALL));
		else
			combined.putAttribute(Allele.NO_CALL, toAdd.getAttribute(Allele.NO_CALL));

	}

	@Override
	public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
		if (!vc.hasAttribute(getRawKeyName()))
			return new HashMap<>();
		String rawMQdata = vc.getAttributeAsString(getRawKeyName(), null);
		if (rawMQdata == null)
			return new HashMap<>();

		ReducibleAnnotationData myData = new ReducibleAnnotationData(rawMQdata);
		parseRawDataString(myData);

		String annotationString = makeFinalizedAnnotationString(vc, myData.getAttributeMap());
		return Collections.singletonMap(getKeyNames().get(0), (Object) annotationString);
	}

	@Override
	public void calculateRawData(VariantContext vc, ReducibleAnnotationData rawAnnotations) {
	}

	@Override
	public Map<String, Object> annotate(RefMetaDataTracker tracker, ChromosomeInformationShare ref, VariantContext vc) {
		return null;
	}

	@Override
	public List<VCFInfoHeaderLine> getDescriptions() {
		final List<VCFInfoHeaderLine> headerLines = new ArrayList<>();

		headerLines.add(GaeaVcfHeaderLines.getInfoLine(getKeyNames().get(0)));
		return headerLines;
	}

	public int getNumOfReads(final VariantContext vc) {
		// don't use the full depth because we don't calculate MQ for reference
		// blocks
		int numOfReads = 0;
		if (vc.hasAttribute(VCFConstants.DEPTH_KEY)) {
			numOfReads += Integer.parseInt(vc.getAttributeAsString(VCFConstants.DEPTH_KEY, "-1"));
			if (vc.hasGenotypes()) {
				for (Genotype gt : vc.getGenotypes()) {
					if (gt.isHomRef()) {
						// site-level DP contribution will come from MIN_DP for
						// gVCF-called reference variants or DP for BP
						// resolution
						if (gt.hasExtendedAttribute("MIN_DP"))
							numOfReads -= Integer.parseInt(gt.getExtendedAttribute("MIN_DP").toString());
						else if (gt.hasDP())
							numOfReads -= gt.getDP();
					}

				}
			}
			return numOfReads;
		}
		return -1;
	}
}
