package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.AlleleSpecificAnnotationData;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.QualityUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public final class AS_RMSMappingQuality extends InfoFieldAnnotation implements AS_StandardAnnotation, ReducibleAnnotation {

    private final String printFormat = "%.2f";

    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";


    @Override
    public Map<String, Object> annotate(final ChromosomeInformationShare ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @Override
    public Map<String, Object> annotateRawData(final ChromosomeInformationShare ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {
        Utils.nonNull(vc);
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new LinkedHashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    /**
     * For AS_RMSMappingQuality annotations, the annotations will simply consist of a list of the total value for
     * every allele computed by parsing all of the individual AS_RMSMappingQuality Raw Key values as doubles
     * and totaling them.
     */
    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<ReducibleAnnotationData<?>>  annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData<Number> combinedData = new AlleleSpecificAnnotationData(vcAlleles, null);

        for (final ReducibleAnnotationData<?> currentValue : annotationList) {
            ReducibleAnnotationData<Number> value = (ReducibleAnnotationData<Number>)currentValue;
            parseRawDataString(value);
            combineAttributeMap(value, combinedData);

        }
        final Map<String, Object> annotations = new HashMap<>();
        String annotationString = makeRawAnnotationString(vcAlleles, combinedData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    public void combineAttributeMap(final ReducibleAnnotationData<Number> toAdd, final ReducibleAnnotationData<Number> combined) {
        //check that alleles match
        for (final Allele currentAllele : combined.getAlleles()){
            //combined is initialized with all alleles, but toAdd might have only a subset
            if(toAdd.getAttribute(currentAllele) != null) {
                if (toAdd.getAttribute(currentAllele) != null && combined.getAttribute(currentAllele) != null) {
                    combined.putAttribute(currentAllele, (double) combined.getAttribute(currentAllele) + (double) toAdd.getAttribute(currentAllele));
                } else {
                    combined.putAttribute(currentAllele, toAdd.getAttribute(currentAllele));
                }
            }
        }
    }

    protected void parseRawDataString(final ReducibleAnnotationData<Number> myData) {
        final String rawDataString = myData.getRawData();
        //get per-allele data by splitting on allele delimiter
        final String[] rawDataPerAllele = rawDataString.split(SPLIT_DELIM);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            final String alleleData = rawDataPerAllele[i];
            myData.putAttribute(myData.getAlleles().get(i), Double.parseDouble(alleleData));
        }
    }


    /**
     * Takes combined raw annotation data sums, and calculates per allele the average root mean squared from the raw data
     * using expected Allele Depth counts data.
     *
     * Will output delineated doubles in the format: sqrt(TotalAllele1RMS/Allele1Depth)|sqrt(TotalAllele1RMS/Allele1Depth)|...
     */
    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getRawKeyName())) {
            return new HashMap<>();
        }
        final String rawMQdata = vc.getAttributeAsString(getRawKeyName(),null);
        if (rawMQdata == null) {
            return new HashMap<>();
        }

        final Map<String,Object> annotations = new HashMap<>();
        final ReducibleAnnotationData myData = new AlleleSpecificAnnotationData<Double>(originalVC.getAlleles(), rawMQdata);
        parseRawDataString(myData);

        final String annotationString = makeFinalizedAnnotationString(vc, myData.getAttributeMap());
        annotations.put(getKeyNames().get(0), annotationString);
        return annotations;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData myData){
        //For the raw data here, we're only keeping track of the sum of the squares of our values
        //When we go to reduce, we'll use the AD info to get the number of reads

        //must use likelihoods for allele-specific annotations
        if (likelihoods == null) {
            return;
        }
        getRMSDataFromLikelihoods(likelihoods, myData);
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GaeaVCFConstants.AS_RMS_MAPPING_QUALITY_KEY); }

    @Override
    public String getRawKeyName() { return GaeaVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY; }

    private void getRMSDataFromLikelihoods(final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData<Number> myData) {
        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles() ) {
            if (bestAllele.isInformative()) {
                final int mq = bestAllele.read.getMappingQuality();
                if ( mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE ) {
                    final double currSquareSum = myData.hasAttribute(bestAllele.allele) ? (double) myData.getAttribute(bestAllele.allele) : 0;
                    myData.putAttribute(bestAllele.allele, currSquareSum + mq * mq);
                }
            }
        }
    }

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Number> perAlleleValues) {
        String annotationString = "";
        for (final Allele current : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += PRINT_DELIM;
            }
            if(perAlleleValues.get(current) != null) {
                annotationString += String.format(printFormat, perAlleleValues.get(current));
            } else {
                annotationString += String.format(printFormat, 0.0);
            }
        }
        return annotationString;
    }

    private String makeFinalizedAnnotationString(final VariantContext vc, final Map<Allele, Number> perAlleleValues) {
        final Map<Allele, Integer> variantADs = getADcounts(vc);
        String annotationString = "";
        for (final Allele current : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty()) {
                annotationString += ",";
            }
            if (perAlleleValues.containsKey(current)) {
                annotationString += String.format(printFormat, Math.sqrt((double) perAlleleValues.get(current) / variantADs.get(current)));
            }
        }
        return annotationString;
    }

    private Map<Allele, Integer> getADcounts(final VariantContext vc) {
        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 ) {
            return null;
        }

        final Map<Allele, Integer> variantADs = new HashMap<>();
        for(final Allele a : vc.getAlleles())
            variantADs.put(a,0);

        for (final Genotype gt : vc.getGenotypes()) {
            if(gt.hasAD()) {
                final int[] ADs = gt.getAD();
                for (int i = 1; i < vc.getNAlleles(); i++) {
                    variantADs.put(vc.getAlternateAllele(i - 1), variantADs.get(vc.getAlternateAllele(i - 1)) + ADs[i]); //here -1 is to reconcile allele index with alt allele index
                }
            }
        }
        return variantADs;
    }
}

