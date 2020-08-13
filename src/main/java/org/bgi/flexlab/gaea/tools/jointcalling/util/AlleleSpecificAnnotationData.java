package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;

import htsjdk.variant.variantcontext.Allele;

public class AlleleSpecificAnnotationData<T> extends ReducibleAnnotationData<T>{
    final private List<Allele> alleleList;
    private Allele refAllele;

    public AlleleSpecificAnnotationData(final List<Allele> inputAlleles, final String inputData) {
        super(inputData);
        attributeMap = new HashMap<>();
        for(final Allele a : inputAlleles) {
            attributeMap.put(a, null);
        }
        alleleList = inputAlleles;
        for(Allele a : alleleList) {
            if(a.isReference()) {
                refAllele = a;
            }
        }
    }

    @Override
    public List<Allele> getAlleles() {return Collections.unmodifiableList(alleleList);}

    /**
     * Get the reference allele for this allele-specific data.
     * (Used in cases where annotations compare some attribute of the alt alleles to that of the reference.)
     * @return  the reference allele for this data
     */
    public Allele getRefAllele() {return refAllele;}

    public void setAttributeMap(Map<Allele, T> inputMap) {
        super.setAttributeMap(inputMap);
        checkRefAlleles();
    }

    private void checkRefAlleles() {
        boolean foundRef = false;
        for (Allele a : alleleList) {
            if (a.isReference()) {
                if (foundRef)
                    throw new UserException("ERROR: multiple reference alleles found in annotation data\n");
                foundRef = true;
            }
        }
        if (!foundRef)
            throw new UserException("ERROR: no reference alleles found in annotation data\n");
    }

    public String makeRawAnnotationString(String printDelim) {
        String annotationString = "";
        for (final Allele current : alleleList) {
            if (!annotationString.isEmpty())
                annotationString += printDelim;
            if(attributeMap.get(current) != null)
                annotationString += attributeMap.get(current).toString();
        }
        return annotationString.replaceAll("[\\[\\]\\s]", "");
    }
}
