package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.exception.UserException;

import htsjdk.variant.variantcontext.Allele;

public class ReducibleAnnotationData<T> {
    protected String rawData;
    protected Map<Allele, T> attributeMap;

    /**
     * Create a new ReducibleAnnotationData using the raw data string from a VCF
     * @param inputData the raw data as read in from a VCF
     */
    public ReducibleAnnotationData(final String inputData) {
        rawData = inputData; attributeMap = new HashMap<>();
        attributeMap.put(Allele.NO_CALL, null);
    }

    /**
     *
     * @return the string of raw data as represented in the VCF
     */
    public String getRawData() {return rawData;}

    /**
     * Note: parent class ReducibleAnnotationData is non-allele specific and stores all values with the no-call allele
     * @return the list of alleles for which we have raw annotation data
     */
    public List<Allele> getAlleles() {
        List ret = new ArrayList<Allele>();
        ret.addAll(attributeMap.keySet());
        return ret;
    }

    /**
     *
     * @param key   the allele of interest
     * @return  do we have data for the allele of interest?
     */
    public boolean hasAttribute(Allele key) {
        return attributeMap.containsKey(key);
    }

    /**
     *
     * @param key the allele of interest
     * @return  data for the allele of interest
     */
    public T getAttribute(Allele key) {
        return attributeMap.get(key);
    }

    /**
     *
     * @param key   the allele of interest
     * @param value raw data corresponding to the allele of interest
     */
    public void putAttribute(Allele key, T value) {
        attributeMap.put(key, value);
    }

    /**
     * Assign all of the per-allele raw data at once
     * @param inputMap  the pre-calculated per-allele data
     */
    public void setAttributeMap(Map<Allele, T> inputMap) {
        attributeMap = inputMap;
    }

    /**
     * Get the stored raw per-allele data
     * @return
     */
    public Map<Allele, T> getAttributeMap() {return Collections.unmodifiableMap(attributeMap);}

    public void validateAllelesList() {
        boolean foundRef = false;
        for (final Allele a : this.getAlleles()) {
            if (a.isReference()) {
                if (foundRef) {
                    throw new UserException("ERROR: multiple reference alleles found in annotation data\n");
                }
                foundRef = true;
            }
        }
        if (!foundRef) {
            throw new UserException("ERROR: no reference alleles found in annotation data\n");
        }
    }
}

