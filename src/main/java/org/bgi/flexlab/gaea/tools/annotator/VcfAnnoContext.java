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
package org.bgi.flexlab.gaea.tools.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.ArrayUtils;
import org.bgi.flexlab.gaea.tools.annotator.db.AnnoDBQuery;
import org.bgi.flexlab.gaea.tools.annotator.interval.Chromosome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Variant;
import org.bgi.flexlab.gaea.tools.annotator.realignment.VcfRefAltAlign;

import java.util.*;
import java.util.stream.Collectors;

public class VcfAnnoContext {

    public static String ANN_FIELD_NAMES[] = { //
            "FLKSEQ", //
            "CHROM", //
            "REF", //
            "START", //
            "END", //
            "POS", //
    };

    private String contig;
    private String refStr;
    private int start;
    private int end;
    private List<String> alts;
    private Map<String, List<VariantContext>> variantContextMap;
    protected LinkedList<Variant> variants;
    private List<AnnotationContext> annotationContexts = null;
    private Map<String, SampleAnnotationContext> sampleAnnoContexts;
    private Map<String, String> annoItems;
    private String annoStr;
    private String flankSeq;

    public VcfAnnoContext(){
        variants = new LinkedList<>();
        sampleAnnoContexts = new HashMap<>();
    }

    public VcfAnnoContext(VariantContext variantContext){
        variants = new LinkedList<>();
        sampleAnnoContexts = new HashMap<>();
        variantContextMap = new HashMap<>();
        init(variantContext);
    }

    public VcfAnnoContext(VariantContext variantContext, String filename){
        variants = new LinkedList<>();
        sampleAnnoContexts = new HashMap<>();
        variantContextMap = new HashMap<>();
        List<VariantContext> variantContexts = new ArrayList<>();
        variantContexts.add(variantContext);
        variantContextMap.put(filename, variantContexts);
        init(variantContext);
    }

    public void init(VariantContext variantContext){
        contig = variantContext.getContig();
        refStr = variantContext.getReference().getBaseString();
        start = variantContext.getStart();
        end = variantContext.getEnd();
        alts = new ArrayList<>();
        for (Allele allele : variantContext.getAlternateAlleles()) {
            if(! Allele.wouldBeStarAllele(allele.getBases()) )
                alts.add(allele.getBaseString());
        }
        addSampleContext(variantContext);
    }

    public void add(VariantContext variantContext, String filename){
        if(variantContextMap.containsKey(filename))
            variantContextMap.get(filename).add(variantContext);
        else {
            List<VariantContext> variantContexts = new ArrayList<>();
            variantContexts.add(variantContext);
            variantContextMap.put(filename, variantContexts);
        }
        for(Allele allele: variantContext.getAlternateAlleles()){
            if(Allele.wouldBeStarAllele(allele.getBases()))
                continue;
            if(!alts.contains(allele.getBaseString()))
                alts.add(allele.getBaseString());
        }
        addSampleContext(variantContext);
    }


    private void addSampleContext(VariantContext variantContext){
        for(String sampleName : variantContext.getSampleNamesOrderedByName())
        {
            Genotype gt = variantContext.getGenotype(sampleName);

            if(!isVar(gt))
                continue;

            if(hasSample(sampleName)){
                SampleAnnotationContext sampleAnnoContext = sampleAnnoContexts.get(sampleName);
                sampleAnnoContext.add(variantContext);
            }else {
                SampleAnnotationContext sampleAnnoContext = new SampleAnnotationContext(sampleName, variantContext);
                sampleAnnoContexts.put(sampleName, sampleAnnoContext);
            }
        }
    }

    public boolean hasSample(String sampleName){
        return sampleAnnoContexts.containsKey(sampleName);
    }

    /**
     * Create a list of variants from this variantContext
     */
    public List<Variant> variants(Genome genome) {
//        if (!variants.isEmpty()) return variants;
        Chromosome chr = genome.getChromosome(contig);

        // interval 使用 0-base 方式建立，应使用start - 1创建variant对象
        for (String alt : alts) {
            Variant variant = createVariant(chr, (int)start - 1, refStr, alt, "");
            if(variant == null)
                continue;
            variants.add(variant);
        }
        return variants;
    }

    /**
     * Create a variant
     */
    Variant createVariant(Chromosome chromo, int start, String reference, String alt, String id) {
        Variant variant = null;
        if (alt != null) alt = alt.toUpperCase();

        if (alt == null || alt.isEmpty() || alt.equals(reference)) {
            // Non-variant
            variant = Variant.create(chromo, start, reference, null, id);
        } else if (alt.charAt(0) == '<') {
            // TODO Structural variants
            System.err.println("Cann't annotate Structural variants! ");
        } else if ((alt.indexOf('[') >= 0) || (alt.indexOf(']') >= 0)) {
            // TODO Translocations
            System.err.println("Cann't annotate Translocations: ALT has \"[\" or \"]\" info!");

        } else if (reference.length() == alt.length()) {
            // Case: SNP, MNP
            if (reference.length() == 1) {
                // SNPs
                // 20     3 .         C      G       .   PASS  DP=100
                variant = Variant.create(chromo, start, reference, alt, id);
            } else {
                // MNPs
                // 20     3 .         TC     AT      .   PASS  DP=100
                // Sometimes the first bases are the same and we can trim them
                int startDiff = Integer.MAX_VALUE;
                for (int i = 0; i < reference.length(); i++)
                    if (reference.charAt(i) != alt.charAt(i)) startDiff = Math.min(startDiff, i);

                // MNPs
                // Sometimes the last bases are the same and we can trim them
                int endDiff = 0;
                for (int i = reference.length() - 1; i >= 0; i--)
                    if (reference.charAt(i) != alt.charAt(i)) endDiff = Math.max(endDiff, i);

                String newRef = reference.substring(startDiff, endDiff + 1);
                String newAlt = alt.substring(startDiff, endDiff + 1);
                variant = Variant.create(chromo, start + startDiff, newRef, newAlt, id);
            }
        } else {
            // Short Insertions, Deletions or Mixed Variants (substitutions)
            VcfRefAltAlign align = new VcfRefAltAlign(alt, reference);
            align.align();
            int startDiff = align.getOffset();

            switch (align.getVariantType()) {
                case DEL:
                    // Case: Deletion
                    // 20     2 .         TC      T      .   PASS  DP=100
                    // 20     2 .         AGAC    AAC    .   PASS  DP=100
                    String ref = "";
                    String ch = align.getAlignment();
                    if (!ch.startsWith("-")) throw new RuntimeException("Deletion '" + ch + "' does not start with '-'. This should never happen!");
                    variant = Variant.create(chromo, start + startDiff, ref, ch, id);
                    break;

                case INS:
                    // Case: Insertion of A { tC ; tCA } tC is the reference allele
                    // 20     2 .         TC      TCA    .   PASS  DP=100
                    ch = align.getAlignment();
                    ref = "";
                    if (!ch.startsWith("+")) throw new RuntimeException("Insertion '" + ch + "' does not start with '+'. This should never happen!");
                    variant = Variant.create(chromo, start + startDiff, ref, ch, id);
                    break;

                case MIXED:
                    // Case: Mixed variant (substitution)
                    reference = reference.substring(startDiff);
                    alt = alt.substring(startDiff);
                    variant = Variant.create(chromo, start + startDiff, reference, alt, id);
                    break;

                default:
                    // Other change type?
                    throw new RuntimeException("Unsupported VCF change type '" + align.getVariantType() + "'\n\tRef: " + reference + "'\n\tAlt: '" + alt + "'\n\tVcfEntry: " + this);
            }
        }

        //---
        // Add original 'ALT' field as genotype
        //---
        if (variant == null) return null;
        variant.setGenotype(alt);

        return variant;
    }

    public boolean isVar(Genotype gt){
        return  gt.isCalled() && !gt.isHomRef();
    }

    public List<String> getAlts() {
        return alts;
    }

    public Set<String> getGenes() {
        Set<String> genes = new HashSet<>();
        if(annotationContexts == null || annotationContexts.isEmpty())
            return null;
        for (AnnotationContext ac : annotationContexts) {
            if (!ac.getGeneName().equals("")) {
                genes.add(ac.getGeneName());
            }
        }
        return genes;
    }

    public Set<String> getTranscriptIds() {
        Set<String> transcriptIds = new HashSet<>();
        if(annotationContexts == null || annotationContexts.isEmpty())
            return null;
        for (AnnotationContext ac : annotationContexts) {
            if (ac.getTranscriptId().startsWith("NR"))
                continue;
            if (!ac.getTranscriptId().equals("")) {
                transcriptIds.add(ac.getTranscriptId());
            }
        }
        return transcriptIds;
    }

    public List<AnnotationContext> getAnnotationContexts() {
        return annotationContexts;
    }

    public List<AnnotationContext> getAnnotationContexts(String alt) {
        return annotationContexts.stream().filter(ac -> ac.getAlt().equals(alt)).collect(Collectors.toList());
    }

    public void setAnnotationContexts(List<AnnotationContext> annotationContexts) {
        this.annotationContexts = annotationContexts;
    }


    public String getContig() {
        return contig;
    }

    public void setContig(String contig) {
        this.contig = contig;
    }

    public String getRefStr() {
        return refStr;
    }

    public void setRefStr(String refStr) {
        this.refStr = refStr;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public void setAlts(List<String> alts) {
        this.alts = alts;
    }

    public Map<String, SampleAnnotationContext> getSampleAnnoContexts() {
        return sampleAnnoContexts;
    }

    public void setSampleAnnoContexts(Map<String, SampleAnnotationContext> sampleAnnoContexts) {
        this.sampleAnnoContexts = sampleAnnoContexts;
    }

    public String getChromeNoChr(){
        if(getContig().startsWith("chr")){
            return getContig().substring(3);
        }
        return getContig();
    }
    public String getChrome(){
        if(!getContig().startsWith("chr")){
            return "chr"+getContig();
        }
        return getContig();
    }

    public String getAnnoStr() {
        return annoStr;
    }

    public Map<String, String> getAnnoItems() {
        return annoItems;
    }

    public String getAnnoItem(String key) {
        if(!annoItems.containsKey(key))
            return ".";
        return annoItems.get(key);
    }

    public void parseAnnotationStrings(String annoLine, List<String> header){
        annoItems = new HashMap<>();
        String[] fields = annoLine.split("\tSI:;");
        annoStr = fields[0];
        String[] annoFields = annoStr.split("\t");
        int i = 1;
        for(String k: header){
            annoItems.put(k,annoFields[i]);
            i ++;
        }

        if(fields.length > 1 && !fields[1].equals("")){
            for(String sampleInfo: fields[1].split(";")){
                SampleAnnotationContext sac = new SampleAnnotationContext();
                sac.parseAlleleString(sampleInfo);
                sampleAnnoContexts.put(sac.getSampleName(), sac);
            }
        }
    }

    public String getFieldByName(String field) {
        switch (field) {
            case "CHROM":
                return contig;

            case "REF":
                return getRefStr();

            case "START":
                return Integer.toString(getStart()-1);

            case "END":
                return Integer.toString(getEnd());

            case "FLKSEQ":
                return getFlankSeq();

            default:
                return ".";
        }
    }

    public List<String> toAnnotationStrings(List<String> fields) {
        List<String> annoStrings = new ArrayList<>();
        for(AnnotationContext ac : annotationContexts){
            StringBuilder sb = new StringBuilder();
            for (String field : fields) {
                if(ArrayUtils.contains(ANN_FIELD_NAMES, field)){
                    sb.append(getFieldByName(field));
                }else {
                    sb.append(ac.getAnnoItemAsString(field, "."));
                }
                sb.append("\t");
            }

            sb.append("SI:");
            for(SampleAnnotationContext sac: sampleAnnoContexts.values()){
                if(sac.hasAlt(ac.getAlt()) || !sac.isCalled()){
                    sb.append(";");
                    sb.append(sac.toAlleleString(ac.getAlt()));
                };
            }
            annoStrings.add(sb.toString());
        }
        return annoStrings;
    }

    public List<Map<String, String>> toAnnotationMaps(List<String> fields) {
        List<Map<String, String>> annoResults = new ArrayList<>();
        for(AnnotationContext ac : annotationContexts){
            Map<String, String> annoResult = new HashMap<>();
            for (String field : fields){
                String value = ac.getAnnoItemAsString(field, ".");
                annoResult.put(field, value);
            }
            annoResult.put(AnnoDBQuery.INDEX_ALT_COLUMN_NAME,ac.getGenotype());
            annoResults.add(annoResult);
        }
        return annoResults;
    }

    public Map<String, List<VariantContext>> toAnnotationVariantContexts(List<String> fields) {
        Map<String, List<VariantContext>> annotationVariantContextMap = new HashMap<>();
        for (Map.Entry<String, List<VariantContext>> entry: variantContextMap.entrySet()){
            String filename = entry.getKey();
            List<VariantContext>  variantContextList = entry.getValue();
            List<VariantContext> annotationVariantContextList = new ArrayList<>();
            for(VariantContext vc : variantContextList){
                VariantContextBuilder variantContextBuilder = new VariantContextBuilder(vc);
                List<String> annoList = new ArrayList<>();
                for(AnnotationContext ac : annotationContexts){
                    Allele allele = ac.getAllele();
                    if(vc.hasAlternateAllele(allele)) {
                        annoList.add(ac.toAnnoString(fields));
                    }
                }
                variantContextBuilder.attribute("ANNO", String.join(",", annoList));
                annotationVariantContextList.add(variantContextBuilder.make());
            }
            annotationVariantContextMap.put(filename, annotationVariantContextList);
        }
        return annotationVariantContextMap;
    }

    public String toString() {
        return getContig() +
                "\t" +  start +
                "\t" +  refStr;
    }

    public void setFlankSeq(String flankSeq) {
        this.flankSeq = flankSeq;
    }

    public String getFlankSeq() {
        return flankSeq;
    }
}
