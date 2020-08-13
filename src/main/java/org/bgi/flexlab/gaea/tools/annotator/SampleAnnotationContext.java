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
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *  The special annotation context for each sample
 */
public class SampleAnnotationContext{

	public enum FilterTag {
		FAIL, PASS, DUBIOUS
	}

	private String sampleName;
	private List<String> alts;   // variants at chr:pos
	private int depth;
	private int alleleDepthSum;
	private Map<String, Integer> alleleDepths = null;
	private Map<String, String> alleleRatios = null;
	private Map<String, String> zygosity = null;
	private Map<String, String> filter = null;
	private boolean hasNearVar = false;
	private boolean isCalled;
	private String singleAlt;
	private String qual;

	public SampleAnnotationContext() {}

	public SampleAnnotationContext(String sampleName) {
		this.sampleName = sampleName;
	}

	public SampleAnnotationContext(String sampleName, VariantContext variantContext) {
		this.sampleName = sampleName;
		init(variantContext);
	}

	void init(VariantContext variantContext){
		DecimalFormat df = new DecimalFormat("0.000");
		df.setRoundingMode(RoundingMode.HALF_UP);
		Genotype gt = variantContext.getGenotype(sampleName);
		Map<String, Integer> alleleDepths = new HashMap<>();
		Map<String, String> zygosity = new HashMap<>();
		int i = 0;
		List<String> alts = new ArrayList<>();
		for(Allele allele: variantContext.getAlleles()){
			if(gt.hasAD()){
				alleleDepths.put(allele.getBaseString(), gt.getAD()[i]);
			}
			zygosity.put(allele.getBaseString(), getZygosityType(gt));
			if(i > 0)
				alts.add(allele.getBaseString());
			i++;
		}
		setQual(df.format(variantContext.getPhredScaledQual()));
		setCalled(gt.isCalled());
		setAlleleDepths(alleleDepths);
		setDepth(gt.getDP());
		setAlts(alts);
		setAlleleDepthSum(gt);
		setZygosity(zygosity);
		setFilter(variantContext, gt);
	}

	void add(VariantContext variantContext){
		Genotype gt = variantContext.getGenotype(sampleName);
		int[] AD = gt.getAD();
		int i = 0;
		for(Allele allele: variantContext.getAlternateAlleles()){
			i++;
			if(alts.contains(allele.getBaseString()))
				continue;
			alts.add(allele.getBaseString());
			alleleDepths.put(allele.getBaseString(), AD[i]);
			zygosity.put(allele.getBaseString(), getZygosityType(gt));
			filter.put(allele.getBaseString(), calculateAlleleFilterTag(variantContext, gt, AD[i]).toString());
		}
	}

	private String getZygosityType(Genotype gt){
		if(gt.isNoCall())
			return "noCall";
		if(gt.isHetNonRef())
			return "het-alt";
		if(gt.isHet())
			return "het-ref";
		if(gt.isHomVar())
			return "hom-alt";
		return ".";
	}

	public String getFieldByName(String fieldName, String allele) {
		switch (fieldName) {
			case "SAMPLE":
				return getSampleName();

			case "NbGID":
				if(hasNearVar) return "1";
				else return "0";

			case "A.Depth":
				return Integer.toString(getAlleleDepth(allele));

			case "A.Ratio":
				return getAlleleRatio(allele);

			case "Zygosity":
				return getAlleleZygosity(allele);

			case "Filter":
				return getAlleleFilter(allele);

			case "QUAL":
				return getQual();

			case "DP":
				return Integer.toString(getDepth());

			default:
				return null;
		}
	}

	public String getSampleName() {
		return sampleName;
	}

	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}

	public boolean isHasNearVar() {
		return hasNearVar;
	}

	public void setHasNearVar() {
		this.hasNearVar = true;
		updateFilterTag();
	}

	private void updateFilterTag() {
		for (String key: filter.keySet()){
			if(filter.get(key).equals(FilterTag.PASS.toString()))
				filter.put(key, FilterTag.DUBIOUS.toString());
		}
	}

	private void setAlleleRatios(){
		if(alleleDepths.isEmpty() || getDepth() == -1) return;

		DecimalFormat df = new DecimalFormat("0.00");
		df.setRoundingMode(RoundingMode.HALF_UP);
		alleleRatios = new HashMap<>();
		for(String alt: getAlts()){
			double ratio = getDepth() == 0 ? 0 : getAlleleDepth(alt)*1.0 / getAlleleDepthSum();
			alleleRatios.put(alt, df.format(ratio));
		}
	}

	public String getAlleleZygosity(String allele) {
		return zygosity.get(allele);
	}

	public Map<String, String> getZygosity() {
		return zygosity;
	}

	public void setZygosity(Map<String, String> zygosity) {
		this.zygosity = zygosity;
	}

	public void setAlleleDepthSum(Genotype gt) {
		alleleDepthSum = 0;
		for(int altDP: gt.getAD()){
			alleleDepthSum += altDP;
		}
	}

	public int getDepth() {
		return depth;
	}

	public int getAlleleDepthSum() {
		return alleleDepthSum;
	}

	public void setDepth(int depth) {
		this.depth = depth;
	}

	public Map<String, Integer> getAlleleDepths() {
		return alleleDepths;
	}

	public void setAlleleDepths(Map<String, Integer> alleleDepths) {
		this.alleleDepths = alleleDepths;
	}

	public List<String> getAlts() {
		return alts;
	}

	public void setAlts(List<String> alts) {
		this.alts = alts;
	}

	public String getSingleAlt() {
		return singleAlt;
	}

	public boolean hasAlt(String alt){
		return alts.contains(alt);
	}

	public String getAlleleRatio(String allele){
		if(alleleRatios == null || alleleRatios.isEmpty())
			return ".";
		return alleleRatios.get(allele);
	}

	public int getAlleleDepth(String allele){
		if(alleleDepths == null || alleleDepths.isEmpty())
			return -1;
		return alleleDepths.get(allele);
	}

	public String toAlleleString(String allele){
		if(null == alleleRatios)
			setAlleleRatios();
		StringBuilder sb = new StringBuilder();
		sb.append(getSampleName());
		sb.append("|");
		sb.append(allele);
		sb.append("|");
		sb.append(getAlleleRatio(allele));
		sb.append("|");
		sb.append(getAlleleDepth(allele) == -1 ? "." : getAlleleDepth(allele));
		sb.append("|");
		sb.append(isHasNearVar() ? 1 : 0);
		sb.append("|");
		sb.append(getAlleleZygosity(allele));
		sb.append("|");
		sb.append(getAlleleFilter(allele));
		sb.append("|");
		sb.append(getQual());
		sb.append("|");
		sb.append(getDepth());
		return sb.toString();
	}

	public void parseAlleleString(String alleleString){
		String[] fields = alleleString.split("\\|");
		setSampleName(fields[0]);
		singleAlt = fields[1];
		alleleRatios = new HashMap<>();
		alleleDepths = new HashMap<>();
		zygosity = new HashMap<>();
		filter = new HashMap<>();
		if(!fields[2].equals("."))
			alleleRatios.put(singleAlt, fields[2]);
		if(!fields[3].equals("."))
			alleleDepths.put(singleAlt, Integer.parseInt(fields[3]));
		if(fields[4].equals("1"))
			setHasNearVar();
		zygosity.put(singleAlt, fields[5]);
		filter.put(singleAlt, fields[6]);
		setQual(fields[7]);
		setDepth(Integer.parseInt(fields[8]));
	}

	public boolean isCalled() {
		return isCalled;
	}

	public void setCalled(boolean called) {
		isCalled = called;
	}

	public String getAlleleFilter(String allele) {
		return filter.get(allele);
	}

	public void setFilter(VariantContext vc, Genotype gt) {
		Map<String, String> filter = new HashMap<>();
		for(Allele allele : gt.getAlleles()){
			if(allele.isNonReference()){
				int alleleDepth = getAlleleDepth(allele.getBaseString());
				FilterTag ft = calculateAlleleFilterTag(vc, gt, alleleDepth);
				filter.put(allele.getBaseString(), ft.toString());
			}
		}
		this.filter = filter;
	}

	private FilterTag calculateAlleleFilterTag(VariantContext vc, Genotype gt, int alleleDepth) {
		int plIndicator = getAllelePLindicator(vc, gt);
		if(alleleDepth <= 2 || plIndicator < 0)
			return FilterTag.FAIL;
		if(alleleDepth >= 8 && plIndicator > 0 && vc.isNotFiltered())
			return FilterTag.PASS;
		return FilterTag.DUBIOUS;
	}

	public static int getAllelePLindicator(VariantContext vc, Genotype gt){
		if(!gt.hasPL())
			return 1;
		List<Integer> index = vc.getAlleleIndices(gt.getAlleles());
		int plIndex = GenotypeLikelihoods.calculatePLindex(index.get(0), index.get(1));
		int[] pls = gt.getPL();
		int plValue = pls[plIndex];
		if(plValue > 3 || hasOtherZeroPL(pls, plIndex)){
			return -1;
		}else if(plValue > 0 || hasOtherLessTenPL(pls, plIndex)) {
			return 0;
		}
		return 1;
	}

	public static boolean hasOtherZeroPL(int[] pls, int skipIndex){
		for (int i = 0; i < pls.length; i++) {
			if(skipIndex == i)
				continue;
			if(pls[i] == 0)
				return true;
		}
		return false;
	}

	public static boolean hasOtherLessTenPL(int[] pls, int skipIndex){
		for (int i = 0; i < pls.length; i++) {
			if(skipIndex == i)
				continue;
			if(pls[i] < 10)
				return true;
		}
		return false;
	}

	public String getQual() {
		return qual;
	}

	public void setQual(String qual) {
		this.qual = qual;
	}
}
