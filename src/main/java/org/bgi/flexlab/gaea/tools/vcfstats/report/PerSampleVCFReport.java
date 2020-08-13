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
package org.bgi.flexlab.gaea.tools.vcfstats.report;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.vcfstats.VariantType;
import org.bgi.flexlab.gaea.util.Histogram;
import org.bgi.flexlab.gaea.util.Pair;
import org.bgi.flexlab.gaea.util.StatsUtils;

import java.util.ArrayList;
import java.util.List;

public class PerSampleVCFReport {

    private static final String[] VARIANT_TYPE_COUNT_LENGTH = {
            "SNP", "MNP", "INS", "DEL"
    };

    private static final String ALLELE_LENGTH_TAG = "AlleleLen:";
    private static final int ALLELE_LENGTH_MAX_PRINT = 200;
    private ReferenceShare genomeShare;
    private boolean countVarLength = false;

    private String sampleName;

    private long mTotalUnchanged = 0;
    private long mHeterozygous = 0;
    private long mHomozygous = 0;
    private long mNoCall = 0;
    private long mDeNovo = 0;
    private long mPhased = 0;

    private long mTotalSnps = 0;
    private long mTransitions = 0;
    private long mTransversions = 0;
    private long mHeterozygousSnps = 0;
    private long mHomozygousSnps = 0;

    private long mTotalMnps = 0;
    private long mHeterozygousMnps = 0;
    private long mHomozygousMnps = 0;

    private long mTotalMixeds = 0;
    private long mHeterozygousMixeds = 0;
    private long mHomozygousMixeds = 0;

    private long mTotalInsertions = 0;
    private long mHeterozygousInsertions = 0;
    private long mHomozygousInsertions = 0;

    private long mTotalDeletions = 0;
    private long mHeterozygousDeletions = 0;
    private long mHomozygousDeletions = 0;

    private long mTotalBreakends = 0;
    private long mHeterozygousBreakends = 0;
    private long mHomozygousBreakends = 0;

    private long mTotalFailedFilters = 0;
    private long mTotalPassedFilters = 0;

    private long mTotalDbSnp = 0;

    private long mTotalInDels = 0;
    private long mHeterozygousInDels = 0;
    private long mHomozygousInDels = 0;

    private long mKnownTransitions = 0;
    private long mKnownTransversions = 0;
    private long mNovelTransitions = 0;
    private long mNovelTransversions = 0;


    private Histogram[] mAlleleLengths = null;


    public PerSampleVCFReport(ReferenceShare genomeShare, boolean countVarLength){
        sampleName = null;
        this.genomeShare = genomeShare;
        this.countVarLength = countVarLength;

        if(countVarLength) {
            mAlleleLengths = new Histogram[VARIANT_TYPE_COUNT_LENGTH.length];
            for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
                // i from SNP as we don't care about NO_CALL/UNCHANGED
                mAlleleLengths[i] = new Histogram();
            }
        }
    }

    // reducerStr: without sampleName
    public void parseReducerString(String reducerStr){

        if(countVarLength && reducerStr.startsWith(ALLELE_LENGTH_TAG)){
            String[] fields = reducerStr.split("\t", 3);
            VariantType type = VariantType.valueOf(fields[1]);
            if(type.ordinal() < mAlleleLengths.length)
                mAlleleLengths[type.ordinal()].addHistogram(fields[2]);
            return;
        }

        String[] fields = reducerStr.split("\t");
        int i = 0;

        mTotalUnchanged += Integer.valueOf(fields[i++]);
        mHeterozygous += Integer.valueOf(fields[i++]);
        mHomozygous += Integer.valueOf(fields[i++]);
        mNoCall += Integer.valueOf(fields[i++]);
        mDeNovo += Integer.valueOf(fields[i++]);
        mPhased += Integer.valueOf(fields[i++]);

        mTotalSnps += Integer.valueOf(fields[i++]);
        mTransitions += Integer.valueOf(fields[i++]);
        mTransversions += Integer.valueOf(fields[i++]);
        mKnownTransitions += Integer.valueOf(fields[i++]);
        mKnownTransversions += Integer.valueOf(fields[i++]);
        mNovelTransitions += Integer.valueOf(fields[i++]);
        mNovelTransversions += Integer.valueOf(fields[i++]);
        mHeterozygousSnps += Integer.valueOf(fields[i++]);
        mHomozygousSnps += Integer.valueOf(fields[i++]);

        mTotalMnps += Integer.valueOf(fields[i++]);
        mHeterozygousMnps += Integer.valueOf(fields[i++]);
        mHomozygousMnps += Integer.valueOf(fields[i++]);

        mTotalMixeds += Integer.valueOf(fields[i++]);
        mHeterozygousMixeds += Integer.valueOf(fields[i++]);
        mHomozygousMixeds += Integer.valueOf(fields[i++]);

        mTotalInsertions += Integer.valueOf(fields[i++]);
        mHeterozygousInsertions += Integer.valueOf(fields[i++]);
        mHomozygousInsertions += Integer.valueOf(fields[i++]);

        mTotalDeletions += Integer.valueOf(fields[i++]);
        mHeterozygousDeletions += Integer.valueOf(fields[i++]);
        mHomozygousDeletions += Integer.valueOf(fields[i++]);

        mTotalBreakends += Integer.valueOf(fields[i++]);
        mHeterozygousBreakends += Integer.valueOf(fields[i++]);
        mHomozygousBreakends += Integer.valueOf(fields[i++]);

        mTotalFailedFilters += Integer.valueOf(fields[i++]);
        mTotalPassedFilters += Integer.valueOf(fields[i++]);

        mTotalDbSnp += Integer.valueOf(fields[i++]);

        mTotalInDels += Integer.valueOf(fields[i++]);
        mHeterozygousInDels += Integer.valueOf(fields[i++]);
        mHomozygousInDels += Integer.valueOf(fields[i]);
    }

    public String toReducerString(){
        StringBuilder sb = new StringBuilder();
        sb.append(sampleName);
        sb.append("\t");
        String value = String.join("\t", getStatistics());
        sb.append(value);
        sb.append("\n");

        if(!countVarLength)
            return sb.toString();

        for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
            Histogram histogram = mAlleleLengths[i];
            sb.append(sampleName);
            sb.append("\t");
            sb.append(ALLELE_LENGTH_TAG);
            sb.append("\t");
            sb.append(VARIANT_TYPE_COUNT_LENGTH[i]);
            sb.append("\t");
            sb.append(histogram.toString());
            sb.append("\n");
        }

        return sb.toString();
    }

    public void add(VariantContext vc, String sample) {
        setSampleName(sample);
        Genotype gt = vc.getGenotype(sample);
        VariantType type = VariantType.determineType(vc, sample);

        if(vc.isFiltered()) {
            mTotalFailedFilters++;
            return;
        }else
            mTotalPassedFilters++;

        if(gt.isNoCall()) {
            mNoCall++;
            return;
        }

        if(gt.isHomRef()) {
            mTotalUnchanged++;
            return;
        }

        if (gt.isHet())
            mHeterozygous++;
        else
            mHomozygous++;

//        if(genomeShare.getChromosomeInfo(vc.getContig()).isSNP(vc.getStart()))
//            mTotalDbSnp++;

        if(vc.hasID())
            mTotalDbSnp++;

        for(Allele allele: gt.getAlleles()){
            if(allele.isReference() || Allele.wouldBeStarAllele(allele.getBases()))
                continue;
            if(countVarLength)
                tallyAlleleLengths(vc.getReference(), allele);
            if(type == VariantType.SNP){
                tallyTransitionTransversionRatio(vc.getReference().getBaseString(), allele.getBaseString(), vc.hasID());
            }
            if(type == VariantType.MIXED){
                VariantType type2 =VariantType.typeOfBiallelicVariant(vc.getReference(), allele);
                count(type2, gt.isHet());
            }
        }

        count(type, gt.isHet());
    }

    private void count(VariantType type, boolean isHet){
        switch (type) {
            case SNP:
                if (isHet) mHeterozygousSnps++;
                else mHomozygousSnps++;
                mTotalSnps++;
                break;
            case MNP:
                if (isHet) mHeterozygousMnps++;
                else mHomozygousMnps++;
                mTotalMnps++;
                break;
            case MIXED:
                if (isHet) mHeterozygousMixeds++;
                else mHomozygousMixeds++;
                mTotalMixeds++;
                break;
            case INS:
                if (isHet) mHeterozygousInsertions++;
                else mHomozygousInsertions++;
                mTotalInsertions++;
                break;
            case DEL:
                if (isHet) mHeterozygousDeletions++;
                else mHomozygousDeletions++;
                mTotalDeletions++;
                break;
            case InDel:
                if (isHet) mHeterozygousInDels++;
                else mHomozygousInDels++;
                mTotalInDels++;
                break;
            case BND:
                if (isHet) mHeterozygousBreakends++;
                else mHomozygousBreakends++;
                mTotalBreakends++;
                break;
            default:
                break;
        }
    }


    private void tallyAlleleLengths(Allele ref, Allele alt) {
        VariantType type = VariantType.typeOfBiallelicVariant(ref, alt);
        int len = Math.max(ref.length(), alt.length());
        if(type.ordinal() < mAlleleLengths.length)
            mAlleleLengths[type.ordinal()].increment(len, 1);
    }

    List<String> getStatistics() {
        final List<String> values = new ArrayList<>();

        values.add(Long.toString(mTotalUnchanged));
        values.add(Long.toString(mHeterozygous));
        values.add(Long.toString(mHomozygous));
        values.add(Long.toString(mNoCall));
        values.add(Long.toString(mDeNovo));
        values.add(Long.toString(mPhased));

        values.add(Long.toString(mTotalSnps));
        values.add(Long.toString(mTransitions));
        values.add(Long.toString(mTransversions));
        values.add(Long.toString(mKnownTransitions));
        values.add(Long.toString(mKnownTransversions));
        values.add(Long.toString(mNovelTransitions));
        values.add(Long.toString(mNovelTransversions));
        values.add(Long.toString(mHeterozygousSnps));
        values.add(Long.toString(mHomozygousSnps));

        values.add(Long.toString(mTotalMnps));
        values.add(Long.toString(mHeterozygousMnps));
        values.add(Long.toString(mHomozygousMnps));

        values.add(Long.toString(mTotalMixeds));
        values.add(Long.toString(mHeterozygousMixeds));
        values.add(Long.toString(mHomozygousMixeds));

        values.add(Long.toString(mTotalInsertions));
        values.add(Long.toString(mHeterozygousInsertions));
        values.add(Long.toString(mHomozygousInsertions));

        values.add(Long.toString(mTotalDeletions));
        values.add(Long.toString(mHeterozygousDeletions));
        values.add(Long.toString(mHomozygousDeletions));

        values.add(Long.toString(mTotalBreakends));
        values.add(Long.toString(mHeterozygousBreakends));
        values.add(Long.toString(mHomozygousBreakends));

        values.add(Long.toString(mTotalFailedFilters));
        values.add(Long.toString(mTotalPassedFilters));

        values.add(Long.toString(mTotalDbSnp));

        values.add(Long.toString(mTotalInDels));
        values.add(Long.toString(mHeterozygousInDels));
        values.add(Long.toString(mHomozygousInDels));

        return values;
    }

    private Pair<List<String>, List<String>> getReportResult() {
        final List<String> names = new ArrayList<>();
        final List<String> values = new ArrayList<>();
        names.add("SampleName");
        values.add(sampleName);
        names.add("Failed Filters");
        values.add(Long.toString(mTotalFailedFilters));
        names.add("Passed Filters");
        values.add(Long.toString(mTotalPassedFilters));
        names.add("SNPs");
        values.add(Long.toString(mTotalSnps));
        names.add("MNPs");
        values.add(Long.toString(mTotalMnps));
        names.add("Insertions");
        values.add(Long.toString(mTotalInsertions));
        names.add("Deletions");
        values.add(Long.toString(mTotalDeletions));
        names.add("InDels");
        values.add(Long.toString(mTotalInDels));
        names.add("MIXEDs");
        values.add(Long.toString(mTotalMixeds));
//        names.add("Structural variant breakends");
//        values.add(mTotalBreakends > 0 ? Long.toString(mTotalBreakends) : "-");
        names.add("Same as reference");
        values.add(Long.toString(mTotalUnchanged));
        names.add("Missing Genotype");
        values.add(Long.toString(mNoCall));
//        names.add("De Novo Genotypes");
//        values.add(mDeNovo > 0 ? Long.toString(mDeNovo) : null);
//        names.add("Phased Genotypes");
//        final long totalNonMissingGenotypes = mTotalSnps + mTotalMnps + mTotalInsertions + mTotalDeletions + mTotalUnchanged;
//        values.add(mPhased > 0 ? StatsUtils.percent(mPhased, totalNonMissingGenotypes) : null);
        names.add("SNP Transitions/Transversions");
        values.add(StatsUtils.divide(mTransitions, mTransversions));
        names.add("Known SNP Transitions/Transversions");
        values.add(StatsUtils.divide(mKnownTransitions, mKnownTransversions));
        names.add("Novel SNP Transitions/Transversions");
        values.add(StatsUtils.divide(mNovelTransitions, mNovelTransversions));

        names.add("Total Het/Hom ratio");
        values.add(StatsUtils.divide(mHeterozygous, mHomozygous));
        names.add("SNP Het/Hom ratio");
        values.add(StatsUtils.divide(mHeterozygousSnps, mHomozygousSnps));
        names.add("MNP Het/Hom ratio");
        values.add(StatsUtils.divide(mHeterozygousMnps, mHomozygousMnps));
        names.add("Insertion Het/Hom ratio");
        values.add(StatsUtils.divide(mHeterozygousInsertions, mHomozygousInsertions));
        names.add("Deletion Het/Hom ratio");
        values.add(StatsUtils.divide(mHeterozygousDeletions, mHomozygousDeletions));

        names.add("Insertion/Deletion ratio");
        values.add(StatsUtils.divide(mTotalInsertions, mTotalDeletions));
        names.add("Indel/SNP+MNP ratio");
        values.add(StatsUtils.divide(mTotalInsertions + mTotalDeletions, mTotalSnps + mTotalMnps));
        names.add("dbSNP ratio");
        values.add(StatsUtils.divide(mTotalDbSnp, mTotalPassedFilters - mNoCall - mTotalUnchanged));
        return Pair.create(names, values);
    }

    private void tallyTransitionTransversionRatio(String ref, String pred) {
        final boolean transition = "AG".contains(ref) && "AG".contains(pred) || "CT".contains(ref) && "CT".contains(pred);
        if (transition) {
            mTransitions++;
        } else {
            mTransversions++;
        }
    }

    private void tallyTransitionTransversionRatio(String ref, String pred, boolean hasID) {
        final boolean transition = "AG".contains(ref) && "AG".contains(pred) || "CT".contains(ref) && "CT".contains(pred);
        if (transition) {
            mTransitions++;
            if(hasID)
                mKnownTransitions++;
            else
                mNovelTransitions++;
        } else {
            mTransversions++;
            if(hasID)
                mKnownTransversions++;
            else
                mNovelTransversions++;
        }
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    /**
     * Append per sample histograms to a buffer.
     * @param sb buffer to append to
     * TODO snp长度个数统计有误？
     */
    private void appendHistograms(StringBuilder sb) {
        sb.append("Variant Allele Lengths :").append("\n");
        //sb.append("bin\tSNP\tMNP\tInsert\tDelete\tIndel").append(StringUtils.LS);
        sb.append("length");
        for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
            if (i <= VariantType.DEL.ordinal() || mAlleleLengths[i].getLength() != 0) {
                sb.append("\t").append(VARIANT_TYPE_COUNT_LENGTH[i]);
            }
        }
        sb.append("\n");

        int size = 0;
        final int[] lengths = new int[mAlleleLengths.length];
        for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
            lengths[i] = mAlleleLengths[i].getLength();
            if (lengths[i] > size) {
                size = lengths[i];
            }
        }
        int bin = 1;
        while (bin < size) {
            if(bin < ALLELE_LENGTH_MAX_PRINT) {
                sb.append(bin);
                for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
                    if (i <= VariantType.DEL.ordinal() || mAlleleLengths[i].getLength() != 0) {
                        long sum = 0L;
                        if(bin<mAlleleLengths[i].getLength())
                            sum = mAlleleLengths[i].getValue(bin);
                        sb.append("\t").append(sum);
                    }
                }
                sb.append("\n");
            }else {
                sb.append(">").append(ALLELE_LENGTH_MAX_PRINT);
                for (int i = VariantType.SNP.ordinal(); i < mAlleleLengths.length; ++i) {
                    if (i <= VariantType.DEL.ordinal() || mAlleleLengths[i].getLength() != 0) {
                        long sum = 0L;
                        for (int j = bin; j < size && j<mAlleleLengths[i].getLength(); ++j) {
                            if (j < lengths[i]) {
                                sum += mAlleleLengths[i].getValue(j);
                            }
                        }
                        sb.append("\t").append(sum);
                    }
                }
                sb.append("\n");
                break;
            }
            bin++;
        }
    }

    public String getReport() {
        Pair<List<String>, List<String>> pair = getReportResult();
        List<String> names = pair.getFirst();
        List<String> values = pair.getSecond();
        StringBuilder outString = new StringBuilder();

        for (int i = 0; i < names.size(); i++) {
            outString.append(String.format("%-37s", names.get(i)));
            outString.append(":  ");
            outString.append(values.get(i));
            outString.append("\n");
        }

        if(!countVarLength)
            return outString.toString();

        outString.append("\n");
        appendHistograms(outString);

        return outString.toString();
    }
}
