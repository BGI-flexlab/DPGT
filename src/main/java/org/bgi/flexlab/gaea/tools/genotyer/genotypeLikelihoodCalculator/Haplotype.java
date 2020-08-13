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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator;

import htsjdk.samtools.Cigar;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.ArrayUtils;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.util.ReadUtils;

import java.io.Serializable;
import java.util.*;

/**
 * Created by zhangyong on 2017/2/6.
 * mainly came from GATK 2.3.9-lite
 */
public class Haplotype {
    protected final byte[] bases;
    protected final double[] quals;
    private GenomeLocation genomeLocation = null;
    private HashMap<Integer, VariantContext> eventMap = null;
    private boolean isRef = false;
    private Cigar cigar;
    private int alignmentStartHapwrtRef;
    public int leftBreakPoint = 0;
    public int rightBreakPoint = 0;
    private Allele artificialAllele = null;
    private int artificialAllelePosition = -1;

    /**
     * Create a simple consensus sequence with provided bases and a uniform quality over all bases of qual
     *
     * @param bases bases
     * @param qual  qual
     */
    public Haplotype( final byte[] bases, final int qual ) {
        this.bases = bases.clone();
        quals = new double[bases.length];
        Arrays.fill(quals, (double)qual);
    }

    public Haplotype( final byte[] bases, final double[] quals ) {
        this.bases = bases.clone();
        this.quals = quals.clone();
    }

    public Haplotype( final byte[] bases ) {
        this(bases, 0);
    }

    protected Haplotype( final byte[] bases, final Allele artificialAllele, final int artificialAllelePosition ) {
        this(bases, 0);
        this.artificialAllele = artificialAllele;
        this.artificialAllelePosition = artificialAllelePosition;
    }

    public Haplotype( final byte[] bases, final GenomeLocation loc ) {
        this(bases);
        this.genomeLocation = loc;
    }

    @Override
    public boolean equals( Object h ) {
        return h instanceof Haplotype && Arrays.equals(bases, ((Haplotype) h).bases);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(bases);
    }

    public HashMap<Integer, VariantContext> getEventMap() {
        return eventMap;
    }

    public void setEventMap( final HashMap<Integer, VariantContext> eventMap ) {
        this.eventMap = eventMap;
    }

    public boolean isReference() {
        return isRef;
    }

    public void setIsReference( boolean isRef ) {
        this.isRef = isRef;
    }

    public double getQualitySum() {
        double s = 0;
        for (int k=0; k < bases.length; k++) {
            s += quals[k];
        }
        return s;
    }

    @Override
    public String toString() {
        return new String(bases);
    }

    public double[] getQuals() {
        return quals.clone();
    }
    public byte[] getBases() {
        return bases.clone();
    }

    public long getStartPosition() {
        return genomeLocation.getStart();
    }

    public long getStopPosition() {
        return genomeLocation.getStop();
    }

    public int getAlignmentStartHapwrtRef() {
        return alignmentStartHapwrtRef;
    }

    public void setAlignmentStartHapwrtRef( final int alignmentStartHapwrtRef ) {
        this.alignmentStartHapwrtRef = alignmentStartHapwrtRef;
    }

    public Cigar getCigar() {
        return cigar;
    }

    public void setCigar( final Cigar cigar ) {
        this.cigar = cigar;
    }

    public boolean isArtificialHaplotype() {
        return artificialAllele != null;
    }

    public Allele getArtificialAllele() {
        return artificialAllele;
    }

    public int getArtificialAllelePosition() {
        return artificialAllelePosition;
    }

    public void setArtificialAllele(final Allele artificialAllele, final int artificialAllelePosition) {
        this.artificialAllele = artificialAllele;
        this.artificialAllelePosition = artificialAllelePosition;
    }

    //@Requires({"refInsertLocation >= 0"})
    public Haplotype insertAllele( final Allele refAllele, final Allele altAllele, final int refInsertLocation, final int genomicInsertLocation ) {
        // refInsertLocation is in ref haplotype offset coordinates NOT genomic coordinates
        final int haplotypeInsertLocation = ReadUtils.getReadCoordinateForReferenceCoordinate(alignmentStartHapwrtRef, cigar, refInsertLocation, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if( haplotypeInsertLocation == -1 || haplotypeInsertLocation + refAllele.length() >= bases.length ) { // desired change falls inside deletion so don't bother creating a new haplotype
            return null;
        }
        byte[] newHaplotypeBases = new byte[]{};
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, ArrayUtils.subarray(bases, 0, haplotypeInsertLocation)); // bases before the variant
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, altAllele.getBases()); // the alt allele of the variant
        newHaplotypeBases = ArrayUtils.addAll(newHaplotypeBases, ArrayUtils.subarray(bases, haplotypeInsertLocation + refAllele.length(), bases.length)); // bases after the variant
        return new Haplotype(newHaplotypeBases, altAllele, genomicInsertLocation);
    }

    public static class HaplotypeBaseComparator implements Comparator<Haplotype>, Serializable {
        @Override
        public int compare( final Haplotype hap1, final Haplotype hap2 ) {
            final byte[] arr1 = hap1.getBases();
            final byte[] arr2 = hap2.getBases();
            // compares byte arrays using lexical ordering
            final int len = Math.min(arr1.length, arr2.length);
            for( int iii = 0; iii < len; iii++ ) {
                final int cmp = arr1[iii] - arr2[iii];
                if (cmp != 0) { return cmp; }
            }
            return arr2.length - arr1.length;
        }
    }

    public static LinkedHashMap<Allele,Haplotype> makeHaplotypeListFromAlleles(final List<Allele> alleleList,
                                                                               final int startPos,
                                                                               final ChromosomeInformationShare ref,
                                                                               final GenomeLocation refWindows,
                                                                               final GenomeLocationParser locationParser,
                                                                               final int haplotypeSize,
                                                                               final int numPrefBases) {
        LinkedHashMap<Allele,Haplotype> haplotypeMap = new LinkedHashMap<Allele, Haplotype>();

        Allele refAllele = null;

        for (Allele a:alleleList) {
            if (a.isReference()) {
                refAllele = a;
                break;
            }
        }

        if (refAllele == null)
            throw new RuntimeException("BUG: no ref alleles in input to makeHaplotypeListfrom Alleles at loc: "+ startPos);

        //byte[] refBases = ref.getBases();

        final int startIdxInReference = 1 + startPos - numPrefBases;
        final String basesBeforeVariant = ref.getGA4GHBaseSequence(startIdxInReference, startIdxInReference + numPrefBases - 1);//new String(Arrays.copyOfRange(refBases, startIdxInReference, startIdxInReference + numPrefBases));

        // protect against long events that overrun available reference context
        final int startAfter = Math.min(startIdxInReference + numPrefBases + refAllele.getBases().length - 1, refWindows.getStop());
        final String basesAfterVariant = ref.getGA4GHBaseSequence(startAfter, refWindows.getStop());//new String(Arrays.copyOfRange(refBases, startAfter, refBases.length));

        //System.err.println("start index in reference:" + startIdxInReference + "\tstart after:" + startAfter);
        //System.err.println("before bases:" + basesBeforeVariant);
        //System.err.println("after bases:" + basesAfterVariant);

        // Create location for all haplotypes
        final int startLoc = startIdxInReference;
        final int stopLoc = startLoc + haplotypeSize - 1;

        final GenomeLocation locus = locationParser.createGenomeLocation(ref.getChromosomeName(),startLoc,stopLoc);

        for (final Allele a : alleleList) {
            //System.err.println("mk hp allele:" + a.getBaseString());

            final byte[] alleleBases = a.getBases();
            if(alleleBases.length <= 0)
                continue;
            // use string concatenation
            String haplotypeString = basesBeforeVariant + new String(Arrays.copyOfRange(alleleBases, 1, alleleBases.length)) + basesAfterVariant;
            haplotypeString = haplotypeString.substring(0,haplotypeSize);
            //System.err.println("haplotype bases:" + haplotypeString);

            haplotypeMap.put(a,new Haplotype(haplotypeString.getBytes(), locus));
        }

        return haplotypeMap;
    }
}
