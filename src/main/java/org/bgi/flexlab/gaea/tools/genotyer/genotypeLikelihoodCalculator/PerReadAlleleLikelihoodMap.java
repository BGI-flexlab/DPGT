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

import htsjdk.variant.variantcontext.Allele;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;

import java.lang.reflect.Constructor;
import java.util.*;

public  class PerReadAlleleLikelihoodMap {

    public static final double INFORMATIVE_LIKELIHOOD_THRESHOLD = 0.1;
    private Allele refAllele;
    protected List<Allele> alleles;
    protected Map<AlignmentsBasic, Map<Allele, Double>> likelihoodReadMap;

    public PerReadAlleleLikelihoodMap() {
        likelihoodReadMap = new LinkedHashMap<AlignmentsBasic, Map<Allele,Double>>();
        alleles = new ArrayList<Allele>();
    }

    // not implemented in the standard version
    public void performPerAlleleDownsampling(final double downsamplingFraction) {}
    //public abstract ReadBackedPileup createPerAlleleDownsampledBasePileup(final ReadBackedPileup pileup, final double downsamplingFraction, final PrintStream log);

    public void add(AlignmentsBasic read, Allele a, Double likelihood) {
        Map<Allele,Double> likelihoodMap;
        if (likelihoodReadMap.containsKey(read)){
            // seen pileup element before
            likelihoodMap = likelihoodReadMap.get(read);
        }
        else {
            likelihoodMap = new HashMap<Allele, Double>();
            likelihoodReadMap.put(read,likelihoodMap);
        }
        likelihoodMap.put(a,likelihood);

        if (!alleles.contains(a))
            alleles.add(a);

    }

    public int size() {
        return likelihoodReadMap.size();
    }

    public void add(PileupReadInfo p, Allele a, Double likelihood) {
        add(p.getReadInfo(), a, likelihood);
    }

    public boolean containsPileupElement(PileupReadInfo p) {
        return likelihoodReadMap.containsKey(p.getReadInfo());
    }

    public boolean isEmpty() {
        return likelihoodReadMap.isEmpty();
    }

    public Map<AlignmentsBasic,Map<Allele,Double>> getLikelihoodReadMap() {
        return likelihoodReadMap;
    }
    public void clear() {
        alleles.clear();
        likelihoodReadMap.clear();
    }

    public Set<AlignmentsBasic> getStoredElements() {
        return likelihoodReadMap.keySet();
    }

    public Collection<Map<Allele,Double>> getLikelihoodMapValues() {
        return likelihoodReadMap.values();
    }

    public int getNumberOfStoredElements() {
        return likelihoodReadMap.size();
    }

    public Map<Allele,Double> getLikelihoodsAssociatedWithPileupElement(PileupReadInfo p) {
        if (!likelihoodReadMap.containsKey(p.getReadInfo()))
            return null;

        return likelihoodReadMap.get(p.getReadInfo());
    }

    public static Allele getMostLikelyAllele( final Map<Allele,Double> alleleMap ) {
        double maxLike = Double.NEGATIVE_INFINITY;
        double prevMaxLike = Double.NEGATIVE_INFINITY;
        Allele mostLikelyAllele = Allele.NO_CALL;

        for (final Map.Entry<Allele,Double> el : alleleMap.entrySet()) {
            if (el.getValue() > maxLike) {
                prevMaxLike = maxLike;
                maxLike = el.getValue();
                mostLikelyAllele = el.getKey();
            } else if( el.getValue() > prevMaxLike ) {
                prevMaxLike = el.getValue();
            }
        }
        return (maxLike - prevMaxLike > INFORMATIVE_LIKELIHOOD_THRESHOLD ? mostLikelyAllele : Allele.NO_CALL );
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
	public static PerReadAlleleLikelihoodMap getBestAvailablePerReadAlleleLikelihoodMap() {
        //final Class PerReadAlleleLikelihoodMapClass = GATKLiteUtils.getProtectedClassIfAvailable(PerReadAlleleLikelihoodMap.class);
    	final Class PerReadAlleleLikelihoodMapClass = PerReadAlleleLikelihoodMap.class;
    	try {
            Constructor constructor = PerReadAlleleLikelihoodMapClass.getDeclaredConstructor((Class[])null);
            constructor.setAccessible(true);
            return (PerReadAlleleLikelihoodMap)constructor.newInstance();
        }
        catch (Exception e) {
            throw new UserException("Unable to create RecalibrationEngine class instance " + PerReadAlleleLikelihoodMapClass.getSimpleName());
        }
    }
	/**
	 * @return the refAllele
	 */
	public Allele getRefAllele() {
		return refAllele;
	}
	/**
	 * @param refAllele the refAllele to set
	 */
	public void setRefAllele(Allele refAllele) {
		this.refAllele = refAllele;
	}
}
