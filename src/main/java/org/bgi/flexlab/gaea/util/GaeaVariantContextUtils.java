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

package org.bgi.flexlab.gaea.util;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.*;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.header.VCFConstants;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;

import java.io.Serializable;
import java.util.*;

public class GaeaVariantContextUtils {
    /// if sum(gl) is bigger than this threshold, we treat GL's as non-informative and will force a no-call.
    public static final double SUM_GL_THRESH_NOCALL = -0.1;
    public final static String MERGE_INTERSECTION = "Intersection";
    public final static String MERGE_FILTER_IN_ALL = "FilteredInAll";
    public final static String MERGE_REF_IN_ALL = "ReferenceInAll";
    public final static String MERGE_FILTER_PREFIX = "filterIn";

    public enum GenotypeMergeType {
        /**
         * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
         */
        UNIQUIFY,
        /**
         * Take genotypes in priority order (see the priority argument).
         */
        PRIORITIZE,
        /**
         * Take the genotypes in any order.
         */
        UNSORTED,
        /**
         * Require that all samples/genotypes be unique between all inputs.
         */
        REQUIRE_UNIQUE
    }

    public enum FilteredRecordMergeType {
        /**
         * Union - leaves the record if any record is unfiltered.
         */
        KEEP_IF_ANY_UNFILTERED,
        /**
         * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
         */
        KEEP_IF_ALL_UNFILTERED,
        /**
         * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
         */
        KEEP_UNCONDITIONAL
    }

    static private void verifyUniqueSampleNames(Collection<VariantContext> unsortedVCs) {
        Set<String> names = new HashSet<String>();
        for (VariantContext vc : unsortedVCs) {
            for (String name : vc.getSampleNames()) {
                //System.out.printf("Checking %s %b%n", name, names.contains(name));
                if (names.contains(name))
                    throw new UserException("REQUIRE_UNIQUE sample names is true but duplicate names were discovered " + name);
            }

            names.addAll(vc.getSampleNames());
        }
    }

    public static List<VariantContext> sortVariantContextsByPriority(Collection<VariantContext> unsortedVCs, List<String> priorityListOfVCs, GenotypeMergeType mergeOption) {
        if (mergeOption == GenotypeMergeType.PRIORITIZE && priorityListOfVCs == null)
            throw new IllegalArgumentException("Cannot merge calls by priority with a null priority list");

        if (priorityListOfVCs == null || mergeOption == GenotypeMergeType.UNSORTED)
            return new ArrayList<VariantContext>(unsortedVCs);
        else {
            ArrayList<VariantContext> sorted = new ArrayList<VariantContext>(unsortedVCs);
            Collections.sort(sorted, new CompareByPriority(priorityListOfVCs));
            return sorted;
        }
    }

    static class CompareByPriority implements Comparator<VariantContext>, Serializable {
        /**
         *
         */
        private static final long serialVersionUID = 6195505964868678001L;
        List<String> priorityListOfVCs;

        public CompareByPriority(List<String> priorityListOfVCs) {
            this.priorityListOfVCs = priorityListOfVCs;
        }

        private int getIndex(VariantContext vc) {
            int i = priorityListOfVCs.indexOf(vc.getSource());
            if (i == -1)
                throw new UserException.BadArgumentValueException(Utils.join(",", priorityListOfVCs), "Priority list " + priorityListOfVCs + " doesn't contain variant context " + vc.getSource());
            return i;
        }

        public int compare(VariantContext vc1, VariantContext vc2) {
            return Integer.valueOf(getIndex(vc1)).compareTo(getIndex(vc2));
        }
    }
    
    /**
	 * Determines the common reference allele
	 *
	 * @param VCs
	 *            the list of VariantContexts
	 * @param loc
	 *            if not null, ignore records that do not begin at this start
	 *            location
	 * @return possibly null Allele
	 */
	public static Allele determineReferenceAllele(final List<VariantContext> VCs, final GenomeLocation loc) {
		Allele ref = null;

		for (final VariantContext vc : VCs) {
			if (contextMatchesLocation(vc, loc)) {
				final Allele myRef = vc.getReference();
				if (ref == null || ref.length() < myRef.length())
					ref = myRef;
				else if (ref.length() == myRef.length() && !ref.equals(myRef)) {
				    throw new TribbleException(String.format("The provided variant file(s) have inconsistent references for the same position(s) at %s:%d, %s vs. %s", vc.getContig(), vc.getStart(), ref, myRef));
                }
			}
		}

		return ref;
	}
	
	public static boolean contextMatchesLocation(final VariantContext vc, final GenomeLocation loc) {
		return loc == null || loc.getStart() == vc.getStart();
	}

    public static Allele determineReferenceAllele(List<VariantContext> VCs) {
        return determineReferenceAllele(VCs,null);
    }

    public static class AlleleMapper {
        private VariantContext vc = null;
        private Map<Allele, Allele> map = null;

        public AlleleMapper(VariantContext vc) {
            this.vc = vc;
        }

        public AlleleMapper(Map<Allele, Allele> map) {
            this.map = map;
        }

        public boolean needsRemapping() {
            return this.map != null;
        }

        public Collection<Allele> values() {
            return map != null ? map.values() : vc.getAlleles();
        }

        public Allele remap(Allele a) {
            return map != null && map.containsKey(a) ? map.get(a) : a;
        }

        public List<Allele> remap(List<Allele> as) {
            List<Allele> newAs = new ArrayList<Allele>();
            for (Allele a : as) {
                //System.out.printf("  Remapping %s => %s%n", a, remap(a));
                newAs.add(remap(a));
            }
            return newAs;
        }
        
        public List<Allele> getUniqueMappedAlleles() {
			if (map == null)
				return Collections.emptyList();
			return new ArrayList<>(new HashSet<>(map.values()));
		}
    }

    /**
     * create a genome location, given a variant context
     *
     * @param genomeLocParser parser
     * @param vc              the variant context
     * @return the genomeLoc
     */
    public static final GenomeLocation getLocation(GenomeLocationParser genomeLocParser, VariantContext vc) {
        return genomeLocParser.createGenomeLocation(vc.getChr(), vc.getStart(), vc.getEnd(), true);
    }


    static public AlleleMapper resolveIncompatibleAlleles(Allele refAllele, VariantContext vc, Set<Allele> allAlleles) {
        if (refAllele.equals(vc.getReference()))
            return new AlleleMapper(vc);
        else {
            // we really need to do some work.  The refAllele is the longest reference allele seen at this
            // start site.  So imagine it is:
            //
            // refAllele: ACGTGA
            // myRef:     ACGT
            // myAlt:     A
            //
            // We need to remap all of the alleles in vc to include the extra GA so that
            // myRef => refAllele and myAlt => AGA
            //

            Allele myRef = vc.getReference();
            if (refAllele.length() <= myRef.length())
                throw new UserException("BUG: myRef=" + myRef + " is longer than refAllele=" + refAllele);
            byte[] extraBases = Arrays.copyOfRange(refAllele.getBases(), myRef.length(), refAllele.length());

            // System.out.printf("Remapping allele at %s%n", vc);
            /// System.out.printf("ref   %s%n", refAllele);
            // System.out.printf("myref %s%n", myRef );
            // System.out.printf("extrabases %s%n", new String(extraBases));

            Map<Allele, Allele> map = new HashMap<Allele, Allele>();
            for (Allele a : vc.getAlleles()) {
                if (a.isReference())
                    map.put(a, refAllele);
                else {
                    Allele extended = Allele.extend(a, extraBases);
                    for (Allele b : allAlleles)
                        if (extended.equals(b))
                            extended = b;
                    //  System.out.printf("  Extending %s => %s%n", a, extended);
                    map.put(a, extended);
                }
            }

            // debugging
//            System.out.printf("mapping %s%n", map);

            return new AlleleMapper(map);
        }
    }

    protected static final boolean hasPLIncompatibleAlleles(final Collection<Allele> alleleSet1, final Collection<Allele> alleleSet2) {
        final Iterator<Allele> it1 = alleleSet1.iterator();
        final Iterator<Allele> it2 = alleleSet2.iterator();

        while (it1.hasNext() && it2.hasNext()) {
            final Allele a1 = it1.next();
            final Allele a2 = it2.next();
            if (!a1.equals(a2))
                return true;
        }

        // by this point, at least one of the iterators is empty.  All of the elements
        // we've compared are equal up until this point.  But it's possible that the
        // sets aren't the same size, which is indicated by the test below.  If they
        // are of the same size, though, the sets are compatible
        return it1.hasNext() || it2.hasNext();
    }

    public static GenotypesContext stripPLs(GenotypesContext genotypes) {
        GenotypesContext newGs = GenotypesContext.create(genotypes.size());

        for (final Genotype g : genotypes) {
            newGs.add(g.hasLikelihoods() ? removePLs(g) : g);
        }

        return newGs;
    }


    public static Genotype removePLs(Genotype g) {
        if (g.hasLikelihoods())
            return new GenotypeBuilder(g).noPL().make();
        else
            return g;
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the
     * proper chromosome-based VCF tags
     *
     * @param vc                the VariantContext
     * @param attributes        the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     * @return the attributes map provided as input, returned for programming convenience
     */
    public static Map<String, Object> calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues) {
        return calculateChromosomeCounts(vc, attributes, removeStaleValues, new HashSet<String>(0));
    }
    
    /**
     * Update the attributes of the attributes map in the VariantContextBuilder to reflect the proper
     * chromosome-based VCF tags based on the current VC produced by builder.make()
     *
     * @param builder     the VariantContextBuilder we are updating
     * @param removeStaleValues should we remove stale values from the mapping?
     */
    public static void calculateChromosomeCounts(VariantContextBuilder builder, boolean removeStaleValues) {
        VariantContext vc = builder.make();
        builder.attributes(calculateChromosomeCounts(vc, new HashMap<>(vc.getAttributes()), removeStaleValues, new HashSet<>(0)));
    }

    /**
     * Update the attributes of the attributes map given the VariantContext to reflect the
     * proper chromosome-based VCF tags
     *
     * @param vc                the VariantContext
     * @param attributes        the attributes map to populate; must not be null; may contain old values
     * @param removeStaleValues should we remove stale values from the mapping?
     * @param founderIds        - Set of founders Ids to take into account. AF and FC will be calculated over the founders.
     *                          If empty or null, counts are generated for all samples as unrelated individuals
     * @return the attributes map provided as input, returned for programming convenience
     */
    public static Map<String, Object> calculateChromosomeCounts(VariantContext vc, Map<String, Object> attributes, boolean removeStaleValues, final Set<String> founderIds) {
        final int AN = vc.getCalledChrCount();

        // if everyone is a no-call, remove the old attributes if requested
        if (AN == 0 && removeStaleValues) {
            if (attributes.containsKey(VCFConstants.ALLELE_COUNT_KEY))
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
            if (attributes.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY))
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            if (attributes.containsKey(VCFConstants.ALLELE_NUMBER_KEY))
                attributes.remove(VCFConstants.ALLELE_NUMBER_KEY);
            return attributes;
        }

        if (vc.hasGenotypes()) {
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, AN);

            // if there are alternate alleles, record the relevant tags
            if (vc.getAlternateAlleles().size() > 0) {
                ArrayList<Double> alleleFreqs = new ArrayList<Double>();
                ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
                ArrayList<Integer> foundersAlleleCounts = new ArrayList<Integer>();
                double totalFoundersChromosomes = (double) vc.getCalledChrCount(founderIds);
                int foundersAltChromosomes;
                for (Allele allele : vc.getAlternateAlleles()) {
                    foundersAltChromosomes = vc.getCalledChrCount(allele, founderIds);
                    alleleCounts.add(vc.getCalledChrCount(allele));
                    foundersAlleleCounts.add(foundersAltChromosomes);
                    if (AN == 0) {
                        alleleFreqs.add(0.0);
                    } else {
                        final Double freq = (double) foundersAltChromosomes / totalFoundersChromosomes;
                        alleleFreqs.add(freq);
                    }
                }

                attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs.size() == 1 ? alleleFreqs.get(0) : alleleFreqs);
            } else {
                // if there's no alt AC and AF shouldn't be present
                attributes.remove(VCFConstants.ALLELE_COUNT_KEY);
                attributes.remove(VCFConstants.ALLELE_FREQUENCY_KEY);
            }
        }

        return attributes;
    }

    protected static void mergeGenotypes(GenotypesContext mergedGenotypes, VariantContext oneVC, AlleleMapper alleleMapping, boolean uniqifySamples) {
        for (Genotype g : oneVC.getGenotypes()) {
            String name = mergedSampleName(oneVC.getSource(), g.getSampleName(), uniqifySamples);
            if (!mergedGenotypes.containsSample(name)) {
                // only add if the name is new
                Genotype newG = g;

                if (uniqifySamples || alleleMapping.needsRemapping()) {
                    final List<Allele> alleles = alleleMapping.needsRemapping() ? alleleMapping.remap(g.getAlleles()) : g.getAlleles();
                    newG = new GenotypeBuilder(g).name(name).alleles(alleles).make();
                }

                mergedGenotypes.add(newG);
            }
        }
    }

    public static String mergedSampleName(String trackName, String sampleName, boolean uniqify) {
        return uniqify ? sampleName + "." + trackName : sampleName;
    }

    /**
     * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
     * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
     * the sample name
     *
     * @param genomeLocParser         loc parser
     * @param unsortedVCs             collection of unsorted VCs
     * @param priorityListOfVCs       priority list detailing the order in which we should grab the VCs
     * @param filteredRecordMergeType merge type for filtered records
     * @param genotypeMergeOptions    merge option for genotypes
     * @param annotateOrigin          should we annotate the set it came from?
     * @param printMessages           should we print messages?
     * @param setKey                  the key name of the set
     * @param filteredAreUncalled     are filtered records uncalled?
     * @param mergeInfoWithMaxAC      should we merge in info from the VC with maximum allele count?
     * @return new VariantContext       representing the merge of unsortedVCs
     */
    public static VariantContext simpleMerge(final GenomeLocationParser genomeLocParser,
                                             final Collection<VariantContext> unsortedVCs,
                                             final List<String> priorityListOfVCs,
                                             final FilteredRecordMergeType filteredRecordMergeType,
                                             final GenotypeMergeType genotypeMergeOptions,
                                             final boolean annotateOrigin,
                                             final boolean printMessages,
                                             final String setKey,
                                             final boolean filteredAreUncalled,
                                             final boolean mergeInfoWithMaxAC) {
        if (unsortedVCs == null || unsortedVCs.size() == 0)
            return null;

        if (annotateOrigin && priorityListOfVCs == null)
            throw new IllegalArgumentException("Cannot merge calls and annotate their origins without a complete priority list of VariantContexts");

        if (genotypeMergeOptions == GenotypeMergeType.REQUIRE_UNIQUE)
            verifyUniqueSampleNames(unsortedVCs);

        final List<VariantContext> preFilteredVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
        // Make sure all variant contexts are padded with reference base in case of indels if necessary
        final List<VariantContext> VCs = new ArrayList<VariantContext>();

        for (final VariantContext vc : preFilteredVCs) {
            if (!filteredAreUncalled || vc.isNotFiltered())
                VCs.add(vc);
        }
        if (VCs.size() == 0) // everything is filtered out and we're filteredAreUncalled
            return null;

        // establish the baseline info from the first VC
        final VariantContext first = VCs.get(0);
        final String name = first.getSource();
        final Allele refAllele = determineReferenceAllele(VCs);
        //Byte referenceBaseForIndel = null;

        final Set<Allele> alleles = new LinkedHashSet<Allele>();
        final Set<String> filters = new HashSet<String>();
        final Map<String, Object> attributes = new LinkedHashMap<String, Object>();
        final Set<String> inconsistentAttributes = new HashSet<String>();
        final Set<String> variantSources = new HashSet<String>(); // contains the set of sources we found in our set of VCs that are variant
        final Set<String> rsIDs = new LinkedHashSet<String>(1); // most of the time there's one id

        GenomeLocation loc = getLocation(genomeLocParser, first);
        int depth = 0;
        int maxAC = -1;
        final Map<String, Object> attributesWithMaxAC = new LinkedHashMap<String, Object>();
        double log10PError = CommonInfo.NO_LOG10_PERROR;
        VariantContext vcWithMaxAC = null;
        GenotypesContext genotypes = GenotypesContext.create();

        // counting the number of filtered and variant VCs
        int nFiltered = 0;

        boolean remapped = false;

        // cycle through and add info from the other VCs, making sure the loc/reference matches

        for (final VariantContext vc : VCs) {
            if (loc.getStart() != vc.getStart())
                throw new UserException("BUG: attempting to merge VariantContexts with different start sites: first=" + first.toString() + " second=" + vc.toString());

            if (getLocation(genomeLocParser, vc).size() > loc.size())
                loc = getLocation(genomeLocParser, vc); // get the longest location

            nFiltered += vc.isFiltered() ? 1 : 0;
            if (vc.isVariant()) variantSources.add(vc.getSource());

            AlleleMapper alleleMapping = resolveIncompatibleAlleles(refAllele, vc, alleles);
            remapped = remapped || alleleMapping.needsRemapping();

            alleles.addAll(alleleMapping.values());

            mergeGenotypes(genotypes, vc, alleleMapping, genotypeMergeOptions == GenotypeMergeType.UNIQUIFY);

            // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
            if (log10PError == CommonInfo.NO_LOG10_PERROR)
                log10PError = vc.getLog10PError();

            filters.addAll(vc.getFilters());

            //
            // add attributes
            //
            // special case DP (add it up) and ID (just preserve it)
            //
            if (vc.hasAttribute(VCFConstants.DEPTH_KEY))
                depth += vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0);
            if (vc.hasID()) rsIDs.add(vc.getID());
            if (mergeInfoWithMaxAC && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                String rawAlleleCounts = vc.getAttributeAsString(VCFConstants.ALLELE_COUNT_KEY, null);
                // lets see if the string contains a , separator
                if (rawAlleleCounts.contains(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) {
                    List<String> alleleCountArray = Arrays.asList(rawAlleleCounts.substring(1, rawAlleleCounts.length() - 1).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR));
                    for (String alleleCount : alleleCountArray) {
                        final int ac = Integer.valueOf(alleleCount.trim());
                        if (ac > maxAC) {
                            maxAC = ac;
                            vcWithMaxAC = vc;
                        }
                    }
                } else {
                    final int ac = Integer.valueOf(rawAlleleCounts);
                    if (ac > maxAC) {
                        maxAC = ac;
                        vcWithMaxAC = vc;
                    }
                }
            }

            for (final Map.Entry<String, Object> p : vc.getAttributes().entrySet()) {
                String key = p.getKey();
                // if we don't like the key already, don't go anywhere
                if (!inconsistentAttributes.contains(key)) {
                    final boolean alreadyFound = attributes.containsKey(key);
                    final Object boundValue = attributes.get(key);
                    final boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                    if (alreadyFound && !boundValue.equals(p.getValue()) && !boundIsMissingValue) {
                        // we found the value but we're inconsistent, put it in the exclude list
                        //System.out.printf("Inconsistent INFO values: %s => %s and %s%n", key, boundValue, p.getValue());
                        inconsistentAttributes.add(key);
                        attributes.remove(key);
                    } else if (!alreadyFound || boundIsMissingValue) { // no value
                        //if ( vc != first ) System.out.printf("Adding key %s => %s%n", p.getKey(), p.getValue());
                        attributes.put(key, p.getValue());
                    }
                }
            }
        }

        // if we have more alternate alleles in the merged VC than in one or more of the
        // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF
        for (final VariantContext vc : VCs) {
            if (vc.getAlleles().size() == 1)
                continue;
            if (hasPLIncompatibleAlleles(alleles, vc.getAlleles())) {
                if (!genotypes.isEmpty())
                    System.out.println(String.format("Stripping PLs at %s due incompatible alleles merged=%s vs. single=%s",
                            genomeLocParser.createGenomeLocation(vc), alleles, vc.getAlleles()));
                genotypes = stripPLs(genotypes);
                // this will remove stale AC,AF attributed from vc
                calculateChromosomeCounts(vc, attributes, true);
                break;
            }
        }

        // take the VC with the maxAC and pull the attributes into a modifiable map
        if (mergeInfoWithMaxAC && vcWithMaxAC != null) {
            attributesWithMaxAC.putAll(vcWithMaxAC.getAttributes());
        }

        // if at least one record was unfiltered and we want a union, clear all of the filters
        if ((filteredRecordMergeType == FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size()) || filteredRecordMergeType == FilteredRecordMergeType.KEEP_UNCONDITIONAL)
            filters.clear();


        if (annotateOrigin) { // we care about where the call came from
            String setValue;
            if (nFiltered == 0 && variantSources.size() == priorityListOfVCs.size()) // nothing was unfiltered
                setValue = MERGE_INTERSECTION;
            else if (nFiltered == VCs.size())     // everything was filtered out
                setValue = MERGE_FILTER_IN_ALL;
            else if (variantSources.isEmpty())    // everyone was reference
                setValue = MERGE_REF_IN_ALL;
            else {
                final LinkedHashSet<String> s = new LinkedHashSet<String>();
                for (final VariantContext vc : VCs)
                    if (vc.isVariant())
                        s.add(vc.isFiltered() ? MERGE_FILTER_PREFIX + vc.getSource() : vc.getSource());
                setValue = Utils.join("-", s);
            }

            if (setKey != null) {
                attributes.put(setKey, setValue);
                if (mergeInfoWithMaxAC && vcWithMaxAC != null) {
                    attributesWithMaxAC.put(setKey, setValue);
                }
            }
        }

        if (depth > 0)
            attributes.put(VCFConstants.DEPTH_KEY, String.valueOf(depth));

        final String ID = rsIDs.isEmpty() ? VCFConstants.EMPTY_ID_FIELD : Utils.join(",", rsIDs);

        final VariantContextBuilder builder = new VariantContextBuilder().source(name).id(ID);
        builder.loc(loc.getContig(), loc.getStart(), loc.getStop());
        builder.alleles(alleles);
        builder.genotypes(genotypes);
        builder.log10PError(log10PError);
        builder.filters(filters.isEmpty() ? filters : new TreeSet<String>(filters));
        builder.attributes(new TreeMap<String, Object>(mergeInfoWithMaxAC ? attributesWithMaxAC : attributes));

        // Trim the padded bases of all alleles if necessary
        final VariantContext merged = builder.make();
        if (printMessages && remapped) System.out.printf("Remapped => %s%n", merged);
        return merged;
    }


    protected static final List<Allele> NO_CALL_ALLELES = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
    public static final int DEFAULT_PLOIDY = 2;

    /**
     * subset the Variant Context to the specific set of alleles passed in (pruning the PLs appropriately)
     *
     * @param vc              variant context with genotype likelihoods
     * @param allelesToUse    which alleles from the vc are okay to use; *** must be in the same relative order as those in the original VC ***
     * @param assignGenotypes true if we should update the genotypes based on the (subsetted) PLs
     * @return genotypes
     */
    @SuppressWarnings("deprecation")
    public static GenotypesContext subsetDiploidAlleles(final VariantContext vc,
                                                        final List<Allele> allelesToUse,
                                                        final boolean assignGenotypes) {

        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // which PLs should be carried forward?
        ArrayList<Integer> likelihoodIndexesToUse = null;

        // an optimization: if we are supposed to use all (or none in the case of a ref call) of the alleles,
        // then we can keep the PLs as is; otherwise, we determine which ones to keep
        if (numNewAltAlleles != numOriginalAltAlleles && numNewAltAlleles > 0) {
            likelihoodIndexesToUse = new ArrayList<Integer>(30);

            final boolean[] altAlleleIndexToUse = new boolean[numOriginalAltAlleles];
            for (int i = 0; i < numOriginalAltAlleles; i++) {
                if (allelesToUse.contains(vc.getAlternateAllele(i)))
                    altAlleleIndexToUse[i] = true;
            }

            // numLikelihoods takes total # of alleles. Use default # of chromosomes (ploidy) = 2
            final int numLikelihoods = GenotypeLikelihoods.numLikelihoods(1 + numOriginalAltAlleles, DEFAULT_PLOIDY);
            for (int PLindex = 0; PLindex < numLikelihoods; PLindex++) {
                final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);
                // consider this entry only if both of the alleles are good
                if ((alleles.alleleIndex1 == 0 || altAlleleIndexToUse[alleles.alleleIndex1 - 1]) && (alleles.alleleIndex2 == 0 || altAlleleIndexToUse[alleles.alleleIndex2 - 1]))
                    likelihoodIndexesToUse.add(PLindex);
            }
        }

        // create the new genotypes
        for (int k = 0; k < oldGTs.size(); k++) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            if (!g.hasLikelihoods()) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;
            if (likelihoodIndexesToUse == null) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = new double[likelihoodIndexesToUse.size()];
                int newIndex = 0;
                for (int oldIndex : likelihoodIndexesToUse)
                    newLikelihoods[newIndex++] = originalLikelihoods[oldIndex];

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if (MathUtils.sum(newLikelihoods) > SUM_GL_THRESH_NOCALL) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), NO_CALL_ALLELES));
            } else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if (numNewAltAlleles == 0)
                    gb.noPL();
                else
                    gb.PL(newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if (!assignGenotypes || MathUtils.sum(newLikelihoods) > SUM_GL_THRESH_NOCALL) {
                    gb.alleles(NO_CALL_ALLELES);
                } else {
                    // find the genotype with maximum likelihoods
                    int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
                    GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);

                    gb.alleles(Arrays.asList(allelesToUse.get(alleles.alleleIndex1), allelesToUse.get(alleles.alleleIndex2)));
                    if (numNewAltAlleles != 0)
                        gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
                }
                newGTs.add(gb.make());
            }
        }

        return newGTs;
    }

    public static VariantContext reverseTrimAlleles(final VariantContext inputVC) {

        // TODO - this function doesn't work with mixed records or records that started as mixed and then became non-mixed

        // see whether we need to trim common reference base from all alleles

        final int trimExtent = computeReverseClipping(inputVC.getAlleles(), inputVC.getReference().getDisplayString().getBytes(), 0, false);
        if (trimExtent <= 0 || inputVC.getAlleles().size() <= 1)
            return inputVC;

        final List<Allele> alleles = new ArrayList<Allele>();
        final GenotypesContext genotypes = GenotypesContext.create();
        final Map<Allele, Allele> originalToTrimmedAlleleMap = new HashMap<Allele, Allele>();

        for (final Allele a : inputVC.getAlleles()) {
            if (a.isSymbolic()) {
                alleles.add(a);
                originalToTrimmedAlleleMap.put(a, a);
            } else {
                // get bases for current allele and create a new one with trimmed bases
                final byte[] newBases = Arrays.copyOfRange(a.getBases(), 0, a.length() - trimExtent);
                final Allele trimmedAllele = Allele.create(newBases, a.isReference());
                alleles.add(trimmedAllele);
                originalToTrimmedAlleleMap.put(a, trimmedAllele);
            }
        }

        // now we can recreate new genotypes with trimmed alleles
        for (final Genotype genotype : inputVC.getGenotypes()) {
            final List<Allele> originalAlleles = genotype.getAlleles();
            final List<Allele> trimmedAlleles = new ArrayList<Allele>();
            for (final Allele a : originalAlleles) {
                if (a.isCalled())
                    trimmedAlleles.add(originalToTrimmedAlleleMap.get(a));
                else
                    trimmedAlleles.add(Allele.NO_CALL);
            }
            genotypes.add(new GenotypeBuilder(genotype).alleles(trimmedAlleles).make());
        }

        return new VariantContextBuilder(inputVC).stop(inputVC.getStart() + alleles.get(0).length() - 1).alleles(alleles).genotypes(genotypes).make();
    }

    public static int computeReverseClipping(final List<Allele> unclippedAlleles,
                                             final byte[] ref,
                                             final int forwardClipping,
                                             final boolean allowFullClip) {
        int clipping = 0;
        boolean stillClipping = true;

        while (stillClipping) {
            for (final Allele a : unclippedAlleles) {
                if (a.isSymbolic())
                    continue;

                // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
                // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
                if (a.length() - clipping == 0)
                    return clipping - (allowFullClip ? 0 : 1);

                if (a.length() - clipping <= forwardClipping || a.length() - forwardClipping == 0) {
                    stillClipping = false;
                } else if (ref.length == clipping) {
                    if (allowFullClip)
                        stillClipping = false;
                    else
                        return -1;
                } else if (a.getBases()[a.length() - clipping - 1] != ref[ref.length - clipping - 1]) {
                    stillClipping = false;
                }
            }
            if (stillClipping)
                clipping++;
        }
        return clipping;
    }

    /**
     *
     * @param vc
     * @param refBasesStartingAtVCWithPad
     * @return
     */
    //@Requires({"vc != null", "refBasesStartingAtVCWithPad != null && refBasesStartingAtVCWithPad.length > 0"})
    public static Pair<List<Integer>,byte[]> getNumTandemRepeatUnits(final VariantContext vc, final byte[] refBasesStartingAtVCWithPad) {
        final boolean VERBOSE = false;
        final String refBasesStartingAtVCWithoutPad = new String(refBasesStartingAtVCWithPad).substring(1);
        if ( ! vc.isIndel() ) // only indels are tandem repeats
            return null;

        final Allele refAllele = vc.getReference();
        final byte[] refAlleleBases = Arrays.copyOfRange(refAllele.getBases(), 1, refAllele.length());

        byte[] repeatUnit = null;
        final ArrayList<Integer> lengths = new ArrayList<Integer>();

        for ( final Allele allele : vc.getAlternateAlleles() ) {
            Pair<int[],byte[]> result = getNumTandemRepeatUnits(refAlleleBases, Arrays.copyOfRange(allele.getBases(), 1, allele.length()), refBasesStartingAtVCWithoutPad.getBytes());

            final int[] repetitionCount = result.first;
            // repetition count = 0 means allele is not a tandem expansion of context
            if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
                return null;

            if (lengths.size() == 0) {
                lengths.add(repetitionCount[0]); // add ref allele length only once
            }
            lengths.add(repetitionCount[1]);  // add this alt allele's length

            repeatUnit = result.second;
            if (VERBOSE) {
                System.out.println("RefContext:"+refBasesStartingAtVCWithoutPad);
                System.out.println("Ref:"+refAllele.toString()+" Count:" + String.valueOf(repetitionCount[0]));
                System.out.println("Allele:"+allele.toString()+" Count:" + String.valueOf(repetitionCount[1]));
                System.out.println("RU:"+new String(repeatUnit));
            }
        }

        return new Pair<List<Integer>, byte[]>(lengths,repeatUnit);
    }

    protected static Pair<int[],byte[]> getNumTandemRepeatUnits(final byte[] refBases, final byte[] altBases, final byte[] remainingRefContext) {
         /* we can't exactly apply same logic as in basesAreRepeated() to compute tandem unit and number of repeated units.
           Consider case where ref =ATATAT and we have an insertion of ATAT. Natural description is (AT)3 -> (AT)5.
         */

        byte[] longB;
        // find first repeat unit based on either ref or alt, whichever is longer
        if (altBases.length > refBases.length)
            longB = altBases;
        else
            longB = refBases;

        // see if non-null allele (either ref or alt, whichever is longer) can be decomposed into several identical tandem units
        // for example, -*,CACA needs to first be decomposed into (CA)2
        final int repeatUnitLength = findRepeatedSubstring(longB);
        final byte[] repeatUnit = Arrays.copyOf(longB, repeatUnitLength);

        final int[] repetitionCount = new int[2];
//        repetitionCount[0] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(refBases, remainingRefContext));
//        repetitionCount[1] = findNumberofRepetitions(repeatUnit, ArrayUtils.addAll(altBases, remainingRefContext));
        int repetitionsInRef = findNumberofRepetitions(repeatUnit,refBases);
        repetitionCount[0] = findNumberofRepetitions(repeatUnit, org.apache.commons.lang.ArrayUtils.addAll(refBases, remainingRefContext))-repetitionsInRef;
        repetitionCount[1] = findNumberofRepetitions(repeatUnit, org.apache.commons.lang.ArrayUtils.addAll(altBases, remainingRefContext))-repetitionsInRef;

        return new Pair<int[], byte[]>(repetitionCount, repeatUnit);

    }

    /**
     * Find out if a string can be represented as a tandem number of substrings.
     * For example ACTACT is a 2-tandem of ACT,
     * but ACTACA is not.
     *
     * @param bases                 String to be tested
     * @return                      Length of repeat unit, if string can be represented as tandem of substring (if it can't
     *                              be represented as one, it will be just the length of the input string)
     */
    protected static int findRepeatedSubstring(byte[] bases) {

        int repLength;
        for (repLength=1; repLength <=bases.length; repLength++) {
            final byte[] candidateRepeatUnit = Arrays.copyOf(bases,repLength);
            boolean allBasesMatch = true;
            for (int start = repLength; start < bases.length; start += repLength ) {
                // check that remaining of string is exactly equal to repeat unit
                final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
                if (!Arrays.equals(candidateRepeatUnit, basePiece)) {
                    allBasesMatch = false;
                    break;
                }
            }
            if (allBasesMatch)
                return repLength;
        }

        return repLength;
    }

    /**
     * Helper routine that finds number of repetitions a string consists of.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnit             Substring
     * @param testString             String to test
     * @return                       Number of repetitions (0 if testString is not a concatenation of n repeatUnit's
     */
    protected static int findNumberofRepetitions(byte[] repeatUnit, byte[] testString) {
        int numRepeats = 0;
        for (int start = 0; start < testString.length; start += repeatUnit.length) {
            int end = start + repeatUnit.length;
            byte[] unit = Arrays.copyOfRange(testString,start, end);
            if(Arrays.equals(unit,repeatUnit))
                numRepeats++;
            else
                return numRepeats;
        }
        return numRepeats;
    }
}