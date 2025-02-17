/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;


public class VariantSiteFinderSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String prefix;
    public String chrom;
    public int start;
    public int end;

    static {
        NativeLibraryLoader.load();
    }

    /**
     * variant site finder spark function
     * @param chrom chromosome
     * @param start 0-based start
     * @param end 0-based end
     */
    public VariantSiteFinderSparkFunc(final String prefix, final String chrom, int start, int end) {
        this.prefix = prefix;
        this.chrom = chrom;
        this.start = start;
        this.end = end;
    }
    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        if (vcfpaths.isEmpty()) {
            // no vcf in this part, return a null output path
            ArrayList<String> result = new ArrayList<>();
            result.add("null");
            return result.iterator();
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);
        String outpath = prefix + idx + JointCallingSparkConsts.VARIANT_SITE_SUFFIX;
        VariantSiteFinder vf = new VariantSiteFinder();
        // call native c++ function to find variant site, return variant site bitset as byte array
        vf.FindVariantSite(vcfpathsArray, outpath, this.chrom, (long)this.start, (long)this.end);
        ArrayList<String> result = new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
