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
import java.lang.String;

public class VariantSiteFinder {
    /**
     * JNI for finding variant site of input vcf files and input genome region
     * @param vcfpaths input vcf files
     * @param outpath output variant site bit set bytes in binary format
     * @param chrom chromosome
     * @param start 0-based start of the region
     * @param end   0-based end of the region
     * @return
     */
    public native void FindVariantSite(
        String[] vcfpaths, String outpath, String chrom, long start, long end);
}
