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
package org.bgi.flexlab.gaea.data.structure.sequenceplatform;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

public enum NGSPlatform {
    ILLUMINA("ILLUMINA", "SLX", "SOLEXA"),
    SOLID("SOLID"),
    LS454("454"),
    COMPLETE_GENOMICS("COMPLETE"),
    PACBIO("PACBIO"),
    ION_TORRENT("IONTORRENT"),
    UNKNOWN("UNKNOWN");

    /**
     * Array of the prefix names in a BAM file for each of the platforms.
     */
    private final String[] BAM_PL_NAMES;

    private NGSPlatform(final String... BAM_PL_NAMES) {
        for ( int i = 0; i < BAM_PL_NAMES.length; i++ )
            BAM_PL_NAMES[i] = BAM_PL_NAMES[i].toUpperCase();
        this.BAM_PL_NAMES = BAM_PL_NAMES;
    }

    /**
     * Returns a representative PL string for this platform
     * @return
     */
    public final String getDefaultPlatform() {
        return BAM_PL_NAMES[0];
    }

    /**
     * Convenience constructor -- calculates the NGSPlatfrom from a SAMRecord.
     * Note you should not use this function if you have a GATKSAMRecord -- use the
     * accessor method instead.
     */
    public static final NGSPlatform fromRead(SAMRecord read) {
        return fromReadGroup(read.getReadGroup());
    }

    /**
     * Returns the NGSPlatform corresponding to the PL tag in the read group
     * @param rg
     * @return an NGSPlatform object matching the PL field of the header, of UNKNOWN if there was no match
     */
    public static final NGSPlatform fromReadGroup(SAMReadGroupRecord rg) {
        if ( rg == null ) return UNKNOWN;
        return fromReadGroupPL(rg.getPlatform());
    }

    /**
     * Returns the NGSPlatform corresponding to the PL tag in the read group
     */
    public static final NGSPlatform fromReadGroupPL(final String plFromRG) {
        if ( plFromRG == null ) return UNKNOWN;

        final String pl = plFromRG.toUpperCase();
        for ( final NGSPlatform ngsPlatform : NGSPlatform.values() ) {
            for ( final String bamPLName : ngsPlatform.BAM_PL_NAMES ) {
                if ( pl.contains(bamPLName) )
                    return ngsPlatform;
            }
        }

        return UNKNOWN;
    }

    /**
     * checks whether or not the requested platform is listed in the set (and is not unknown)
     */
    public static final boolean isKnown (final String platform) {
        return fromReadGroupPL(platform) != UNKNOWN;
    }
}
