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

import org.broadinstitute.hellbender.utils.SimpleInterval;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class VariantSiteSet {
    private static final Logger logger = LoggerFactory.getLogger(VariantSiteSet.class);
    public SimpleInterval interval;
    public BitSet data;

    public VariantSiteSet(final SimpleInterval interval) {
        this.interval = interval;
        this.data = new BitSet(this.interval.size());
    }

    

}
