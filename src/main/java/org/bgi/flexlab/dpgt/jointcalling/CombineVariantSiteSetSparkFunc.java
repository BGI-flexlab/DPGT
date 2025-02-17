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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.BitSet;
import org.apache.spark.api.java.function.Function2;
import org.bgi.flexlab.dpgt.utils.VariantSiteSetUtils;


public class CombineVariantSiteSetSparkFunc implements Function2<Integer, Iterator<String>, Iterator<BitSet>> {
    public String prefix;
    public CombineVariantSiteSetSparkFunc(final String prefix) {
        this.prefix = prefix;
    }

    @Override public Iterator<BitSet> call(Integer idx, Iterator<String> vcfpathIter) {
        BitSet bitSet = VariantSiteSetUtils.combine(vcfpathIter);
        ArrayList<BitSet> result = new ArrayList<>();
        result.add(bitSet);
        return result.iterator();
    }
}
