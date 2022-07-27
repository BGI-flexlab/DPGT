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
