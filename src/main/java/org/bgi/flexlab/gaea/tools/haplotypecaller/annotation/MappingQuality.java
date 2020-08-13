package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.List;
import java.util.OptionalInt;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.GaeaProtectedMathUtils;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.primitives.Ints;

import htsjdk.variant.variantcontext.VariantContext;

public class MappingQuality extends PerAlleleAnnotation implements StandardMutectAnnotation {
    public static final String KEY = "MMQ";

    @Override
    protected int aggregate(final List<Integer> values) {
        return values.isEmpty() ? 0 : GaeaProtectedMathUtils.median(Ints.toArray(values));
    }

    @Override
    protected String getVcfKey() { return KEY; }

    @Override
    protected String getDescription() { return "median mapping quality"; }

    @Override
    protected OptionalInt getValueForRead(final GaeaSamRecord read, final VariantContext vc) {
        Utils.nonNull(read);
        return OptionalInt.of(read.getMappingQuality());
    }
}

