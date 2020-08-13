package org.bgi.flexlab.gaea.util;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntConsumer;
import java.util.function.IntPredicate;
import java.util.function.IntToDoubleFunction;
import java.util.function.IntUnaryOperator;

import org.bgi.flexlab.gaea.tools.jointcalling.util.JointCallingUtils;

public final class IndexRange {

    /**
     * First index in the range.
     * <p>
     *     It won't ever be negative nor greater than {@link #to}.
     * </p>
     */
    public final int from;

    /**
     * Index following the last index included in the range.
     *
     * <p>
     *     It won't ever be negative nor less than {@link #from}.
     * </p>
     */
    public final int to;

    /**
     * Creates a new range given its {@code from} and {@code to} indices.
     *
     * @param fromIndex the {@code from} index value.
     * @param toIndex   the {@code to} index value.
     * @throws IllegalArgumentException if {@code fromIndex} is larger than {@code toIndex} or either is
     *                                  negative.
     */
    public IndexRange(final int fromIndex, final int toIndex) {
    	JointCallingUtils.validateArg(fromIndex <= toIndex, "the range size cannot be negative");
    	JointCallingUtils.validateArg(fromIndex >= 0, "the range cannot contain negative indices");
        from = fromIndex;
        to = toIndex;
    }
    
    public boolean isValidLength(final int length) {
        return to <= length;
    }

    /**
     * Returns number indexes expanded by this range.
     *
     * @return 0 or greater.
     */
    public int size() {
        return to - from;
    }

    /**
     * Iterate through all indexes in the range in ascending order to be processed by the
     * provided {@link IntConsumer integer consumer} lambda function.
     *
     * <p>
     *     Exceptions thrown by the execution of the index consumer {@code lambda}
     *     will be propagated to the caller immediately thus stopping early and preventing
     *     further indexes to be processed.
     * </p>
     * @param lambda the index consumer lambda.
     * @throws IllegalArgumentException if {@code lambda} is {@code null}.
     * @throws RuntimeException if thrown by {@code lambda} for some index.
     * @throws Error if thrown by {@code lambda} for some index.
     */
    public void forEach(final IntConsumer lambda) {
    	JointCallingUtils.nonNull(lambda, "the lambda function cannot be null");
        for (int i = from; i < to; i++) {
            lambda.accept(i);
        }
    }
    
    /**
     * Apply an int -> double function to this range, producing a double[]
     *
     * @param lambda the int -> double function
     */
    public double[] mapToDouble(final IntToDoubleFunction lambda) {
    	JointCallingUtils.nonNull(lambda, "the lambda function cannot be null");
        final double[] result = new double[size()];
        for (int i = from; i < to; i++) {
            result[i - from] = lambda.applyAsDouble(i);
        }
        return result;
    }

    /**
     * Sums the values of an int -> double function applied to this range
     *
     * @param lambda the int -> double function
     */
    public double sum(final IntToDoubleFunction lambda) {
    	JointCallingUtils.nonNull(lambda, "the lambda function cannot be null");
        double result = 0;
        for (int i = from; i < to; i++) {
            result += lambda.applyAsDouble(i);
        }
        return result;
    }

    /**
     * Apply an int -> int function to this range, producing an int[]
     *
     * @param lambda the int -> int function
     */
    public int[] mapToInteger(final IntUnaryOperator lambda) {
    	JointCallingUtils.nonNull(lambda, "the lambda function cannot be null");
        final int[] result = new int[size()];
        for (int i = from; i < to; i++) {
            result[i - from] = lambda.applyAsInt(i);
        }
        return result;
    }
    
    /**
     * Find the elements of this range for which an int -> boolean predicate is true
     *
     * @param predicate the int -> boolean predicate
     * @return
     */
    public List<Integer> filter(final IntPredicate predicate) {
    	JointCallingUtils.nonNull(predicate, "predicate may not be null");
        final List<Integer> result = new ArrayList<>();
        forEach(i -> {if (predicate.test(i)) result.add(i); } );
        return result;
    }

    @Override
    public boolean equals(final Object other) {
        if (other == this) {
            return true;
        } else if (!(other instanceof IndexRange)) {
            return false;
        } else {
            final IndexRange otherCasted = (IndexRange) other;
            return otherCasted.from == this.from && otherCasted.to == this.to;
        }
    }

    @Override
    public int hashCode() {
        // Inspired on {@link Arrays#hashCode(Object[])}.
        return (( 31 + Integer.hashCode(from) ) * 31 ) + Integer.hashCode(to);
    }
    
    @Override
    public String toString() {
        return String.format("%d-%d",from,to);
    }
}
