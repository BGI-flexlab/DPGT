package org.bgi.flexlab.gaea.tools.jointcalling.genotypelikelihood;

import java.util.Arrays;

public class IntMaxHeap {
	private int size;

    private int[] values;

    /**
     * Creates a new empty heap indicating its initial capacity.
     * @param initialCapacity number of elements you expect to have at most in the heap.
     *
     * @throws IllegalArgumentException if {@code initialCapacity} is negative.
     */
    public IntMaxHeap(final int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException();
        // We force it to have at least length 1 so that the capacity expansion works when adding;
        // it doubles current length and twice 0 = 0.
        values = new int[initialCapacity == 0 ? 1 : initialCapacity];
    }

    /**
     * Adds a new element to the heap.
     *
     * <p>The heap with grow if it runs out of capacity to hold the new element</p>
     *
     * @param v the new element.
     */
    public void add(final int v) {
        // Double capacity if overflow:
        ensureCapacity(size + 1);
        addWithoutCheckingCapacity(v);
    }

    /**
     * Implements the heap addition floating up the value.
     * @param v the value to add.
     */
    private void addWithoutCheckingCapacity(final int v) {
        int p;
        values[p = size++] = v;

        // Float up the recently added element:
        while (p > 0) {
            final int q = (p - 1) >> 1;      // parent index.
            final int u = values[q];         // parent value.

            //Finish check and update:
            if (u >= v)
                break;
            values[p] = u;
            values[q] = v;
            p = q;
        }
    }

    /**
     * Add several integers into the heap.
     * @param v values to add.
     */
    public void add(final int ... v) {
        if (v == null)
            throw new IllegalArgumentException("the input array cannot be null");
        ensureCapacity(v.length + size);
        for (int i : v)
            addWithoutCheckingCapacity(i);
    }

    private void ensureCapacity(final int newSize) {
        if (newSize > values.length)
            values = Arrays.copyOf(values,Math.max(newSize,10 + values.length << 1));
    }

    /**
     * Returns the current minimum element.
     *
     * @throws IllegalStateException if the heap is empty.
     *
     * @return the minimum element in the heap.
     */
    public int peek() {
        if (size == 0)
            throw new IllegalStateException("the heap is empty");
        return values[0];
    }

    /**
     * Returns the minimum element of the heap and removes it.
     *
     * @throws IllegalStateException if the heap is empty.
     *
     * @return the minimum element in the heap before removing it.
     */
    public int remove() {
        if (size == 0)
            throw new IllegalArgumentException("the heap is empty");
        final int result = values[0];
        removeUpdate();
        return result;
    }

    /**
     * Updates the heap after a removal, sinking the last element from the top-down.
     */
    private void removeUpdate() {
        // if the remove make the heap to be empty there is nothing to do.
        if (--size == 0)
            return;

        final int v = values[size]; // the last value.

        int p;
        values[p = 0] = v;

        // limit := first index in the heap that does not have any descendants within the heap.
        final int limit = (size >> 1);

        // Sorry! for the big loop but doesn't seem to be any other *practical* option that would reduce its size.
        while (p < limit) {

            // Initialize variables:
            final int r = (p + 1) << 1; // left descendant index.
            final int l = r - 1;        // right descendant index (no guarantee to be in the heap).
            int u = v;                  // will contain min(v,values[l],values[r]).
            int q = p;                  // wilL contain argmin_x(values[x], x in {p,l,r}).

            // Check left descendant:
            int lv = values[l];   // left descendant value.
            if (lv > u) {         // is the left descendant'v value more than v.
                u = lv;
                q = l;
            }

            // Check right descendant:
            if (r < size) {        // make sure that r is within the heap.
                int rv = values[r];
                if (rv > u) {      // is the right descendant's value less than v or left's
                    u = rv;
                    q = r;
                }
            }

            // Finish check and update:
            if (p == q)           // q == p if neither left or right descendants are less than v.
                break;

            values[p] = u;
            values[q] = v;
            p = q;
        }
    }

    /**
     * Checks whether the heap is empty.
     *
     * @return {@code true} iff the heap is empty.
     */
    public boolean isEmpty() {
        return size == 0;
    }

    /**
     * Returns the current size of the heap.
     *
     * @return 0 or greater.
     */
    public int size() {
        return size;
    }

    /**
     * Removes all elements from the heap.
     */
    public void clear() {
        size = 0;
    }
}
