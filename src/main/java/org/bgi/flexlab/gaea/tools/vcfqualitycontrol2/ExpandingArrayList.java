package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2;

import java.util.ArrayList;
import java.util.Collection;

public class ExpandingArrayList<E> extends ArrayList<E> {
    private static final long serialVersionUID = 1L;

    public ExpandingArrayList() { super(); }
    public ExpandingArrayList(Collection<? extends E> c) { super(c); }
    public ExpandingArrayList(int initialCapacity) { super(initialCapacity); }

    /**
     * Returns the element at the specified position in this list.  If index > size,
     * returns null.  Otherwise tries to access the array
     * @param index
     * @return
     * @throws IndexOutOfBoundsException in index < 0
     */
    @Override
    public E get(int index) throws IndexOutOfBoundsException {
        if ( index < size() )
            return super.get(index);
        else
            return null;
    }

    public E expandingGet(int index, E default_value) throws IndexOutOfBoundsException {
        maybeExpand(index, default_value);
        return super.get(index);
    }

    @Override
    public E set(int index, E element) {
        maybeExpand(index, null);
        return super.set(index, element);
    }
    private void maybeExpand(int index, E value) {
        if ( index >= size() ) {
            ensureCapacity(index+1); // make sure we have space to hold at least index + 1 elements
            // We need to add null items until we can safely set index to element
            for ( int i = size(); i <= index; i++ )
                add(value);
        }
    }

}
