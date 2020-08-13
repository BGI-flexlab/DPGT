package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.lang.reflect.Constructor;
import java.util.Collection;
import java.util.Collections;

import htsjdk.tribble.Feature;

public final class RodBindingCollection<T extends Feature> {

    /** The Java class expected for this RodBinding.  Must correspond to the type emitted by Tribble */
    final private Class<T> type;

    private Collection<RodBinding<T>> rodBindings;

    public RodBindingCollection(final Class<T> type, final Collection<RodBinding<T>> rodBindings) {
        this.type = type;
        this.rodBindings = Collections.unmodifiableCollection(rodBindings);
    }

    /**
     * @return the collection of RodBindings
     */
    final public Collection<RodBinding<T>> getRodBindings() {
        return rodBindings;
    }

    /**
     * @return the string name of the tribble type, such as vcf, bed, etc.
     */
    final public Class<T> getType() {
        return type;
    }

    @Override
    public String toString() {
        return String.format("(RodBindingCollection %s)", getRodBindings());
    }

    /**
     * Utility method to help construct a RodBindingCollection of the given Feature type
     *
     * @param type         the Feature type
     * @param rodBindings  the rod bindings to put into the collection
     * @return a new RodBindingCollection object
     */
    public static Object createRodBindingCollectionOfType(final Class<? extends Feature> type, final Collection<RodBinding> rodBindings) {
        try {
            final Constructor ctor = RodBindingCollection.class.getConstructor(Class.class, Collection.class);
            return ctor.newInstance(type, rodBindings);
        } catch (final Exception e) {
            throw new IllegalStateException("Failed to create a RodBindingCollection for type " + type);
        }
    }
}

