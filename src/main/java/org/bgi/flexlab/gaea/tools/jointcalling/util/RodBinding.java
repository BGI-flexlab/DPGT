package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.HashMap;
import java.util.Map;

import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.tribble.Feature;

public final class RodBinding<T extends Feature> {
    protected final static String UNBOUND_VARIABLE_NAME = "";
    protected final static String UNBOUND_SOURCE = "UNBOUND";
    protected final static String UNBOUND_TRIBBLE_TYPE = "";

    /**
     * Create an unbound Rodbinding of type.  This is the correct programming
     * style for an optional RodBinding<T>
     *
     *     At Input()
     *     RodBinding<T> x = RodBinding.makeUnbound(T.class)
     *
     * The unbound binding is guaranteed to never match any binding.  It uniquely
     * returns false to isBound().
     *
     * @param type the Class type produced by this unbound object
     * @param <T> any class extending Tribble Feature
     * @return the UNBOUND RodBinding producing objects of type T
     */
    protected final static <T extends Feature> RodBinding<T> makeUnbound(Class<T> type) {
        return new RodBinding<T>(type);
    }

    /** The name of this binding.  Often the name of the field itself, but can be overridden on cmdline */
    final private String name;
    /** where the data for this ROD is coming from.  A file or special value if coming from stdin */
    final private String source;
    /** the string name of the tribble type, such as vcf, bed, etc. */
    final private String tribbleType;
    /** The command line tags associated with this RodBinding */
    final private Tags tags;
    /** The Java class expected for this RodBinding.  Must correspond to the type emitted by Tribble */
    final private Class<T> type;
    /** True for all RodBindings except the special UNBOUND binding, which is the default for optional arguments */
    final private boolean bound;

    /**
     * The name counter.  This is how we create unique names for collections of RodBindings
     * on the command line.  If you have provide the GATK with -X file1 and -X file2 to a
     * RodBinding argument as List<RodBinding<T>> then each binding will receive automatically
     * the name of X and X2.
     */
    final private static Map<String, Integer> nameCounter = new HashMap<String, Integer>();

    /** for UnitTests */
    final public static void resetNameCounter() {
        nameCounter.clear();
    }

    final private static synchronized String countedVariableName(final String rawName) {
    	Utils.nonNull(rawName, "raw name cann't be null!");
        Integer count = nameCounter.get(rawName);
        if ( count == null ) {
            nameCounter.put(rawName, 1);
            return rawName;
        } else {
            nameCounter.put(rawName, count + 1);
            return rawName + (count + 1);
        }
    }

    public RodBinding(Class<T> type, final String rawName, final String source, final String tribbleType, final Tags tags) {
    	Utils.nonNull(type, "type cann't be null!");
    	Utils.nonNull(rawName, "raw name cann't be null!");
    	Utils.nonNull(source, "source cann't be null!");
    	Utils.nonNull(tribbleType, "tribbleType cann't be null!");
    	Utils.nonNull(tags, "tags cann't be null!");
        this.type = type;
        this.name = countedVariableName(rawName);
        this.source = source;
        this.tribbleType = tribbleType;
        this.tags = tags;
        this.bound = true;
    }

    /**
     * For testing purposes only.  Creates a RodBinding sufficient for looking up associations to rawName
     * @param type
     * @param rawName
     */
    public RodBinding(Class<T> type, final String rawName) {
        this(type, rawName, "missing", type.getSimpleName(), new Tags());
    }

    /**
     * Make an unbound RodBinding<T>.  Only available for creating the globally unique UNBOUND object
     * @param type class this unbound RodBinding creates
     */
    private RodBinding(Class<T> type) {
    	Utils.nonNull(type, "type cann't be null!");
        this.type = type;
        this.name = UNBOUND_VARIABLE_NAME;  // special value can never be found in RefMetaDataTracker
        this.source = UNBOUND_SOURCE;
        this.tribbleType = UNBOUND_TRIBBLE_TYPE;
        this.tags = new Tags();
        this.bound = false;
    }


   /**
     * @return True for all RodBindings except the special UNBOUND binding, which is the default for optional arguments
     */
    final public boolean isBound() {
        return bound;
    }

    /**
     * @return The name of this binding.  Often the name of the field itself, but can be overridden on cmdline
     */
    final public String getName() {
        return name;
    }

    /**
     * @return the string name of the tribble type, such as vcf, bed, etc.
     */
    final public Class<T> getType() {
        return type;
    }

    /**
     * @return where the data for this ROD is coming from.  A file or special value if coming from stdin
     */
    final public String getSource() {
        return source;
    }

    /**
     * @return The command line tags associated with this RodBinding.  Will include the tags used to
     * determine the name and type of this RodBinding
     */
    final public Tags getTags() {
        return tags;
    }

    /**
     * @return The Java class expected for this RodBinding.  Must correspond to the type emited by Tribble
     */
    final public String getTribbleType() {
        return tribbleType;
    }

    @Override
    public String toString() {
        return String.format("(RodBinding name=%s source=%s)", getName(), getSource());
    }
}

