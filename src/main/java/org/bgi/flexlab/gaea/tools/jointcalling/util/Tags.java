package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Tags {
	 /**
     * Storage for the ordered, unkeyed, positional tags.
     */
    private final List<String> positionalTags = new ArrayList<String>();

    /**
     * Storage for key-value tags of the form <key>=<value>
     */
    private Map<String,String> keyValueTags = new HashMap<String,String>();

    /**
     * Tests to see whether two tag sets are equal.
     * @param other Other object to test for equality.
     * @return True if objects are the same.  False if objects differ.
     */
    @Override
    public boolean equals(Object other) {
        if(other == null)
            return false;

        if(!(other instanceof Tags))
            return false;

        Tags otherTags = (Tags)other;
        return this.positionalTags.equals(otherTags.positionalTags) && this.keyValueTags.equals(otherTags.keyValueTags);
    }

    /**
     * Returns whether any tags are specified on the command-line for this operation.
     * @return True if the tags are empty; false otherwise.
     */
    public boolean isEmpty() {
        return positionalTags.isEmpty() && keyValueTags.isEmpty();
    }

    /**
     * Retrieves the list of all positional tags associated with this argument.
     * @return A list of positional tags.
     */
    public List<String> getPositionalTags() {
        return Collections.unmodifiableList(positionalTags);
    }

    /**
     * Gets the value associated with a given <key>=<value> argument tag.
     * @param key The key for which to retrieve the value.
     * @return The value paired with the given key, or null if no such element exists.
     */
    public String getValue(final String key) {
        return keyValueTags.get(key);
    }

    /**
     * Returns true if tags contains given key
     * @param key The key for which to check existence.
     * @return true if tags contains given key
     */
    public boolean containsKey(final String key) {
        return keyValueTags.containsKey(key);
    }

    /**
     * Adds positional tag(s) to the tag object.
     * @param tags The tag strings to add.
     */
    protected void addPositionalTag(final String... tags) {
        positionalTags.addAll(Arrays.asList(tags));
    }

    /**
     * Adds a <key>-<value> tag to this tag library.
     * @param key key tag to add.
     * @param value value to associate with this key.
     */
    protected void addKeyValueTag(final String key, final String value) {
        keyValueTags.put(key,value);
    }
}
