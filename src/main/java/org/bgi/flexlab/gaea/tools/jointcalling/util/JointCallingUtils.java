package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.function.Supplier;

import org.bgi.flexlab.gaea.util.Utils;

public class JointCallingUtils extends Utils{
	public static <T> T nonNull(final T object) {
        return JointCallingUtils.nonNull(object, "Null object is not allowed here.");
    }
	
	public static <T> T nonNull(final T object, final String message) {
        if (object == null) {
            throw new IllegalArgumentException(message);
        }
        return object;
    }
	
	public static <T> T nonNull(final T object, final Supplier<String> message) {
        if (object == null) {
            throw new IllegalArgumentException(message.get());
        }
        return object;
    }

    public static void validateArg(final boolean condition, final String msg){
        if (!condition){
            throw new IllegalArgumentException(msg);
        }
    }

    public static void validateArg(final boolean condition, final Supplier<String> msg){
        if (!condition){
            throw new IllegalArgumentException(msg.get());
        }
    }
}
