package org.bgi.flexlab.dpgt.jointcalling;


public class NativeInitializer {
    /**
     * Initialize JNI C/C++ environment.
     * 1. check if jemalloc is in use.
     * 2. initialize spdlog logger
     */
    public native void apply();
}
