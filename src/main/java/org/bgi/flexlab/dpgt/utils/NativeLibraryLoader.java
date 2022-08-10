package org.bgi.flexlab.dpgt.utils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public final class NativeLibraryLoader {
    private static final Logger logger = LoggerFactory.getLogger(NativeLibraryLoader.class);
    private static final String LIB_NAME = "cdpgt";
    public static boolean isLoad = false;
    public static boolean isLocal = false;
    
    public static synchronized void load() {
        if (!isLoad || !isLocal) {
            try{
                System.loadLibrary(LIB_NAME);
            }catch(UnsatisfiedLinkError e){
                logger.error("Failed to load Native code library {}. {}", LIB_NAME, e.getMessage());
                System.exit(1);
            }
            isLoad = true;
        }
    }
}
