package org.bgi.flexlab.dpgt.utils;

import java.util.BitSet;
import java.util.Iterator;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.BufferedInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class VariantSiteSetUtils {
    private static final Logger logger = LoggerFactory.getLogger(VariantSiteSetUtils.class);

    public static BitSet loadVariantSiteSet(final String input) {
        try {
            BufferedInputStream bufferedInputStream = new BufferedInputStream(new FileInputStream(input));
            byte[] sizeBuff = new byte[4];  // buffer for read size of the bit set file
            int size = 0;
            int n = 0;
            if ( (n = bufferedInputStream.read(sizeBuff)) == 4 )
            {
                // bytes to long, from little-endian(c/c++) to big-endian(java)
                size = ((sizeBuff[0]&0xFF) | ((sizeBuff[1]&0xFF)<<8) | ((sizeBuff[2]&0xFF)<<16) | ((sizeBuff[3]&0xFF)<<24));
            } else {
                bufferedInputStream.close();
                logger.error("Failed to read bitset bytes size from {}, file may be truncated.", input);
                System.exit(1);
            }

            byte[] bytes = new byte[size];
            if ( (n = bufferedInputStream.read(bytes)) == size)
            {
                bufferedInputStream.close();
                return BitSet.valueOf(bytes);
            } else {
                bufferedInputStream.close();
                logger.error("Failed to read bytes from {}, expected {} bytes but read {} bytes, file may be truncated.", input, size, n);
                System.exit(1);
            }
        } catch (IOException e) {
            logger.error("Failed open variant site file {}, {}.", input, e.getMessage());
            System.exit(1);
        }
        return null;
    }

    public static BitSet combine(Iterator<String> variantSiteFilesIter) {
        if (!variantSiteFilesIter.hasNext()) {
            return null;
        }
        BitSet result = VariantSiteSetUtils.loadVariantSiteSet(variantSiteFilesIter.next());
        while (variantSiteFilesIter.hasNext()) {
            result.or(VariantSiteSetUtils.loadVariantSiteSet(variantSiteFilesIter.next()));
        }
        return result;
    }

}
