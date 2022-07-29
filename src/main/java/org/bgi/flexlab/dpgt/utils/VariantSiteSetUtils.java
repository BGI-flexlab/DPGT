package org.bgi.flexlab.dpgt.utils;

import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.BufferedOutputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import org.broadinstitute.hellbender.utils.SimpleInterval;

public class VariantSiteSetUtils {
    private static final Logger logger = LoggerFactory.getLogger(VariantSiteSetUtils.class);

    public static BitSet loadVariantSiteSet(final String input) {
        try {
            BufferedInputStream bufferedInputStream = new BufferedInputStream(new FileInputStream(input));
            byte[] sizeBuff = new byte[4];  // buffer for read size of the bit set file in little-endian order
            int size = 0;
            int n = 0;
            if ( (n = bufferedInputStream.read(sizeBuff)) == 4 )
            {
                // bytes to int, from little-endian(c/c++) to big-endian(java)
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
            logger.error("Failed to read variant site file {}, {}.", input, e.getMessage());
            System.exit(1);
        }
        return null;
    }

    public static void writeVariantSiteSet(final BitSet bitSet, final String output) {
        try {
            byte[] bytes = bitSet.toByteArray();
            BufferedOutputStream outputStream = new BufferedOutputStream(new FileOutputStream(output));
            byte[] sizeBuff = new byte[4];  // buffer for write size of the bit set in little-endian order
            int size = bytes.length;
            sizeBuff[0] = (byte) (size & 0xFF);
            sizeBuff[1] = (byte) ((size >> 8) & 0xFF);
            sizeBuff[2] = (byte) ((size >> 16) & 0xFF);
            sizeBuff[3] = (byte) ((size >> 24) & 0xFF);
            outputStream.write(sizeBuff);
            outputStream.write(bytes);
            outputStream.close();
        } catch (IOException e) {
            logger.error("Failed to write variant site file {}, {}.", output, e.getMessage());
            System.exit(1);
        }
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

    public static ArrayList<BitSet> getSubBitSets(final BitSet largeBitSet, final SimpleInterval largeInterval, final List<SimpleInterval> subIntervals)
    {
        ArrayList<BitSet> result = new ArrayList<>();
        for (SimpleInterval interval: subIntervals) {
            int i = interval.getStart() - largeInterval.getStart();
            int j = interval.getEnd() - largeInterval.getStart() + 1;
            result.add(largeBitSet.get(i, j));
        }
        return result;
    }
}
