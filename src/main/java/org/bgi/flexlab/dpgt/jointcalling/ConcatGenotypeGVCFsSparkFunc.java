/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
package org.bgi.flexlab.dpgt.jointcalling;

import java.util.ArrayList;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.Iterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.spark.api.java.function.Function2;


public class ConcatGenotypeGVCFsSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    private static final Logger logger = LoggerFactory.getLogger(ConcatGenotypeGVCFsSparkFunc.class);
    public String headerPath;
    public String output;
    public long offset;
    public ConcatGenotypeGVCFsSparkFunc(final String headerPath, final String output, long offset) {
        this.headerPath = headerPath;
        this.output = output;
        this.offset = offset;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
        File outputFile = new File(output);
        if (headerPath != null && outputFile.exists()) {
            outputFile.delete();
        }
        BufferedOutputStream bufferedOutputStream = null;
        try {
            FileOutputStream outputStream = new FileOutputStream(outputFile, true);
            FileChannel outputChannel = outputStream.getChannel();
            outputChannel.truncate(offset);
            bufferedOutputStream = new BufferedOutputStream(outputStream);
        } catch (IOException e) {
            logger.error("{}", output, e.getMessage());
            System.exit(1);
        }
        byte[] buffer = new byte[1024*1024*4];
        if (headerPath != null) {
            copyBytes(headerPath, bufferedOutputStream, buffer);
        }
        while(vcfpathIter.hasNext()) {
            String vcfpath = vcfpathIter.next();
            if (vcfpath == null || vcfpath.equals("null")) {
                continue;
            }
            copyBytes(vcfpath, bufferedOutputStream, buffer);
        }
        try {
            bufferedOutputStream.close();
        } catch (Exception e) {
            logger.error("Failed to close {}. {}", output, e.getMessage());
            System.exit(1);
        }
        ArrayList<String> result=new ArrayList<>();
        result.add(output);
        return result.iterator();
    }

    private void copyBytes(final String input, BufferedOutputStream outputStream, byte[] buffer) {
        try {
            BufferedInputStream bufferedInputStream = new BufferedInputStream(new FileInputStream(input));
            int n = 0;
            while ( (n = bufferedInputStream.read(buffer)) > 0) {
                outputStream.write(buffer, 0, n);
            }
            bufferedInputStream.close();
        } catch (IOException e) {
            logger.error("Failed to open {} for read. {}", input, e.getMessage());
            System.exit(1);
        }
    }
}
