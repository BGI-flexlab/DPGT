/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.vcf;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.RuntimeEOFException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFHdfsLoader;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.EnumSet;

public class VCFLocalWriter extends VCFFileWriter {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 3915150472678717360L;
	
	public enum OutputType {
		VCF,
		BLOCK_COMPRESSED_VCF
	}
	
	public final EnumSet<OutputType> FILE_TYPES = EnumSet.allOf(OutputType.class);
		
	public VCFLocalWriter(String filePath, boolean doNotWriteGenotypes,
			final boolean allowMissingFieldsInHeader) throws IOException {
		super(filePath, doNotWriteGenotypes, allowMissingFieldsInHeader, null);
	}
/**
 * htsjdk for now only support BlockCompressedOutputStream for local file system rather than hdfs.
 */
	@Override
	public void initOutputStream(String filePath, Configuration conf) {
		// TODO Auto-generated method stub
		
		OutputType typeTobuild = determineOutputTypeFromFile(filePath);
		try {
			switch (typeTobuild) {
			case VCF:
				os = new FileOutputStream(new File(filePath));
				break;
			case BLOCK_COMPRESSED_VCF:
				os = new BlockCompressedOutputStream(new File(filePath));
			}
		} catch (IOException e) {
			throw new RuntimeEOFException(e);
		}
	}
	
    public static OutputType determineOutputTypeFromFile(final String f) {
        if (isCompressedVCF(f)) {
            return OutputType.BLOCK_COMPRESSED_VCF;
        } else {
            return OutputType.VCF;
        }
    }
    
    private static boolean isCompressedVCF(final String outFile) {
        if (outFile == null)
            return false;

        return VCFHdfsLoader.hasBlockCompressedExtension(outFile);
    }
}
