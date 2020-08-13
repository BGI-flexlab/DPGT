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
 * Copyright (c) 2010 Aalto University 
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
package org.bgi.flexlab.gaea.util;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;
import org.apache.hadoop.mapreduce.lib.partition.TotalOrderPartitioner;
import org.bgi.flexlab.gaea.data.structure.header.MultipleVCFHeader;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.BioJob;
import org.bgi.flexlab.gaea.tools.mapreduce.vcf.sort.VCFSortOptions;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;

public class SortUilts {
	
	public static void configureSampling(Path outPath, BioJob job, VCFSortOptions options) throws IOException{
    	Configuration conf = job.getConfiguration();
		final Path partition = outPath.getFileSystem(conf).makeQualified(new Path(outPath, "_partitioning" + "VCF"));
    	
    	TotalOrderPartitioner.setPartitionFile(conf, partition);
    	try {
    		final URI partitionURI = new URI(partition.toString() + "#" + partition.getName());
    		
    		if(partitionURI.getScheme().equals("file"))
    			return;
    		job.addCacheFile(partitionURI);
//    		ReferenceShare.distributeCache(partitionURI.toString(), job);
    	} catch (URISyntaxException e) { throw new RuntimeException(e); }
    	
	}
	
	public static void merge(MultipleVCFHeader mVcfHeader, VCFSortOptions options, Configuration conf) {
		try {
			System.out.println("vcf-MultiSampleSort :: Merging output...");

			// First, place the VCF or BCF header.
            final Path outpath = new Path(options.getOutputPath());
			final Path wrkPath = new Path(options.getWorkPath());
			final FileSystem srcFS = wrkPath.getFileSystem(conf);
			final FileSystem dstFS = outpath.getFileSystem(conf);

			Map<String, OutputStream> outs = new HashMap<String, OutputStream>();
			Map<Integer, String> multiOutputs = options.getMultiOutputs();
			for(String result : multiOutputs.values()){
				Path sPath = new Path(options.getOutputPath() + "/" + result + ".vcf");
				OutputStream os = dstFS.create(sPath);
				outs.put(result, os);
			}
			
		    final VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
			VariantContextWriter writer;
			Map<Integer, SingleVCFHeader> id2VcfHeader = mVcfHeader.getID2SingleVcfHeader();
            for( int id : multiOutputs.keySet()){
            	VCFHeader newHeader = id2VcfHeader.get(id).getHeader();
            	writer = builder.setOutputStream(
            			new FilterOutputStream(outs.get(multiOutputs.get(id))) {
            				@Override public void close() throws IOException {
            					this.out.flush();
            				}
            			}).setOptions(VariantContextWriterBuilder.NO_OPTIONS).build();
            	
            	writer.writeHeader(newHeader);
            	writer.close();
            	
            	final FileStatus[] parts = srcFS.globStatus(new Path(options.getWorkPath(), multiOutputs.get(id) + "-*-[0-9][0-9][0-9][0-9][0-9]*"));
            	
            	int i = 0;
            	
            	for(final FileStatus part : parts){
            		System.out.printf("sort:: Merging part %d ( size %d)...\n", i++, part.getLen());
            		System.out.flush();
            		
            		final FSDataInputStream ins = srcFS.open(part.getPath());
            		IOUtils.copyBytes(ins, outs.get(multiOutputs.get(id)), conf, false);
            		ins.close();            		
            	}
               for (final FileStatus part : parts)
					srcFS.delete(part.getPath(), false);  			

			    outs.get(multiOutputs.get(id)).close();

            }
		} catch (IOException e) {
			System.err.printf("vcf-MultiSampleSort :: Output merging failed: %s\n", e);
		}
	}

}
