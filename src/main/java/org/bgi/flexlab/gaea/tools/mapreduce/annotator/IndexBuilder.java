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
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.bgi.flexlab.gaea.tools.annotator.util.BufferedRandomAccessFile;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by 黄凯文 on 2017/8/7.
 */
public class IndexBuilder {
    static public String getChromeNoChr(VariantContext context){
        if(context.getContig().startsWith("chr")){
            return context.getContig().substring(3);
        }
        return context.getContig();
    }

    public static void main(String[] args) throws IOException{
        final int K = Integer.valueOf(args[1]);//pos索引间隔
        VCFFileReader vcfReader = new VCFFileReader(new File(args[0]));

        CloseableIterator<VariantContext> vcfIter =  vcfReader.iterator();

        Set<String> keys = new HashSet<String>();
        while(vcfIter.hasNext()){
            VariantContext vc = vcfIter.next();
            for(Map.Entry<String, Object> entry: vc.getAttributes().entrySet()){
                keys.add(entry.getKey());
            }
        }
        vcfIter.close();

        BufferedRandomAccessFile contentWriter = new BufferedRandomAccessFile(args[0] + ".tsv", "rw");
        BufferedRandomAccessFile indexWriter = new BufferedRandomAccessFile(args[0] + ".idx", "rw");

        StringBuffer sb = new StringBuffer();
        sb.append("#CHR\tPOS\tREF\tALT");
        for(String key: keys){
            sb.append("\t"+key);
        }
        sb.append("\r\n");
        contentWriter.write(sb.toString().getBytes());

        vcfIter =  vcfReader.iterator();
        long textSize = contentWriter.getFilePointer();//sb.length()+1;
        List<String> indexList = new ArrayList<>();

        while(vcfIter.hasNext()) {
            VariantContext vc = vcfIter.next();
            for(Allele allele: vc.getAlleles()){
                StringBuffer contentLine = new StringBuffer();

                String chr = getChromeNoChr(vc);
                int pos = vc.getStart();
                String ref = vc.getReference().getDisplayString();
                String alt = allele.getDisplayString();

                if(ref.compareTo(alt) != 0) {
                    contentLine.append(chr + "\t" + pos + "\t" + ref + "\t" + alt);
                    String indexKey = chr + "-" + pos + "-" + vc.getEnd();
                    String indexValue = String.valueOf(textSize);

                    if(!indexList.isEmpty() && indexList.get(indexList.size() - 1).startsWith(indexKey))
                    {
                        indexList.set(indexList.size() - 1, indexList.get(indexList.size() - 1) + "," + indexValue);
                    }
                    else
                    {
                        indexList.add(indexKey + "\t" + indexValue);

                    }
                    for (String key : keys) {
                        if (vc.hasAttribute(key)) {
                            contentLine.append("\t" + String.valueOf(vc.getAttribute(key)));
                        } else {
                            contentLine.append("\t.");
                        }
                    }
                    contentLine.append("\r\n");
                    contentWriter.write(contentLine.toString().getBytes());

                    textSize = contentWriter.getFilePointer();
                }
            }
        }


        textSize = 0;
        String preChr = "";
        long num = 0;
        long numPos = 0;
        long preTextSize = 0;
        String chr = null;
        BufferedRandomAccessFile indexPos = new BufferedRandomAccessFile(args[0] + "-pos.idx", "rw");
        BufferedRandomAccessFile indexChr = new BufferedRandomAccessFile(args[0] + "-chr.idx", "rw");
        indexChr.write((K+"\r\n").getBytes());
        for(int i = 0; i < indexList.size(); i++){
            String str = indexList.get(i);
            chr = str.split("-")[0];
            String pos = str.split("-")[1];
            if(preChr.compareTo(chr) != 0)
            {
                if(!preChr.equals(""))
                    indexChr.write((preChr + "\t" + preTextSize+"\t" + numPos+"\r\n").getBytes());

                num = 0;
                numPos = 0;
                preChr = chr;
                preTextSize = indexPos.getFilePointer();

            }
            if(num % K == 0 || (i < indexList.size() - 1 && indexList.get(i + 1).startsWith(chr+"-"+pos)))
            {
                if(pos.contains("-"))
                    System.err.println(pos);
                indexPos.write((pos + "\t" + textSize + "\r\n").getBytes());
                numPos++;
            }
            indexWriter.write((str+"\r\n").getBytes());
            textSize = indexWriter.getFilePointer();
            num ++;
        }
        indexChr.write((preChr + "\t" + preTextSize+"\t" + numPos).getBytes());
        indexChr.close();
        indexPos.close();

        indexWriter.close();
        contentWriter.close();
    }
}
