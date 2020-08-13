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
package org.bgi.flexlab.gaea.tools.annotator.db;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class VCFAdapter extends DBAdapter {

    private VCFFileReader vcfReader = null;
    private String filepath = null;

    public VCFAdapter(String confDir) {
        filepath = confDir;
    }

    @Override
    public void connection(String tableName) throws IOException{
        String fileName = filepath + "/" + tableName;
        vcfReader = new VCFFileReader(new File(fileName));
    }

    @Override
    public void disconnection() throws IOException{
        if(vcfReader != null)
            vcfReader.close();
    }

    public List<HashMap<String, String>> getResult(String reg, List<String> fields) throws IOException{
        List<HashMap<String, String>> results = new ArrayList<>();
        String[] arr = reg.split("\t");
        String chr = arr[0];
        int start = Integer.valueOf(arr[1]);
        int end = Integer.valueOf(arr[2]);
        CloseableIterator<VariantContext> vcfIter = vcfReader.query(chr, start, end);

        while (vcfIter.hasNext()) {
            VariantContext vc = vcfIter.next();
            if (start != vc.getStart() || end != vc.getEnd())
                continue;

            HashMap<String,String> resultMap = new HashMap<>();
            for(String field: fields){
                String v = vc.getAttributeAsString(field, ".");
                if (v.startsWith("[") && v.endsWith("]"))
                    v = v.substring(1, v.length() - 1);
                resultMap.put(field, v);
            }

            List<String> alts = new ArrayList<>();
            for(Allele allele: vc.getAlternateAlleles()){
                alts.add(allele.getBaseString());
            }
            resultMap.put("ID", vc.getID());
            resultMap.put("ALT", String.join(",", alts));
            results.add(resultMap);
        }
        return results;
    }

    @Override
    public HashMap<String, String> getResult(String tableName,
                                             String rowKey) throws IOException{
        HashMap<String,String> resultMap = new HashMap<>();

        String[] arr = rowKey.split("\t");
        String chr = arr[0];
        int start = Integer.valueOf(arr[1]);
        int end = Integer.valueOf(arr[2]);
        //直接根据chr-start-end查找数据
        CloseableIterator<VariantContext> vcfIter = vcfReader.query(chr, start, end);

        while (vcfIter.hasNext()) {
            VariantContext vc = vcfIter.next();
            if (start == vc.getStart() && end == vc.getEnd()) {//输出start，end完全一致的数据
                Map<String, Object> mp = vc.getAttributes();
                for (Map.Entry<String, Object> entry : mp.entrySet()) {
                    resultMap.put(entry.getKey(), entry.getValue().toString());
                }

                resultMap.put("POS", String.valueOf(vc.getStart()));
                resultMap.put("ID", vc.getID());
                resultMap.put("REF", vc.getReference().getDisplayString());
                StringBuilder altStr = new StringBuilder();
                for(Allele allele: vc.getAlternateAlleles()){
                    altStr.append(allele.getDisplayString()+ ",");
                }
                altStr.deleteCharAt(altStr.length() - 1);
                resultMap.put("ALT", altStr.toString());
            }
        }
        return resultMap;
    }
}
