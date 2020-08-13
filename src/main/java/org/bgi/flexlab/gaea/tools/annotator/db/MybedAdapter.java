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

import htsjdk.tribble.readers.TabixReader;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

public class MybedAdapter extends DBAdapter {

    private TabixReader bedReader;
    private String filePath = null;
    private String[] header = null;
    public MybedAdapter(String confDir) {
        filePath = confDir;
    }

    @Override
    public void connection(String tableName) throws IOException{
        bedReader = new TabixReader(filePath + "/" + tableName);

        InputStream in = new GZIPInputStream(new FileInputStream(filePath + "/" + tableName));
        Scanner sc=new Scanner(in);
        while(sc.hasNextLine()){
            String line = sc.nextLine();
            if(line.startsWith("##"))
                continue;
            if(line.startsWith("#")){
                header = line.substring(1).split("\t");
                break;
            }
            if(!line.startsWith("#"))
                throw new IOException("header error");
        }
        sc.close();
    }

    @Override
    public void disconnection() throws IOException{
        if(bedReader != null)
            bedReader.close();
    }

    @Override
    public HashMap<String, String> getResult(String tableName, String rowKey, List<String> fieldMap) throws IOException{
        return null;
    }

    @Override
    public HashMap<String, String> getResult(String tableName,
                                             String reg) throws IOException{

        HashMap<String,String> resultMap = new HashMap<>();
        String[] region = reg.split("-");

        String[] fields = null;
        String s;
        TabixReader.Iterator iter = bedReader.query(region[0],  Integer.valueOf(region[1]), Integer.valueOf(region[2]));
        while (iter != null && (s = iter.next()) != null){
            if(fields == null)
                fields = s.split("\t");
            else {
                String[] t = s.split("\t");
                for (int i = 3; i < header.length; i++) {
                    if(!fields[i].equals(t[i]))
                        fields[i] = fields[i] + "," + t[i];
                }
            }

        }

        if(fields == null)
            return null;

        for (int i = 3; i < header.length; i++) {
            resultMap.put(header[i], fields[i]);
        }

        return resultMap;

    }
}
