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

import org.bgi.flexlab.gaea.tools.annotator.util.BufferedRandomAccessFile;

import org.bgi.flexlab.gaea.util.Pair;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static java.lang.Integer.min;

public class TSVAdapter extends DBAdapter {

    private BufferedRandomAccessFile tsvReader = null;
    private String filePath = null;
    private HashMap<String ,Pair<Long, Long>> fileIndexChr;
    private long numK;

    private HashMap<String ,ArrayList<Pair<Long, Long>>> fileIndexPos;
    private BufferedRandomAccessFile indexReader = null;
    private List<String> keyList = null;
    public TSVAdapter(String confDir) {
        filePath = confDir;
    }
    @Override
    public void connection(String dbName) throws IOException{

        fileIndexPos = new HashMap<>();
        fileIndexChr = new HashMap<>();
        //chr索引
        BufferedRandomAccessFile br =
                new BufferedRandomAccessFile(new File(filePath + "/" + dbName + "-chr.idx"), "r");
        //读取第一行，pos索引间隔
        numK = Long.valueOf(br.readLine());
        System.out.println("numK:"+numK);
        String line = null;
        while((line = br.readLine())!= null && line.length() > 1)
        {
            String[] arr = line.split("\t");
            fileIndexChr.put(arr[0], new Pair<>(Long.valueOf(arr[1]), Long.valueOf(arr[2])));
        }
        br.close();
        tsvReader = new BufferedRandomAccessFile(new File(filePath + "/" + dbName + ".tsv"), "r");
        keyList = new ArrayList<String>();
        String header = tsvReader.readLine();
        if(!header.contains("#"))
            throw new IOException("header error");
        String[] arr = header.substring(1).split("\t");
        for(String key: arr){
            keyList.add(key);
        }


    }

    @Override
    public void disconnection() throws IOException{
        if(tsvReader != null)
            tsvReader.close();

    }

    @Override
    public HashMap<String, String> getResult(String tableName,
                                             String rowKey) throws IOException{

        HashMap<String,String> resultMap = new HashMap<>();
        if(tableName.contains("vcf"))//查找索引
        {
            String chr = rowKey.split("-")[0];
            long posSearch = Long.valueOf(rowKey.split("-")[1]);
            if(!fileIndexChr.containsKey(chr))
                return null;
            if(!fileIndexPos.containsKey(chr))
            {
                //pos索引
                indexReader = new BufferedRandomAccessFile(filePath + "/" + tableName + "-pos.idx", "r");

                indexReader.seek(fileIndexChr.get(chr).first);

                Long count = fileIndexChr.get(chr).second;
                ArrayList<Pair<Long, Long>> arr = new ArrayList<Pair<Long, Long>>();
                String line = null;

                while( count>0 && (line = indexReader.readLine()) != null) {
                    String[] values = line.split("\t");
                    arr.add(new Pair<>(Long.valueOf(values[0]), Long.valueOf(values[1])));
                    count --;
                }
                fileIndexPos.put(chr ,arr);
                indexReader.close();
            }
            ArrayList<Pair<Long, Long>> array = fileIndexPos.get(chr);
            int i = 0;
            int j = array.size() - 1;
            String value = null;
            long posMin = array.get(i).first;
            long posMax = array.get(j).first;
            long offset = -1;
            //二分查找确定索引范围，再遍历找到索引
            if(posSearch <=  posMax && posSearch >=posMin) {
                while (i <= j) {
                    int mid = (i + j) / 2;
                    long pos = array.get(mid).first;
                    if (pos == posSearch) {
                        offset = array.get(mid).second;
                        break;
                    } else if (pos < posSearch)
                        i = mid + 1;
                    else//(pos > posSearch)
                        j = mid - 1;
                }
                if(offset<0)
                    offset = array.get(j).second;
                indexReader = new BufferedRandomAccessFile(filePath + "/" + tableName  + ".idx", "r");
                indexReader.seek(offset);
                String line = null;
                for(int k = 0; k < numK && (line = indexReader.readLine()) != null; k ++)
                {
                    if(line.startsWith(rowKey))
                    {
                        value = line.split("\t")[1];
                    }
                }
            }

            if(value != null) {
                System.out.println(rowKey + " " + value);
                resultMap.put(rowKey, value);
            }
            indexReader.close();
        }
        else if(tableName.contains("data"))//查找数据
        {
            long pos = Long.valueOf(rowKey);
            tsvReader.seek(pos);
            String line = tsvReader.readLine();
            String[] arr = line.split("\t");
            if(keyList.size() != arr.length)
            {
                throw new IOException("keys num error");
            }
            for(int i = 0; i < min(keyList.size(), arr.length); i ++)
            {
                resultMap.put(keyList.get(i), arr[i]);
            }

        }
        else
            throw  new IOException("table name error!!");

        return resultMap;
    }
}
