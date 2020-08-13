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
package org.bgi.flexlab.gaea.data.mapreduce.writable;


import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

public class VcfLineWritable implements Writable {
    private Text fileName;
    private Text vcfLine;


    public VcfLineWritable()
    {
        fileName = new Text();
        vcfLine = new Text();
    }

    public void set(String fileName, String vcfLine)
    {
        this.fileName.set(fileName);
        this.vcfLine.set(vcfLine);
    }


    public String getFileName()
    {
        return fileName.toString();
    }
    public String getVCFLine()
    {
        return vcfLine.toString();
    }

    @Override
    public void write(DataOutput out) throws IOException {
        fileName.write(out);
        vcfLine.write(out);
    }
    @Override
    public void readFields(DataInput in) throws IOException {
        fileName.readFields(in);
        vcfLine.readFields(in);
    }
}
