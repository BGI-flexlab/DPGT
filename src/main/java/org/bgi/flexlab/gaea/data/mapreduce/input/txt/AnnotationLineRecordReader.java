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
package org.bgi.flexlab.gaea.data.mapreduce.input.txt;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.Seekable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.compress.*;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class AnnotationLineRecordReader extends RecordReader<LongWritable, Text>  {
    private static final Log LOG = LogFactory.getLog(AnnotationLineRecordReader.class);
    public static final String MAX_LINE_LENGTH =
            "mapreduce.input.linerecordreader.line.maxlength";

    private long start;
    private long pos;
    private long end;
    private SplitLineReader in;
    private FSDataInputStream fileIn;
    private Seekable filePosition;
    private int maxLineLength;
    private LongWritable key;
    private Text value;
    private boolean isCompressedInput;
    private Decompressor decompressor;
    private byte[] recordDelimiterBytes;
    private Map<String, Integer> contigIndexs;

    public AnnotationLineRecordReader() {
    }

    public AnnotationLineRecordReader(byte[] recordDelimiter) {
        this.recordDelimiterBytes = recordDelimiter;
    }

    @Override
    public void initialize(InputSplit genericSplit, TaskAttemptContext context) throws IOException, InterruptedException {
        FileSplit split = (FileSplit) genericSplit;
        Configuration job = context.getConfiguration();
        this.maxLineLength = job.getInt(MAX_LINE_LENGTH, Integer.MAX_VALUE);
        start = split.getStart();
        end = start + split.getLength();
        final Path file = split.getPath();

        // open the file and seek to the start of the split
        final FileSystem fs = file.getFileSystem(job);
        fileIn = fs.open(file);

        CompressionCodec codec = new CompressionCodecFactory(job).getCodec(file);
        if (null!=codec) {
            isCompressedInput = true;
            decompressor = CodecPool.getDecompressor(codec);
            if (codec instanceof SplittableCompressionCodec) {
                final SplitCompressionInputStream cIn =
                        ((SplittableCompressionCodec)codec).createInputStream(
                                fileIn, decompressor, start, end,
                                SplittableCompressionCodec.READ_MODE.BYBLOCK);
                in = new CompressedSplitLineReader(cIn, job,
                        this.recordDelimiterBytes);
                start = cIn.getAdjustedStart();
                end = cIn.getAdjustedEnd();
                filePosition = cIn;
            } else {
                in = new SplitLineReader(codec.createInputStream(fileIn,
                        decompressor), job, this.recordDelimiterBytes);
                filePosition = fileIn;
            }
        } else {
            fileIn.seek(start);
            in = new UncompressedSplitLineReader(
                    fileIn, job, this.recordDelimiterBytes, split.getLength());
            filePosition = fileIn;
        }
        // If this is not the first split, we always throw away first record
        // because we always (except the last split) read one extra line in
        // next() method.
        if (start != 0) {
            start += in.readLine(new Text(), 0, maxBytesToConsume(start));
        }
        this.pos = start;

        String[] contigNames = job.getStrings("contigName");
        contigIndexs = new HashMap<>();
        for (int i = 0; i < contigNames.length; i++) {
            contigIndexs.put(contigNames[i], i);
        }

    }

    private int maxBytesToConsume(long pos) {
        return isCompressedInput
                ? Integer.MAX_VALUE
                : (int) Math.max(Math.min(Integer.MAX_VALUE, end - pos), maxLineLength);
    }

    public long getKey(String rec) {
        String[] fields = rec.split("\t", 4);
        long keyout = contigIndexs.get(fields[1]);
        keyout = (keyout << 32) | Integer.parseInt(fields[2]);
        return keyout;
    }

    @Override
    public boolean nextKeyValue() throws IOException, InterruptedException {
        if (key == null) {
            key = new LongWritable();
        }
        if (value == null) {
            value = new Text();
        }
        int newSize = 0;
        // We always read one extra line, which lies outside the upper
        // split limit i.e. (end - 1)
        while (getFilePosition() <= end || in.needAdditionalRecordAfterSplit()) {
            if (pos == 0) {
                newSize = skipUtfByteOrderMark();
            } else {
                newSize = in.readLine(value, maxLineLength, maxBytesToConsume(pos));
                pos += newSize;
            }

            if ((newSize == 0) || (newSize < maxLineLength)) {
                break;
            }

            // line too long. try again
            LOG.info("Skipped line of size " + newSize + " at pos " +
                    (pos - newSize));
        }
        if (newSize == 0) {
            key = null;
            value = null;
            return false;
        } else {
            key.set(getKey(value.toString()));
            return true;
        }
    }

    private int skipUtfByteOrderMark() throws IOException {
        // Strip BOM(Byte Order Mark)
        // Text only support UTF-8, we only need to check UTF-8 BOM
        // (0xEF,0xBB,0xBF) at the start of the text stream.
        int newMaxLineLength = (int) Math.min(3L + (long) maxLineLength,
                Integer.MAX_VALUE);
        int newSize = in.readLine(value, newMaxLineLength, maxBytesToConsume(pos));
        // Even we read 3 extra bytes for the first line,
        // we won't alter existing behavior (no backwards incompat issue).
        // Because the newSize is less than maxLineLength and
        // the number of bytes copied to Text is always no more than newSize.
        // If the return size from readLine is not less than maxLineLength,
        // we will discard the current line and read the next line.
        pos += newSize;
        int textLength = value.getLength();
        byte[] textBytes = value.getBytes();
        if ((textLength >= 3) && (textBytes[0] == (byte)0xEF) &&
                (textBytes[1] == (byte)0xBB) && (textBytes[2] == (byte)0xBF)) {
            // find UTF-8 BOM, strip it.
            LOG.info("Found UTF-8 BOM and skipped it");
            textLength -= 3;
            newSize -= 3;
            if (textLength > 0) {
                // It may work to use the same buffer and not do the copyBytes
                textBytes = value.copyBytes();
                value.set(textBytes, 3, textLength);
            } else {
                value.clear();
            }
        }
        return newSize;
    }

    @Override
    public LongWritable getCurrentKey() {
        return key;
    }

    @Override
    public Text getCurrentValue() {
        return value;
    }

    private long getFilePosition() throws IOException {
        long retVal;
        if (isCompressedInput && null != filePosition) {
            retVal = filePosition.getPos();
        } else {
            retVal = pos;
        }
        return retVal;
    }

    @Override
    public float getProgress() throws IOException, InterruptedException {
        if (start == end) {
            return 0.0f;
        } else {
            return Math.min(1.0f, (getFilePosition() - start) / (float)(end - start));
        }
    }

    public synchronized void close() throws IOException {
        try {
            if (in != null) {
                in.close();
            }
        } finally {
            if (decompressor != null) {
                CodecPool.returnDecompressor(decompressor);
                decompressor = null;
            }
        }
        contigIndexs.clear();
    }
}
