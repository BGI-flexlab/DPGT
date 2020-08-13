package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.util.BlockCompressedInputStream;

public final class XReadLines implements Iterator<String>, Iterable<String>, AutoCloseable {
    private final BufferedReader in;      // The stream we're reading from
    private String nextLine = null;       // Return value of next call to next()
    private final boolean trimWhitespace;
    private final String commentPrefix;

    /**
     * Opens the given file for reading lines.
     * The file may be a text file or a gzipped text file (the distinction is made by the file extension).
     * By default, it will trim whitespaces.
     */
    public XReadLines(final File filename) throws IOException {
        this(filename, true);
    }

    /**
     * Opens the given file for reading lines and optionally trim whitespaces.
     * The file may be a text file or a gzipped text file (the distinction is made by the file extension).
     */
    public XReadLines(final File filename, final boolean trimWhitespace) throws IOException {
        this(makeReaderMaybeGzipped(filename), trimWhitespace, null);
    }
    
    public static Reader makeReaderMaybeGzipped(File file) throws IOException {
        final InputStream in = new BufferedInputStream( new FileInputStream(file));
        return makeReaderMaybeGzipped(in, file.getPath().endsWith(".gz"));
    }

    /**
     * makes a reader for an inputStream wrapping it in an appropriate unzipper if necessary
     * @param zipped is this stream zipped
     */
    public static Reader makeReaderMaybeGzipped(InputStream in, boolean zipped) throws IOException {
        if (zipped) {
            return new InputStreamReader(makeZippedInputStream(in));
        } else {
            return new InputStreamReader(in);
        }
    }
    
    public static InputStream makeZippedInputStream(InputStream in) throws IOException {
        Utils.nonNull(in);
        if (BlockCompressedInputStream.isValidFile(in)) {
                return new BlockCompressedInputStream(in);
        } else {
            return new GZIPInputStream(in);
        }
    }

    /**
     * Creates a new xReadLines object to read lines from an bufferedReader
     *
     * @param reader file name
     * @param trimWhitespace trim whitespace
     * @param commentPrefix prefix for comments or null if no prefix is set
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace, final String commentPrefix) {
        this.in = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
        this.trimWhitespace = trimWhitespace;
        this.commentPrefix = commentPrefix;
        try {
            this.nextLine = readNextLine();
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    /**
     * Reads all of the lines in the file, and returns them as a list of strings
     *
     * @return all of the lines in the file.
     */
    public List<String> readLines() {
        List<String> lines = new LinkedList<>();
        for ( String line : this ) {
            lines.add(line);
        }
        return lines;
    }

    /**
     * I'm an iterator too...
     * @return an iterator
     */
    @Override
    public Iterator<String> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return this.nextLine != null;
    }

    /**
     * Actually reads the next line from the stream, not accessible publicly
     * @return the next line or null
     * @throws IOException if an error occurs
     */
    private String readNextLine() throws IOException {
        String nextLine;
        while ((nextLine = this.in.readLine()) != null) {
            if (this.trimWhitespace) {
                nextLine = nextLine.trim();
                if (nextLine.isEmpty())
                    continue;
            }
            if (this.commentPrefix != null)
                if (nextLine.startsWith(this.commentPrefix))
                    continue;
            break;
        }
        return nextLine;
    }

    /**
     * Returns the next line (optionally minus whitespace)
     * @return the next line
     */
    @Override
    public String next() {
        try {
            String result = this.nextLine;
            this.nextLine = readNextLine();

            // If we haven't reached EOF yet
            if (this.nextLine == null) {
                in.close();             // And close on EOF
            }

            // Return the line we read last time through.
            return result;
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    // The file is read-only; we don't allow lines to be removed.
    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() throws IOException {
        this.in.close();
    }
}