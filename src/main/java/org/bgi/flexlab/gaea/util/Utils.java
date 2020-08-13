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

package org.bgi.flexlab.gaea.util;

import htsjdk.samtools.util.StringUtil;
import org.apache.log4j.Logger;

import java.lang.reflect.Array;
import java.net.InetAddress;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Utils {
	/** our log, which we want to capture anything from this class */
	private static Logger logger = Logger.getLogger(Utils.class);

	public static final float JAVA_DEFAULT_HASH_LOAD_FACTOR = 0.75f;

	/**
	 * Calculates the optimum initial size for a hash table given the maximum
	 * number of elements it will need to hold. The optimum size is the smallest
	 * size that is guaranteed not to result in any rehash/table-resize
	 * operations.
	 *
	 * @param maxElements
	 *            The maximum number of elements you expect the hash table will
	 *            need to hold
	 * @return The optimum initial size for the table, given maxElements
	 */
	public static int optimumHashSize(int maxElements) {
		return (int) (maxElements / JAVA_DEFAULT_HASH_LOAD_FACTOR) + 2;
	}

	/**
	 * Compares two objects, either of which might be null.
	 *
	 * @param lhs
	 *            One object to compare.
	 * @param rhs
	 *            The other object to compare.
	 *
	 * @return True if the two objects are equal, false otherwise.
	 */
	public static boolean equals(Object lhs, Object rhs) {
		if (lhs == null && rhs == null)
			return true;
		else if (lhs == null)
			return false;
		else
			return lhs.equals(rhs);
	}

	public static <T> T nonNull(final T object) {
		return Utils.nonNull(object, "Null object is not allowed here.");
	}

	/**
	 * Checks that an {@link Object} is not {@code null} and returns the same
	 * object or throws an {@link IllegalArgumentException}
	 * 
	 * @param object
	 *            any Object
	 * @param message
	 *            the text message that would be passed to the exception thrown
	 *            when {@code o == null}.
	 * @return the same object
	 * @throws IllegalArgumentException
	 *             if a {@code o == null}
	 */
	public static <T> T nonNull(final T object, final String message) {
		if (object == null) {
			throw new IllegalArgumentException(message);
		}
		return object;
	}

	public static <T> List<T> cons(final T elt, final List<T> l) {
		List<T> l2 = new ArrayList<T>();
		l2.add(elt);
		if (l != null)
			l2.addAll(l);
		return l2;
	}
	
	public static <I, T extends Collection<I>> T nonEmpty(T collection, String message){
	    nonNull(collection, "The collection is null: " + message);
	    if(collection.isEmpty()){
	        throw new IllegalArgumentException("The collection is empty: " + message);
	    } else {
	        return collection;
	    }
	}

	public static void warnUser(final String msg) {
		warnUser(logger, msg);
	}

	public static void warnUser(final Logger logger, final String msg) {
		logger.warn(String.format("********************************************************************************"));
		logger.warn(String.format("* WARNING:"));
		logger.warn(String.format("*"));
		prettyPrintWarningMessage(logger, msg);
		logger.warn(String.format("********************************************************************************"));
	}

	/**
	 * pretty print the warning message supplied
	 *
	 * @param logger
	 *            logger for the message
	 * @param message
	 *            the message
	 */
	private static void prettyPrintWarningMessage(Logger logger, String message) {
		StringBuilder builder = new StringBuilder(message);
		while (builder.length() > 70) {
			int space = builder.lastIndexOf(" ", 70);
			if (space <= 0)
				space = 70;
			logger.warn(String.format("* %s", builder.substring(0, space)));
			builder.delete(0, space + 1);
		}
		logger.warn(String.format("* %s", builder));
	}

	public static ArrayList<Byte> subseq(char[] fullArray) {
		byte[] fullByteArray = new byte[fullArray.length];
		StringUtil.charsToBytes(fullArray, 0, fullArray.length, fullByteArray, 0);
		return subseq(fullByteArray);
	}

	public static ArrayList<Byte> subseq(byte[] fullArray) {
		return subseq(fullArray, 0, fullArray.length - 1);
	}

	public static ArrayList<Byte> subseq(byte[] fullArray, int start, int end) {
		assert end < fullArray.length;
		ArrayList<Byte> dest = new ArrayList<Byte>(end - start + 1);
		for (int i = start; i <= end; i++) {
			dest.add(fullArray[i]);
		}
		return dest;
	}

	public static String baseList2string(List<Byte> bases) {
		byte[] basesAsbytes = new byte[bases.size()];
		int i = 0;
		for (Byte b : bases) {
			basesAsbytes[i] = b;
			i++;
		}
		return new String(basesAsbytes);
	}

	/**
	 * join the key value pairs of a map into one string, i.e. myMap =
	 * [A->1,B->2,C->3] with a call of: joinMap("-","*",myMap) -> returns
	 * A-1*B-2*C-3
	 *
	 * Be forewarned, if you're not using a map that is aware of the ordering
	 * (i.e. HashMap instead of LinkedHashMap) the ordering of the string you
	 * get back might not be what you expect! (i.e. C-3*A-1*B-2 vrs A-1*B-2*C-3)
	 *
	 * @param keyValueSeperator
	 *            the string to seperate the key-value pairs
	 * @param recordSeperator
	 *            the string to use to seperate each key-value pair from other
	 *            key-value pairs
	 * @param map
	 *            the map to draw from
	 * @param <L>
	 *            the map's key type
	 * @param <R>
	 *            the map's value type
	 * @return a string representing the joined map
	 */
	public static <L, R> String joinMap(String keyValueSeperator, String recordSeperator, Map<L, R> map) {
		if (map.size() < 1) {
			return null;
		}
		String joinedKeyValues[] = new String[map.size()];
		int index = 0;
		for (L key : map.keySet()) {
			joinedKeyValues[index++] = String.format("%s%s%s", key.toString(), keyValueSeperator,
					map.get(key).toString());
		}
		return join(recordSeperator, joinedKeyValues);
	}

	/**
	 * Splits a String using indexOf instead of regex to speed things up.
	 *
	 * @param str
	 *            the string to split.
	 * @param delimiter
	 *            the delimiter used to split the string.
	 * @return an array of tokens.
	 */
	public static ArrayList<String> split(String str, String delimiter) {
		return split(str, delimiter, 10);
	}

	/**
	 * Splits a String using indexOf instead of regex to speed things up.
	 *
	 * @param str
	 *            the string to split.
	 * @param delimiter
	 *            the delimiter used to split the string.
	 * @param expectedNumTokens
	 *            The number of tokens expected. This is used to initialize the
	 *            ArrayList.
	 * @return an array of tokens.
	 */
	public static ArrayList<String> split(String str, String delimiter, int expectedNumTokens) {
		final ArrayList<String> result = new ArrayList<String>(expectedNumTokens);

		int delimiterIdx = -1;
		do {
			final int tokenStartIdx = delimiterIdx + 1;
			delimiterIdx = str.indexOf(delimiter, tokenStartIdx);
			final String token = (delimiterIdx != -1 ? str.substring(tokenStartIdx, delimiterIdx)
					: str.substring(tokenStartIdx));
			result.add(token);
		} while (delimiterIdx != -1);

		return result;
	}

	/**
	 * join an array of strings given a seperator
	 * 
	 * @param separator
	 *            the string to insert between each array element
	 * @param strings
	 *            the array of strings
	 * @return a string, which is the joining of all array values with the
	 *         separator
	 */
	public static String join(String separator, String[] strings) {
		return join(separator, strings, 0, strings.length);
	}

	public static String join(String separator, String[] strings, int start, int end) {
		if ((end - start) == 0) {
			return "";
		}
		StringBuilder ret = new StringBuilder(strings[start]);
		for (int i = start + 1; i < end; ++i) {
			ret.append(separator);
			ret.append(strings[i]);
		}
		return ret.toString();
	}

	public static String join(String separator, int[] ints) {
		if (ints == null || ints.length == 0)
			return "";
		else {
			StringBuilder ret = new StringBuilder();
			ret.append(ints[0]);
			for (int i = 1; i < ints.length; ++i) {
				ret.append(separator);
				ret.append(ints[i]);
			}
			return ret.toString();
		}
	}

	public static String join(String separator, double[] doubles) {
		if (doubles == null || doubles.length == 0)
			return "";
		else {
			StringBuilder ret = new StringBuilder();
			ret.append(doubles[0]);
			for (int i = 1; i < doubles.length; ++i) {
				ret.append(separator);
				ret.append(doubles[i]);
			}
			return ret.toString();
		}
	}

	/**
	 * Returns a string of the form elt1.toString() [sep elt2.toString() ... sep
	 * elt.toString()] for a collection of elti objects (note there's no actual
	 * space between sep and the elti elements). Returns "" if collection is
	 * empty. If collection contains just elt, then returns elt.toString()
	 *
	 * @param separator
	 *            the string to use to separate objects
	 * @param objects
	 *            a collection of objects. the element order is defined by the
	 *            iterator over objects
	 * @param <T>
	 *            the type of the objects
	 * @return a non-null string
	 */
	public static <T> String join(final String separator, final Collection<T> objects) {
		if (objects.isEmpty()) { // fast path for empty collection
			return "";
		} else {
			final Iterator<T> iter = objects.iterator();
			final T first = iter.next();

			if (!iter.hasNext()) // fast path for singleton collections
				return first.toString();
			else { // full path for 2+ collection that actually need a join
				final StringBuilder ret = new StringBuilder(first.toString());
				while (iter.hasNext()) {
					ret.append(separator);
					ret.append(iter.next().toString());
				}
				return ret.toString();
			}
		}
	}

	public static String dupString(char c, int nCopies) {
		char[] chars = new char[nCopies];
		Arrays.fill(chars, c);
		return new String(chars);
	}

	public static byte[] dupBytes(byte b, int nCopies) {
		byte[] bytes = new byte[nCopies];
		Arrays.fill(bytes, b);
		return bytes;
	}

	// trim a string for the given character (i.e. not just whitespace)
	public static String trim(String str, char ch) {
		char[] array = str.toCharArray();

		int start = 0;
		while (start < array.length && array[start] == ch)
			start++;

		int end = array.length - 1;
		while (end > start && array[end] == ch)
			end--;

		return str.substring(start, end + 1);
	}

	public static byte listMaxByte(List<Byte> quals) {
		if (quals.size() == 0)
			return 0;
		byte m = quals.get(0);
		for (byte b : quals) {
			m = b > m ? b : m;
		}
		return m;
	}

	/**
	 * Splits expressions in command args by spaces and returns the array of
	 * expressions. Expressions may use single or double quotes to group any
	 * individual expression, but not both.
	 * 
	 * @param args
	 *            Arguments to parse.
	 * @return Parsed expressions.
	 */
	public static String[] escapeExpressions(String args) {
		// special case for ' and " so we can allow expressions
		if (args.indexOf('\'') != -1)
			return escapeExpressions(args, "'");
		else if (args.indexOf('\"') != -1)
			return escapeExpressions(args, "\"");
		else
			return args.trim().split(" +");
	}

	/**
	 * Splits expressions in command args by spaces and the supplied delimiter
	 * and returns the array of expressions.
	 * 
	 * @param args
	 *            Arguments to parse.
	 * @param delimiter
	 *            Delimiter for grouping expressions.
	 * @return Parsed expressions.
	 */
	private static String[] escapeExpressions(String args, String delimiter) {
		String[] command = {};
		String[] split = args.split(delimiter);
		String arg;
		for (int i = 0; i < split.length - 1; i += 2) {
			arg = split[i].trim();
			if (arg.length() > 0) // if the unescaped arg has a size
				command = Utils.concatArrays(command, arg.split(" +"));
			command = Utils.concatArrays(command, new String[] { split[i + 1] });
		}
		arg = split[split.length - 1].trim();
		if (split.length % 2 == 1) // if the command ends with a delimiter
			if (arg.length() > 0) // if the last unescaped arg has a size
				command = Utils.concatArrays(command, arg.split(" +"));
		return command;
	}

	/**
	 * Concatenates two String arrays.
	 * 
	 * @param A
	 *            First array.
	 * @param B
	 *            Second array.
	 * @return Concatenation of A then B.
	 */
	public static String[] concatArrays(String[] A, String[] B) {
		String[] C = new String[A.length + B.length];
		System.arraycopy(A, 0, C, 0, A.length);
		System.arraycopy(B, 0, C, A.length, B.length);
		return C;
	}

	/**
	 * Appends String(s) B to array A.
	 * 
	 * @param A
	 *            First array.
	 * @param B
	 *            Strings to append.
	 * @return A with B(s) appended.
	 */
	public static String[] appendArray(String[] A, String... B) {
		return concatArrays(A, B);
	}

	/**
	 * Returns indices of all occurrences of the specified symbol in the string
	 * 
	 * @param s
	 *            Search string
	 * @param ch
	 *            Character to search for
	 * @return Indices of all occurrences of the specified symbol
	 */
	public static int[] indexOfAll(String s, int ch) {
		int[] pos = new int[64];
		int z = 0;

		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == ch)
				pos[z++] = i;
		}
		return reallocate(pos, z);
	}

	public static int countSetBits(boolean[] array) {
		int counter = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i])
				counter++;
		}
		return counter;
	}

	/**
	 * Returns new (reallocated) integer array of the specified size, with
	 * content of the original array <code>orig</code> copied into it. If
	 * <code>newSize</code> is less than the size of the original array, only
	 * first <code>newSize</code> elements will be copied. If new size is
	 * greater than the size of the original array, the content of the original
	 * array will be padded with zeros up to the new size. Finally, if new size
	 * is the same as original size, no memory reallocation will be performed
	 * and the original array will be returned instead.
	 *
	 * @param orig
	 *            Original size.
	 * @param newSize
	 *            New Size.
	 *
	 * @return New array with length equal to newSize.
	 */
	public static int[] reallocate(int[] orig, int newSize) {
		if (orig.length == newSize)
			return orig;
		int[] new_array = new int[newSize];
		int L = (newSize > orig.length ? orig.length : newSize);
		for (int i = 0; i < L; i++)
			new_array[i] = orig[i];
		return new_array;
	}

	/**
	 * Returns a copy of array a, extended with additional n elements to the
	 * right (if n > 0 ) or -n elements to the left (if n<0), copying the values
	 * form the original array. Newly added elements are filled with value v.
	 * Note that if array a is being padded to the left, first (-n) elements of
	 * the returned array are v's, followed by the content of array a.
	 * 
	 * @param a
	 *            original array
	 * @param n
	 *            number of (v-filled) elements to append to a on the right
	 *            (n>0) or on the left (n<0)
	 * @param v
	 *            element value
	 * @return the extended copy of array a with additional n elements
	 */
	public static byte[] extend(final byte[] a, int n, byte v) {

		byte[] newA;

		if (n > 0) {
			newA = Arrays.copyOf(a, a.length + n);
			if (v != 0) { // java pads with 0's for us, so there is nothing to
							// do if v==0
				for (int i = a.length; i < newA.length; i++)
					newA[i] = v;
			}
			return newA;
		}

		// we are here only if n < 0:
		n = (-n);
		newA = new byte[a.length + n];
		int i;
		if (v != 0) {
			i = 0;
			for (; i < n; i++)
				newA[i] = v;
		} else {
			i = n;
		}
		for (int j = 0; j < a.length; i++, j++)
			newA[i] = a[j];
		return newA;
	}

	/**
	 * Returns a copy of array a, extended with additional n elements to the
	 * right (if n > 0 ) or -n elements to the left (if n<0), copying the values
	 * form the original array. Newly added elements are filled with value v.
	 * Note that if array a is padded to the left, first (-n) elements of the
	 * returned array are v's, followed by the content of array a.
	 * 
	 * @param a
	 *            original array
	 * @param n
	 *            number of (v-filled) elements to append to a on the right
	 *            (n>0) or on the left (n<0)
	 * @param v
	 *            element value
	 * @return the extended copy of array a with additional n elements
	 */
	public static short[] extend(final short[] a, int n, short v) {

		short[] newA;

		if (n > 0) {
			newA = Arrays.copyOf(a, a.length + n);
			if (v != 0) { // java pads with 0's for us, so there is nothing to
							// do if v==0
				for (int i = a.length; i < newA.length; i++)
					newA[i] = v;
			}
			return newA;
		}

		// we are here only if n < 0:
		n = (-n);
		newA = new short[a.length + n];
		int i;
		if (v != 0) {
			i = 0;
			for (; i < n; i++)
				newA[i] = v;
		} else {
			i = n;
		}
		for (int j = 0; j < a.length; i++, j++)
			newA[i] = a[j];
		return newA;
	}

	/**
	 * a helper method. Turns a single character string into a char.
	 *
	 * @param str
	 *            the string
	 *
	 * @return a char
	 */
	public static char stringToChar(String str) {
		if (str.length() != 1)
			throw new IllegalArgumentException("String length must be one");
		return str.charAt(0);
	}

	public static <T extends Comparable<T>> List<T> sorted(Collection<T> c) {
		return sorted(c, false);
	}

	public static <T extends Comparable<T>> List<T> sorted(Collection<T> c, boolean reverse) {
		List<T> l = new ArrayList<T>(c);
		Collections.sort(l);
		if (reverse)
			Collections.reverse(l);
		return l;
	}

	public static <T extends Comparable<T>, V> List<V> sorted(Map<T, V> c) {
		return sorted(c, false);
	}

	public static <T extends Comparable<T>, V> List<V> sorted(Map<T, V> c, boolean reverse) {
		List<T> t = new ArrayList<T>(c.keySet());
		Collections.sort(t);
		if (reverse)
			Collections.reverse(t);

		List<V> l = new ArrayList<V>();
		for (T k : t) {
			l.add(c.get(k));
		}
		return l;
	}

	public static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {
		List<T> t = new ArrayList<T>(c.keySet());
		Collections.sort(t);

		List<String> pairs = new ArrayList<String>();
		for (T k : t) {
			pairs.add(k + "=" + c.get(k));
		}

		return "{" + join(", ", pairs) + "}";
	}

	/**
	 * Reverse a byte array of bases
	 *
	 * @param bases
	 *            the byte array of bases
	 * @return the reverse of the base byte array
	 */
	static public byte[] reverse(byte[] bases) {
		byte[] rcbases = new byte[bases.length];

		for (int i = 0; i < bases.length; i++) {
			rcbases[i] = bases[bases.length - i - 1];
		}

		return rcbases;
	}

	static public final <T> List<T> reverse(final List<T> l) {
		final List<T> newL = new ArrayList<T>(l);
		Collections.reverse(newL);
		return newL;
	}

	/**
	 * Reverse an int array of bases
	 *
	 * @param bases
	 *            the int array of bases
	 * @return the reverse of the base int array
	 */
	static public int[] reverse(int[] bases) {
		int[] rcbases = new int[bases.length];

		for (int i = 0; i < bases.length; i++) {
			rcbases[i] = bases[bases.length - i - 1];
		}

		return rcbases;
	}

	/**
	 * Reverse (NOT reverse-complement!!) a string
	 *
	 * @param bases
	 *            input string
	 * @return the reversed string
	 */
	static public String reverse(String bases) {
		return new String(reverse(bases.getBytes()));
	}

	public static byte[] charSeq2byteSeq(char[] seqIn) {
		byte[] seqOut = new byte[seqIn.length];
		for (int i = 0; i < seqIn.length; i++) {
			seqOut[i] = (byte) seqIn[i];
		}
		return seqOut;
	}

	public static boolean isFlagSet(int value, int flag) {
		return ((value & flag) == flag);
	}

	/**
	 * Helper utility that calls into the InetAddress system to resolve the
	 * hostname. If this fails, unresolvable gets returned instead.
	 *
	 * @return
	 */
	public static final String resolveHostname() {
		try {
			return InetAddress.getLocalHost().getCanonicalHostName();
		} catch (java.net.UnknownHostException uhe) { // [beware typo in code
														// sample -dmw]
			return "unresolvable";
			// handle exception
		}
	}

	public static byte[] arrayFromArrayWithLength(byte[] array, int length) {
		byte[] output = new byte[length];
		for (int j = 0; j < length; j++)
			output[j] = array[(j % array.length)];
		return output;
	}

	public static void fillArrayWithByte(byte[] array, byte value) {
		for (int i = 0; i < array.length; i++)
			array[i] = value;
	}
	
	public static <E> Collection<E> makeCollection(Iterable<E> iter) {
		Collection<E> list = new ArrayList<E>();
		for (E item : iter) {
			list.add(item);
		}
		return list;
	}

	/**
	 * Returns the number of combinations represented by this collection of
	 * collection of options.
	 *
	 * For example, if this is [[A, B], [C, D], [E, F, G]] returns 2 * 2 * 3 =
	 * 12
	 *
	 * @param options
	 * @param <T>
	 * @return
	 */
	public static <T> int nCombinations(final Collection<T>[] options) {
		int nStates = 1;
		for (Collection<T> states : options) {
			nStates *= states.size();
		}
		return nStates;
	}

	public static <T> int nCombinations(final List<List<T>> options) {
		if (options.isEmpty())
			return 0;
		else {
			int nStates = 1;
			for (Collection<T> states : options) {
				nStates *= states.size();
			}
			return nStates;
		}
	}

	/**
	 * Make all combinations of N size of objects
	 *
	 * if objects = [A, B, C] if N = 1 => [[A], [B], [C]] if N = 2 => [[A, A],
	 * [B, A], [C, A], [A, B], [B, B], [C, B], [A, C], [B, C], [C, C]]
	 *
	 * @param objects
	 * @param n
	 * @param <T>
	 * @param withReplacement
	 *            if false, the resulting permutations will only contain unique
	 *            objects from objects
	 * @return
	 */
	public static <T> List<List<T>> makePermutations(final List<T> objects, final int n,
			final boolean withReplacement) {
		final List<List<T>> combinations = new ArrayList<List<T>>();

		if (n <= 0)
			;
		else if (n == 1) {
			for (final T o : objects)
				combinations.add(Collections.singletonList(o));
		} else {
			final List<List<T>> sub = makePermutations(objects, n - 1, withReplacement);
			for (List<T> subI : sub) {
				for (final T a : objects) {
					if (withReplacement || !subI.contains(a))
						combinations.add(Utils.cons(a, subI));
				}
			}
		}

		return combinations;
	}

	/**
	 * Convenience function that formats the novelty rate as a %.2f string
	 *
	 * @param known
	 *            number of variants from all that are known
	 * @param all
	 *            number of all variants
	 * @return a String novelty rate, or NA if all == 0
	 */
	public static String formattedNoveltyRate(final int known, final int all) {
		return formattedPercent(all - known, all);
	}

	/**
	 * Convenience function that formats the novelty rate as a %.2f string
	 *
	 * @param x
	 *            number of objects part of total that meet some criteria
	 * @param total
	 *            count of all objects, including x
	 * @return a String percent rate, or NA if total == 0
	 */
	public static String formattedPercent(final long x, final long total) {
		return total == 0 ? "NA" : String.format("%.2f", (100.0 * x) / total);
	}

	/**
	 * Convenience function that formats a ratio as a %.2f string
	 *
	 * @param num
	 *            number of observations in the numerator
	 * @param denom
	 *            number of observations in the denumerator
	 * @return a String formatted ratio, or NA if all == 0
	 */
	public static String formattedRatio(final long num, final long denom) {
		return denom == 0 ? "NA" : String.format("%.2f", num / (1.0 * denom));
	}

	/**
	 * Create a constant map that maps each value in values to itself
	 * 
	 * @param values
	 * @param <T>
	 * @return
	 */
	public static <T> Map<T, T> makeIdentityFunctionMap(Collection<T> values) {
		Map<T, T> map = new HashMap<T, T>(values.size());
		for (final T value : values)
			map.put(value, value);
		return Collections.unmodifiableMap(map);
	}

	public static void validateArg(final boolean condition, final Supplier<String> msg) {
		if (!condition) {
			throw new IllegalArgumentException(msg.get());
		}
	}

	public static void validateArg(final boolean condition, final String msg) {
		if (!condition) {
			throw new IllegalArgumentException(msg);
		}
	}

	public static int validIndex(final int index, final int length) {
		if (index < 0) {
			throw new IllegalArgumentException("the index cannot be negative: " + index);
		} else if (index >= length) {
			throw new IllegalArgumentException("the index points past the last element of the collection or array: "
					+ index + " > " + (length - 1));
		}
		return index;
	}

	public static int lastIndexOf(final byte[] reference, final byte[] query) {
		int queryLength = query.length;

		// start search from the last possible matching position and search to
		// the left
		for (int r = reference.length - queryLength; r >= 0; r--) {
			int q = 0;
			while (q < queryLength && reference[r + q] == query[q]) {
				q++;
			}
			if (q == queryLength) {
				return r;
			}
		}
		return -1;
	}

	public static void containsNoNull(final Collection<?> collection, final String message) {
		Utils.nonNull(collection, message);
		// cannot use Collection.contains(null) here because this throws a
		// NullPointerException when used with many Sets
		if (collection.stream().anyMatch(v -> v == null)) {
			throw new IllegalArgumentException(message);
		}
	}

	public static <T> T skimArray(final T source, final int sourceOffset, final T dest, final int destOffset,
			final boolean[] remove, final int removeOffset) {
		Utils.nonNull(source, "the source array cannot be null");

		@SuppressWarnings("unchecked")
		final Class<T> sourceClazz = (Class<T>) source.getClass();

		if (!sourceClazz.isArray()) {
			throw new IllegalArgumentException("the source array is not in fact an array instance");
		}
		final int length = Array.getLength(source) - sourceOffset;
		if (length < 0) {
			throw new IllegalArgumentException("the source offset goes beyond the source array length");
		}
		return skimArray(source, sourceOffset, dest, destOffset, remove, removeOffset, length);
	}

	public static <T> T skimArray(final T source, final int sourceOffset, final T dest, final int destOffset,
			final boolean[] remove, final int removeOffset, final int length) {
		Utils.nonNull(source, "the source array cannot be null");
		Utils.nonNull(remove, "the remove array cannot be null");
		if (sourceOffset < 0) {
			throw new IllegalArgumentException("the source array offset cannot be negative");
		}
		if (destOffset < 0) {
			throw new IllegalArgumentException("the destination array offset cannot be negative");
		}
		if (removeOffset < 0) {
			throw new IllegalArgumentException("the remove array offset cannot be negative");
		}
		if (length < 0) {
			throw new IllegalArgumentException("the length provided cannot be negative");
		}

		final int removeLength = Math.min(remove.length - removeOffset, length);

		if (removeLength < 0) {
			throw new IllegalArgumentException("the remove offset provided falls beyond the remove array end");
		}

		@SuppressWarnings("unchecked")
		final Class<T> sourceClazz = (Class<T>) source.getClass();

		if (!sourceClazz.isArray()) {
			throw new IllegalArgumentException("the source array is not in fact an array instance");
		}

		final Class<T> destClazz = skimArrayDetermineDestArrayClass(dest, sourceClazz);

		final int sourceLength = Array.getLength(source);

		if (sourceLength < length + sourceOffset) {
			throw new IllegalArgumentException("the source array is too small considering length and offset");
		}

		// count how many positions are to be removed.

		int removeCount = 0;

		final int removeEnd = removeLength + removeOffset;
		for (int i = removeOffset; i < removeEnd; i++) {
			if (remove[i]) {
				removeCount++;
			}
		}

		final int newLength = length - removeCount;

		@SuppressWarnings("unchecked")
		final T result = skimArrayBuildResultArray(dest, destOffset, destClazz, newLength);
		// No removals, just copy the whole thing.

		if (removeCount == 0) {
			System.arraycopy(source, sourceOffset, result, destOffset, length);
		} else if (length > 0) { // if length == 0 nothing to do.
			int nextOriginalIndex = 0;
			int nextNewIndex = 0;
			int nextRemoveIndex = removeOffset;
			while (nextOriginalIndex < length && nextNewIndex < newLength) {
				while (nextRemoveIndex < removeEnd && remove[nextRemoveIndex++]) {
					nextOriginalIndex++;
				} // skip positions to be spliced.
				// Since we make the nextNewIndex < newLength check in the while
				// condition
				// there is no need to include the following break, as is
				// guaranteed not to be true:
				// if (nextOriginalIndex >= length) break; // we reach the final
				// (last positions are to be spliced.
				final int copyStart = nextOriginalIndex;
				while (++nextOriginalIndex < length && (nextRemoveIndex >= removeEnd || !remove[nextRemoveIndex])) {
					nextRemoveIndex++;
				}
				final int copyEnd = nextOriginalIndex;
				final int copyLength = copyEnd - copyStart;
				System.arraycopy(source, sourceOffset + copyStart, result, destOffset + nextNewIndex, copyLength);
				nextNewIndex += copyLength;
			}
		}
		return result;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> T skimArrayBuildResultArray(final T dest, final int destOffset, final Class<T> destClazz, final int newLength) {
	    final T result;

	    if (dest == null) {
	        result = (T) Array.newInstance(destClazz.getComponentType(), newLength + destOffset);
	    } else if (Array.getLength(dest) < newLength + destOffset) {
	        result = (T) Array.newInstance(destClazz.getComponentType(),newLength + destOffset);
	        if (destOffset > 0) {
	            System.arraycopy(dest, 0, result, 0, destOffset);
	        }
	    } else {
	        result = dest;
	    }
	    return result;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> Class<T> skimArrayDetermineDestArrayClass(final T dest, final Class<T> sourceClazz) {
	    final Class<T> destClazz;
	    if (dest == null) {
	        destClazz = sourceClazz;
	    } else {
	        destClazz = (Class<T>) dest.getClass();
	        if (destClazz != sourceClazz) {
	            if (!destClazz.isArray()) {
	                throw new IllegalArgumentException("the destination array class must be an array");
	            }
	            if (sourceClazz.getComponentType().isAssignableFrom(destClazz.getComponentType())) {
	                throw new IllegalArgumentException("the provided destination array class cannot contain values from the source due to type incompatibility");
	            }
	        }
	    }
	    return destClazz;
	}
	
	public static boolean equalRange(final byte[] left, final int leftOffset, final byte[] right, final int rightOffset, final int length) {
	    Utils.nonNull(left, "left cannot be null");
	    Utils.nonNull(right, "right cannot be null");

	    for (int i = 0; i < length; i++) {
	        if (left[leftOffset + i] != right[rightOffset + i]) {
	            return false;
	        }
	    }
	    return true;
	}
	
	public static <T> Stream<T> stream(final Iterable<T> iterable) {
	    return StreamSupport.stream(iterable.spliterator(), false);
	}
}
