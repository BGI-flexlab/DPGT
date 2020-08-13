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
package org.bgi.flexlab.gaea.tools.recalibrator;

import htsjdk.samtools.SAMReadGroupRecord;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.sequenceplatform.NGSPlatform;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.QualityCovariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.ReadGroupCovariate;
import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.ReadUtils;

public class RecalibratorUtil {
	/**
	 * table name
	 */
	public final static String ARGUMENT_TABLE_NAME = "Argument";
	public final static String QUANTIZED_TABLE_NAME = "Quantized";
	public final static String[] RECALIBRATOR_TABLE_NAME = { "ReadGroupTable", "QualityScoreTable", "CovariateTable" };

	/**
	 * table column name
	 */
	public final static String ARGUMENT_VALUE_COLUMN_NAME = "Value";
	public final static String QUANTIZED_SCORE_COLUMN_NAME = QualityCovariate.class.getSimpleName()
			.split("Covariate")[0];
	public final static String QUANTIZED_VALUE_COLUMN_NAME = "QuantizedScore";
	public static final String QUANTIZED_COUNT_COLUMN_NAME = "Count";
	public final static String COVARIATE_VALUE_COLUMN_NAME = "CovariateValue";
	public final static String COVARIATE_NAME_COLUMN_NAME = "CovariateName";
	public final static String NUMBER_OBSERVATIONS_COLUMN_NAME = "Observations";
	public final static String NUMBER_ERRORS_COLUMN_NAME = "Errors";
	public final static String EVENT_TYPE_COLUMN_NAME = "EventType";
	public final static String EMPIRICAL_QUALITY_COLUMN_NAME = "EmpiricalQuality";
	public final static String ESTIMATED_Q_REPORTED_COLUMN_NAME = "EstimatedQReported";
	public final static String READGROUP_COLUMN_NAME = ReadGroupCovariate.class.getSimpleName().split("Covariate")[0];

	public enum SolidNocallStrategy {
		THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, PURGE_READ;

		public static SolidNocallStrategy nocallStrategyFromString(String nocallStrategy) {
			if (nocallStrategy.equals("THROW_EXCEPTION"))
				return SolidNocallStrategy.THROW_EXCEPTION;
			if (nocallStrategy.equals("LEAVE_READ_UNRECALIBRATED"))
				return SolidNocallStrategy.LEAVE_READ_UNRECALIBRATED;
			if (nocallStrategy.equals("PURGE_READ"))
				return SolidNocallStrategy.PURGE_READ;

			throw new UserException.BadArgumentValueException(nocallStrategy,
					"is not a valid solid nocall strategy value");
		}
	}

	public enum SolidRecallMode {

		DO_NOTHING,
		/**
		 * Set reference inserted bases and the previous base .This is the
		 * default option.
		 */
		SET_Q_ZERO,
		/**
		 * In addition to setting the quality scores to zero, also set the base
		 * itself to 'N'. This is useful to visualize in IGV.
		 */
		SET_Q_ZERO_BASE_N,
		/**
		 * Look at the color quality scores and probabilistically decide to
		 * change the reference inserted base to be the base which is implied by
		 * the original color space instead of the reference.
		 */
		REMOVE_REF_BIAS;

		public static SolidRecallMode recalModeFromString(String recalMode) {
			if (recalMode.equals("DO_NOTHING"))
				return SolidRecallMode.DO_NOTHING;
			if (recalMode.equals("SET_Q_ZERO"))
				return SolidRecallMode.SET_Q_ZERO;
			if (recalMode.equals("SET_Q_ZERO_BASE_N"))
				return SolidRecallMode.SET_Q_ZERO_BASE_N;
			if (recalMode.equals("REMOVE_REF_BIAS"))
				return SolidRecallMode.REMOVE_REF_BIAS;

			throw new UserException.BadArgumentValueException(recalMode, "is not a valid SOLID_RECAL_MODE value");
		}
	}

	private final static String COLOR_SPACE_INCONSISTENCY_TAG = "ZC";
	private final static String COLOR_SPACE_ATTRIBUTE_TAG = "CS";
	private final static byte MINIMUM_COLOR = '0';
	private final static byte MAXIMUM_COLOR = '3';

	private static boolean hasNoCallInColorSpace(final byte[] colorSpace) {
		final int length = colorSpace.length;
		for (int i = 1; i < length; i++) {
			final byte color = colorSpace[i];
			if (color < MINIMUM_COLOR || color > MAXIMUM_COLOR) {
				return true;
			}
		}

		return false;
	}

	private static byte[] colorIndexMap = { 'C', 'A', 'G', 'T', 'G', 'T', 'C', 'A', 'T', 'G', 'A', 'C' };

	private static byte performColor(byte base, int idx) {
		int index = BaseUtils.simpleBaseToBaseIndex(base);

		if (index == -1)
			return base;

		int baseIndex = (idx - MINIMUM_COLOR - 1) * 4 + index;
		return colorIndexMap[baseIndex];
	}

	private static byte getNextBaseFromColor(GaeaSamRecord read, final byte prevBase, final byte color) {
		if (color < MINIMUM_COLOR || color > MAXIMUM_COLOR)
			throw new UserException.MalformedBAM(read,
					"Unrecognized color space in SOLID read, color = " + (char) color);

		if (color == '0')
			return prevBase;

		return performColor(prevBase, color);
	}

	public static void defaultPlatformForRead(GaeaSamRecord read, String forcePlatform, String defaultPlatform) {
		SAMReadGroupRecord readGroup = read.getReadGroup();
		String platform = readGroup.getPlatform();
		if (forcePlatform != null && (platform == null || !platform.equals(forcePlatform)))
			readGroup.setPlatform(forcePlatform);

		if (readGroup.getPlatform() == null) {
			if (defaultPlatform != null)
				readGroup.setPlatform(defaultPlatform);
		}
	}

	public static class Consistent {
		private boolean isConsistent;
		private byte[] inconsistency = null;

		public Consistent(boolean isConsistent, byte[] inconsistency) {
			this.isConsistent = isConsistent;
			this.inconsistency = inconsistency;
		}

		public boolean isColorSpaceConsistent() {
			return this.isConsistent;
		}

		public boolean isColorSpaceConsistent(int offset, boolean negativeStrand) {
			if (inconsistency == null)
				return true;
			if (negativeStrand)
				return inconsistency[inconsistency.length - offset - 1] == (byte) 0;
			return inconsistency[offset] == (byte) 0;
		}
	}

	public static Consistent isColorSpaceConsistent(final SolidNocallStrategy strategy, final GaeaSamRecord read) {
		if (!ReadUtils.isSOLiDRead(read))
			return new Consistent(true, null);

		if (read.getAttribute(COLOR_SPACE_INCONSISTENCY_TAG) == null) {
			final Object attr = read.getAttribute(COLOR_SPACE_ATTRIBUTE_TAG);
			if (attr != null) {
				byte[] colorSpace;
				if (attr instanceof String)
					colorSpace = ((String) attr).getBytes();
				else
					throw new UserException.MalformedBAM(read,
							String.format("Value encoded by %s in %s isn't a string!", COLOR_SPACE_ATTRIBUTE_TAG,
									read.getReadName()));

				final boolean badColor = hasNoCallInColorSpace(colorSpace);
				if (badColor) {
					if (strategy == SolidNocallStrategy.LEAVE_READ_UNRECALIBRATED) {
						return new Consistent(false, null);
					} else if (strategy == SolidNocallStrategy.PURGE_READ) {
						read.setReadFailsVendorQualityCheckFlag(true);
						return new Consistent(false, null);
					}
				}

				byte[] readBases = read.getReadBases();
				if (read.getReadNegativeStrandFlag())
					readBases = BaseUtils.simpleReverseComplement(read.getReadBases());

				final byte[] inconsistency = new byte[readBases.length];
				int i;
				byte prevBase = colorSpace[0];
				for (i = 0; i < readBases.length; i++) {
					final byte thisBase = getNextBaseFromColor(read, prevBase, colorSpace[i + 1]);
					inconsistency[i] = (byte) (thisBase == readBases[i] ? 0 : 1);
					prevBase = readBases[i];
				}
				return new Consistent(true, inconsistency);
			} else if (strategy == SolidNocallStrategy.THROW_EXCEPTION) {
				throw new UserException.MalformedBAM(read,
						"Unable to find color space information in SOLiD read. First observed at read with name = "
								+ read.getReadName());
			} else {
				return new Consistent(false, null);
			}
		}

		return new Consistent(true, null);
	}

	public static ReadCovariates computeCovariates(final GaeaSamRecord read, final Covariate[] covariates) {
		final ReadCovariates readCovariates = new ReadCovariates(read.getReadLength(), covariates.length);
		computeCovariates(read, covariates, readCovariates);
		return readCovariates;
	}

	public static void computeCovariates(final GaeaSamRecord read, final Covariate[] covariates,
			ReadCovariates readCovariates) {
		for (int i = 0; i < covariates.length; i++) {
			readCovariates.setCovariateIndex(i);
			covariates[i].recordValues(read, readCovariates);
		}
	}

	public static boolean isSOLiDRead(GaeaSamRecord read) {
		return NGSPlatform.fromRead(read) == NGSPlatform.SOLID;
	}
}
