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

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReportTable;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.Random;

public class RecalibratorDatum {
	private final static double DEFAULT_DOUBLE = -1.0;
	private final static int CONSTANT = 1;

	/**
	 * number of bases in total
	 */
	private long numBases;

	/**
	 * number of bases that didn't match the reference
	 */
	private long numMismatches;

	/**
	 * estimated reported quality score
	 */
	private double estimatedQuality;

	/**
	 * the empirical quality for datums
	 */
	private double empiricalQuality;

	public RecalibratorDatum(long _numBases, long _numMismatches, double _estimatedQuality) {
		if (_numBases < 0)
			throw new IllegalArgumentException("base number < 0");
		if (_numMismatches < 0)
			throw new IllegalArgumentException("mismatch base number < 0");
		if (_estimatedQuality < 0)
			throw new IllegalArgumentException("estimated quality < 0");

		this.numBases = _numBases;
		this.numMismatches = _numMismatches;
		this.estimatedQuality = _estimatedQuality;
		this.empiricalQuality = DEFAULT_DOUBLE;
	}

	public RecalibratorDatum(RecalibratorDatum clone) {
		this.numBases = clone.getBasesNumber();
		this.numMismatches = clone.numMismatches;
		this.estimatedQuality = clone.getEstimatedQuality();
		this.empiricalQuality = clone.getEmpiricalQuality();
	}

	public static RecalibratorDatum build(byte quality, boolean isError) {
		return new RecalibratorDatum(1, isError ? 1 : 0, quality);
	}

	public static RecalibratorDatum build(int maxBasesNumber, int maxErrors) {
		final Random random = new Random();
		final int nObservations = random.nextInt(maxBasesNumber);
		final int nErrors = random.nextInt(maxErrors);
		final int qual = random.nextInt(QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE);
		return new RecalibratorDatum(nObservations, nErrors, (byte) qual);
	}

	public static RecalibratorDatum build(RecalibratorReportTable reportTable, final int row,
			final boolean hasEstimatedColumn) {
		final long nObservations = (Long) reportTable.get(row, RecalibratorUtil.NUMBER_OBSERVATIONS_COLUMN_NAME);
		final long nErrors = (Long) reportTable.get(row, RecalibratorUtil.NUMBER_ERRORS_COLUMN_NAME);
		final double empiricalQuality = (Double) reportTable.get(row, RecalibratorUtil.EMPIRICAL_QUALITY_COLUMN_NAME);

		final double estimatedQReported = hasEstimatedColumn
				? (Double) reportTable.get(row, RecalibratorUtil.ESTIMATED_Q_REPORTED_COLUMN_NAME)
				: Byte.parseByte((String) reportTable.get(row, RecalibratorUtil.QUANTIZED_SCORE_COLUMN_NAME));

		RecalibratorDatum datum = new RecalibratorDatum(nObservations, nErrors, (byte) 1);
		datum.setEstimatedQReported(estimatedQReported);
		datum.setEmpiricalQuality(empiricalQuality);
		return datum;
	}

	public void combine(RecalibratorDatum other) {
		double combineErrors = this.expectedErrors() + other.expectedErrors();
		addBaseNumber(other.getBasesNumber());
		this.estimatedQuality = -10 * Math.log10(combineErrors / this.numBases);
		addMismatchNumber(other.numMismatches);
	}

	public void increment(RecalibratorDatum other) {
		if (other == null)
			throw new UserException("recalibrator datum cann't be null !");

		addBaseNumber(other.getBasesNumber());
		addMismatchNumber(other.getMismatchNumber());
	}

	public void increment(boolean isError) {
		addBaseNumber(1);
		if (isError)
			addMismatchNumber(1);
	}

	public void addBaseNumber(long number) {
		this.numBases += number;
		this.empiricalQuality = DEFAULT_DOUBLE;
	}

	public long getBasesNumber() {
		return this.numBases;
	}

	public void addMismatchNumber(long mismatchNumber) {
		this.numMismatches += mismatchNumber;
		this.empiricalQuality = DEFAULT_DOUBLE;
	}

	public long getMismatchNumber() {
		return this.numMismatches;
	}

	public void setEstimatedQReported(final double estimatedQuality) {
		if (estimatedQuality < 0 || Double.isInfinite(estimatedQuality) || Double.isNaN(estimatedQuality))
			throw new IllegalArgumentException("estimatedQReported is error");

		this.estimatedQuality = estimatedQuality;
	}

	public double getEstimatedQuality() {
		return this.estimatedQuality;
	}

	public double getEmpiricalQuality() {
		if (this.empiricalQuality == DEFAULT_DOUBLE)
			calcEmpiricalQuality();
		return this.empiricalQuality;
	}

	public void setEmpiricalQuality(double empir) {
		this.empiricalQuality = empir;
	}

	private final void calcEmpiricalQuality() {
		final double empiricalQual = -10 * Math.log10(getEmpiricalErrorRate());
		this.empiricalQuality = Math.min(empiricalQual, (double) QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE);
	}

	public double getEmpiricalErrorRate() {
		if (numBases == 0)
			return 0.0;
		else {
			double mismathes = (double) (numMismatches + CONSTANT);
			double bases = (double) (numBases + (CONSTANT << 1));
			return mismathes / bases;
		}
	}

	public double expectedErrors() {
		return this.numBases * QualityUtils.qualityToErrorProbability(this.estimatedQuality);
	}

	public String toString() {
		return String.format("%d\t%d\t%.4f", getBasesNumber(), getMismatchNumber(), estimatedQuality);
	}
}
