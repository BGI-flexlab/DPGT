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
package org.bgi.flexlab.gaea.tools.recalibrator.report;

import org.apache.commons.lang.math.NumberUtils;
import org.bgi.flexlab.gaea.util.StandardDataType;

import java.util.Arrays;
import java.util.Collection;

public class RecalibratorReportColumn {
	private boolean isRightAlignment = false;
	private final StandardDataType dataType;
	private final String format;
	private int maxWidth;
	private String name;

	public RecalibratorReportColumn(final String columnName, final String format) {
		this.name = columnName;
		this.maxWidth = columnName.length();
		if (format.equals("")) {
			this.format = "%s";
			this.dataType = StandardDataType.Unknown;
		} else {
			this.format = format;
			this.dataType = StandardDataType.fromFormatString(format);
		}
	}

	public StandardDataType getDataType() {
		return dataType;
	}

	public void updateFormatting(final Object value) {
		if (value != null) {
			final String formatted = formatValue(value);
			if (formatted.length() > 0) {
				updateMaxWidth(formatted);
				updateFormat(formatted);
			}
		}
	}

	private String formatValue(final Object obj) {
		String value;
		if (obj == null) {
			value = "null";
		} else if (dataType.equals(StandardDataType.Unknown) && (obj instanceof Double || obj instanceof Float)) {
			value = String.format("%.8f", obj);
		} else
			value = String.format(format, obj);

		return value;
	}

	private void updateMaxWidth(final String formatted) {
		maxWidth = Math.max(formatted.length(), maxWidth);
	}

	private static final Collection<String> RIGHT_ALIGN_SETS = Arrays.asList("null", "NA",
			String.valueOf(Double.POSITIVE_INFINITY), String.valueOf(Double.NEGATIVE_INFINITY),
			String.valueOf(Double.NaN));

	protected static boolean isRightAlign(final String value) {
		return value == null || RIGHT_ALIGN_SETS.contains(value) || NumberUtils.isNumber(value.trim());
	}

	private void updateFormat(final String formatted) {
		if (!this.isRightAlignment)
			isRightAlignment = isRightAlign(formatted) ? true : false;
	}

	public String getColumnName() {
		return this.name;
	}

	public String getFormat() {
		return dataType.equals(StandardDataType.Unknown) ? "%s" : format;
	}

	public String getColumnNameFormat() {
		return "%-" + maxWidth + "s";
	}

	public String getColumnValueFormat() {
		if (isRightAlignment)
			return "%" + maxWidth + "s";
		else
			return "%-" + maxWidth + "s";
	}
}
