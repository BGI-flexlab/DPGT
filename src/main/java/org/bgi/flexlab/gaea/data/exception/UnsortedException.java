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
package org.bgi.flexlab.gaea.data.exception;

public class UnsortedException extends UserException {
	/**
	 * 
	 */
	private static final long serialVersionUID = 9013680156879378293L;

	public UnsortedException(String readName, int lastestStart, int start) {
		super(
				String.format(
						"incoming objects must ordered by alignment start,but saw the reads %s with alignment start %d after lastest start %d",
						readName, start, lastestStart));
	}
}
