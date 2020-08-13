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

public class FileNotExistException extends UserException{

	private static final long serialVersionUID = -8379153613046946464L;

	public FileNotExistException(String path) {
		super(String.format("file %s is not exist!", path));
	}
	
	public static class MissingHeaderException extends UserException{

		/**
		 * 
		 */
		private static final long serialVersionUID = -1223885852313830694L;

		public MissingHeaderException(String path) {
			super(String.format("file %s is not exist header!", path));
		}
		
	}
}
