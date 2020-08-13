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

public class BAMQCException extends RuntimeException{
	/**
	 * 
	 */
	private static final long serialVersionUID = -7654073424139393545L;


	public BAMQCException(String msg) { super(msg); }
    public BAMQCException(String msg, Throwable e) { super(msg, e); }
    public BAMQCException(Throwable e) { super("", e); } // cannot be called, private access

	
	public static class WrongNumOfColException extends BAMQCException {
		/**
		 * 
		 */
		private static final long serialVersionUID = -4596452748966703629L;

		public WrongNumOfColException(int columnNum) {
            super(String.format("column number should be: %s", columnNum));
        }
	}
	
	public static class GenderInformationException extends BAMQCException {
		/**
		 * 
		 */
		private static final long serialVersionUID = -6457897353260936041L;

		public GenderInformationException(String message) {
			super(String.format("Gender Information wrong: %s", message));
		}
	}
}
