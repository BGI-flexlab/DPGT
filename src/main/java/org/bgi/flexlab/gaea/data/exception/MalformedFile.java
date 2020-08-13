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

import java.io.File;

public class MalformedFile extends UserException {
	/**
	 * 
	 */
	private static final long serialVersionUID = -3457749925016903826L;

	public MalformedFile(String message) {
        super(String.format("Unknown file is malformed: %s", message));
    }

    public MalformedFile(String message, Exception e) {
        super(String.format("Unknown file is malformed: %s caused by %s", message, getMessage(e)));
    }

    public MalformedFile(File f, String message) {
        super(String.format("File %s is malformed: %s", f.getAbsolutePath(), message));
    }

    public MalformedFile(File f, String message, Exception e) {
        super(String.format("File %s is malformed: %s caused by %s", f.getAbsolutePath(), message, getMessage(e)));
    }

    public MalformedFile(String name, String message) {
        super(String.format("File associated with name %s is malformed: %s", name, message));
    }

    public MalformedFile(String name, String message, Exception e) {
        super(String.format("File associated with name %s is malformed: %s caused by %s", name, message, getMessage(e)));
    }
}
