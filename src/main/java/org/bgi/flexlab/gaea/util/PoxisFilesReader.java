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
package org.bgi.flexlab.gaea.util;

import java.io.*;
import java.util.ArrayList;

public class PoxisFilesReader extends GaeaFilesReader {
	private ArrayList<File> files = new ArrayList<File>();
	private BufferedReader bufferedReader = null;

	@Override
	public void traversal(String path) {
		File file = new File(path);

		if (file.isDirectory()) {
			File[] fileArray = file.listFiles();
			for (File f : fileArray) {
				if (f.isDirectory())
					traversal(file.getAbsolutePath());
				else {
					if (!filter(f.getName()))
						files.add(f);
				}
			}
		} else {
			if (!filter(file.getName()))
				files.add(file);
		}
		
		if(size() == 0)
			return;
		
		try {
			bufferedReader = new BufferedReader(new FileReader(files.get(0)));
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e.toString());
		}
	}

	protected int size() {
		return files.size();
	}

	@Override
	public boolean hasNext() {
		if (bufferedReader != null) {
			String str;
			try {
				if ((str = bufferedReader.readLine()) != null) {
					currentLine = str;
					return true;
				} else {
					currentFileIndex++;
					bufferedReader.close();
					if (currentFileIndex < size()) {
						bufferedReader = new BufferedReader(new FileReader(files.get(currentFileIndex)));
						if ((str = bufferedReader.readLine()) != null) {
							currentLine = str;
							return true;
						}
					}
				}
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
		}

		currentLine = null;
		return false;
	}

	@Override
	public void clear() {
		if(bufferedReader != null){
			try {
				bufferedReader.close();
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
		}
		files.clear();
		currentLine = null;
		currentFileIndex = 0;
	}
}
