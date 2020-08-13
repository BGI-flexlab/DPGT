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
package org.bgi.flexlab.gaea.tools.annotator.db;

import org.bgi.flexlab.gaea.tools.annotator.AnnotationContext;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;

public class GeneDBQuery extends DBQuery {

	@Override
	public Results query(Condition condition)throws IOException{
		Results results = new Results();

		for (String gene : condition.getGenes()) {
			HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getTable(), gene, condition.getFields());
			if (result ==null || result.isEmpty()) return null;
			results.add(gene, result);
		}
		return results;
	}

	@Override
	public HashMap<String,String> getMergeResult(AnnotationContext ac) {
		return results.getMergeResult(ac.getGeneName());
	}

}
