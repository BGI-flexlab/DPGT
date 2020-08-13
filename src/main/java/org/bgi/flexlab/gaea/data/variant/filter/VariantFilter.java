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
package org.bgi.flexlab.gaea.data.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;

import java.io.IOException;
import java.util.ArrayList;

public abstract class VariantFilter {
	public abstract int filter(VariantContext context, int start, int end);

	public ArrayList<VariantContext> listFilter(ArrayList<VariantContext> contexts, int start, int end) {
		ArrayList<VariantContext> filter = new ArrayList<VariantContext>();

		if (contexts != null) {
			for (VariantContext context : contexts) {
				int result = filter(context, start, end);
				if (result == 0)
					continue;
				if (result == 1)
					filter.add(context);
				else
					break;
			}
		}
		return filter;
	}

	public ArrayList<VariantContext> loadFilter(VCFLocalLoader loader, String referenceName, long startPosition,
			int end) {
		if(startPosition < 0)
			throw new RuntimeException(String.format("start position  %d is less than zero",startPosition));
		try {
			return loader.query(referenceName, startPosition, end);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
}
