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

import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class FragmentCollection<T> {
	Collection<T> singletons;
    Collection<List<T>> overlappingPairs;

    public FragmentCollection(final Collection<T> singletons, final Collection<List<T>> overlappingPairs) {
        this.singletons = singletons == null ? Collections.<T>emptyList() : singletons;
        this.overlappingPairs = overlappingPairs == null ? Collections.<List<T>>emptyList() : overlappingPairs;
    }

    /**
     * Gets the T elements not containing overlapping elements, in no particular order
     */
    public Collection<T> getSingletonReads() {
        return singletons;
    }

    /**
     * Gets the T elements containing overlapping elements, in no particular order
     */
    public Collection<List<T>> getOverlappingPairs() {
        return overlappingPairs;
    }
}
