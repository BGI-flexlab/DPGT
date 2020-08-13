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

import org.bgi.flexlab.gaea.data.exception.UserException;

public enum EventType {
    SNP(0, "M"),
    Insertion(1, "I"),
    Deletion(2, "D");

    public final int index;
    private final String representation;

    private EventType(int index, String representation) {
        this.index = index;
        this.representation = representation;
    }

    public static EventType eventFrom(int index) {
        switch (index) {
            case 0:
                return SNP;
            case 1:
                return Insertion;
            case 2:
                return Deletion;
            default:
                throw new UserException(String.format("Event %d does not exist.", index));
        }        
    }
    
    public static EventType eventFrom(String event) {
        for (EventType eventType : EventType.values())
            if (eventType.representation.equals(event))
                return eventType;
        throw new UserException(String.format("Event %s does not exist.", event));
    }

    @Override
    public String toString() {
        return representation;
    }
}
