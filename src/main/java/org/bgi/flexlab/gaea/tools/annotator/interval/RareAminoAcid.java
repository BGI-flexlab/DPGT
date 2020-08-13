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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.interval;

import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;


/**
 * Rare amino acid annotation:
 * 
 * These are amino acids that occurs very rarely in an organism. For instance, humans 
 * are supposed to use 20 amino acids, but there is also one rare AA. Selenocysteine, 
 * single letter code 'U', appears roughly 100 times in the whole genome. The amino 
 * acid is so rare that usually it does not appear in codon translation tables. It 
 * is encoded as UGA, which , under normal conditions, is a STOP codon. Secondary 
 * RNA structures are assumed to enable this special translation.
 * 
 * @author pcingola
 */
public class RareAminoAcid extends Marker {

	private static final long serialVersionUID = -1926572865764543849L;

	public RareAminoAcid() {
		super();
		type = EffectType.RARE_AMINO_ACID;
	}

	public RareAminoAcid(Marker parent, int start, int end, String id) {
		super(parent, start, end, false, id);
		type = EffectType.RARE_AMINO_ACID;
	}
}
