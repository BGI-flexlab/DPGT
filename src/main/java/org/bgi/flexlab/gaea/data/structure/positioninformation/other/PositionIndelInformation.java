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
package org.bgi.flexlab.gaea.data.structure.positioninformation.other;

import org.bgi.flexlab.gaea.data.structure.bam.SAMInformationBasic;
import org.bgi.flexlab.gaea.data.structure.positioninformation.BooleanPositionInformation;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CalculatePositionInforamtionInterface;
import org.bgi.flexlab.gaea.data.structure.positioninformation.CompoundInformation;
import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.SystemConfiguration;



public class PositionIndelInformation extends BooleanPositionInformation implements CalculatePositionInforamtionInterface<SAMInformationBasic>{

	public PositionIndelInformation(int windowSize) {
		super(windowSize);
		// TODO Auto-generated constructor stub
	}

	@SuppressWarnings("rawtypes")
	@Override
	public void add(CompoundInformation posInfo) {
		// TODO Auto-generated method stub
		CigarState cigarState = posInfo.getCigarState();
		int cigar = cigarState.getCurrentCigar();
		int cigarLength = (cigar >> 4);
		
		if (eligibleCigar(cigarState, cigar, cigarLength, posInfo)) {
			cigar = cigarState.getCigar().get(cigarState.getCigarState()[0] + 1);
			int cigarOP = (cigar & 0xf);
			if(cigarOP == SystemConfiguration.BAM_CDEL || cigarOP == SystemConfiguration.BAM_CINS) {
				info[(int) (posInfo.distBetweenRefPosAndWinStart())] = true;
			}
		}
	}
	
	@SuppressWarnings("rawtypes")
	private boolean eligibleCigar(CigarState cigarState, int cigar, int cigarLength, CompoundInformation posInfo) {
		return cigarState.getCigarState()[1] + cigarLength - 1 == posInfo.getRefPosition() && 
				cigarState.getCigarState()[0] + 1 < cigarState.getCigar().size();
	}

}
