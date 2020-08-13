package org.bgi.flexlab.gaea.data.structure.pileup.filter;

import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;

import java.util.ArrayList;

/**
 * Created by zhangyong on 2017/4/13.
 */
public class PileupMappingBaseQualFilter implements PileupFilter{
    private int minMappingQuality = -1;

    private int minBaseQuality = -1;

    public PileupMappingBaseQualFilter(int minMappingQuality, int minBaseQuality) {
        this.minMappingQuality = minMappingQuality;
        this.minBaseQuality = minBaseQuality;
    }

    @Override
    public ArrayList<PileupReadInfo> filter(Pileup pileup) {
        if(pileup.getTotalPileup().size() == 0)
            return null;

        ArrayList<PileupReadInfo> filterPileup = new ArrayList<>();
        for(PileupReadInfo readInfo : pileup.getTotalPileup()) {
            if(readInfo.getMappingQuality() >= minMappingQuality && (readInfo.isDeletionBase() || readInfo.getBaseQuality() >= minBaseQuality)) {
                filterPileup.add(readInfo);
            }
        }
        return filterPileup;
    }
}
