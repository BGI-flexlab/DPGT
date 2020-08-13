package org.bgi.flexlab.gaea.data.structure.pileup.filter;

import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;

import java.util.ArrayList;

/**
 * Created by zhangyong on 2017/4/13.
 */
public interface PileupFilter {
    ArrayList<PileupReadInfo> filter(Pileup pileup);
}
