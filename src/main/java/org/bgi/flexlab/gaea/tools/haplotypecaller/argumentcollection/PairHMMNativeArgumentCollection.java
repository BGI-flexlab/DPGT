package org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;

public class PairHMMNativeArgumentCollection {

    private int pairHmmNativeThreads = 4;

    private boolean useDoublePrecision = false;

    public PairHMMNativeArguments getPairHMMArgs(){
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        args.maxNumberOfThreads = pairHmmNativeThreads;
        args.useDoublePrecision = useDoublePrecision;
        return args;
    }

}
