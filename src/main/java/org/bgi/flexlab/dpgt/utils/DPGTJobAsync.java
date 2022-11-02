package org.bgi.flexlab.dpgt.utils;

import java.util.List;

import org.apache.spark.api.java.JavaFutureAction;

public abstract class DPGTJobAsync<T, R> extends DPGTJobBase<R>  {
    // private static final Logger logger = LoggerFactory.getLogger(DPGTJob.class);

    protected List<JavaFutureAction<T>> futures = null;

    /**
     * submit the job, non-blocking
     * @return
     */
    public abstract void submit();

    /**
     * get result from futures and set jobstate, blocking
     * @return
     */
    public abstract R get();

    // public R run() {
    //     readStateFile();
    //     R result = null;
    //     if (!jobState.isSuccess()) {
    //         // job is not success, work on it
    //         result = work();
    //         jobState.jobState = DPGTJobState.State.SUCCESS;
    //         writeStateFile();
    //     } else {
    //         // job is success, load result
    //         result = load();
    //     }
    //     return result;
    // }

}
