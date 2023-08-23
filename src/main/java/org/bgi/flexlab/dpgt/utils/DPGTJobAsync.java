package org.bgi.flexlab.dpgt.utils;

import java.util.List;
import java.util.concurrent.TimeUnit;
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

    /**
     * Waits if necessary for at most the given time for the computation
     * to complete, and then retrieves its result, if available.
     * 
     * @param timeout the maximum time to wait
     * @param unit the time unit of the timeout argument
     * @return the compute result
     */
    public abstract R get(long timeout, TimeUnit unit);

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
