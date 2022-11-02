package org.bgi.flexlab.dpgt.utils;


/**
 * a job that checks job state(success, not run, failed) before runing the job
 */
public abstract class DPGTJob<R> extends DPGTJobBase<R> {
    // private static final Logger logger = LoggerFactory.getLogger(DPGTJob.class);

    /**
     * run the job, blocking
     * @return
     */
    public abstract R work();

    public R run() {
        readStateFile();
        R result = null;
        if (!jobState.isSuccess()) {
            // job is not success, work on it
            result = work();
            jobState.jobState = DPGTJobState.State.SUCCESS;
            writeStateFile();
        } else {
            // job is success, load result
            result = load();
        }
        return result;
    }

}
