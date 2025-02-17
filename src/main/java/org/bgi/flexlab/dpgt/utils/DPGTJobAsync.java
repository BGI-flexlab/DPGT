/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
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
