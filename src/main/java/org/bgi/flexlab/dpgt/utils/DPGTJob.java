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
        jobState = readStateFile(stateFile);
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
