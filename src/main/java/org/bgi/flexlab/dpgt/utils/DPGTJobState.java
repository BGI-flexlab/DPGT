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
import java.io.File;
import java.util.TreeMap;
import com.alibaba.fastjson.annotation.JSONField;


public class DPGTJobState {

    public enum State {
        NOT_RUN,
        SUCCESS,
        FAIL;
    }

    @JSONField(name = "job_state")
    public State jobState = State.NOT_RUN;

    @JSONField(name = "output")
    public TreeMap<Integer, List<String>> outPutFiles = new TreeMap<>();

    @JSONField(name = "metadata")
    public TreeMap<String, String> metaData = new TreeMap<>();

    /**
     * if
     * 1. jobState is SUCCESS
     * 2. outPutFiles are exits
     * then return true, otherwise return false
     * @return
     */
    @JSONField(serialize = false)
    public boolean isSuccess() {
        if (jobState != State.SUCCESS) return false;
        for (final List<String> files: outPutFiles.values()) {
            for (final String f: files) {
                File file = new File(f);
                if (!file.exists()) {
                    return false;
                }
            }
        }
        return true;
    }
}
