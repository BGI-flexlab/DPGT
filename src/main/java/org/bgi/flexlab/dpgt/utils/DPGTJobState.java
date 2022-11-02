package org.bgi.flexlab.dpgt.utils;

import java.util.List;
import java.util.ArrayList;
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
