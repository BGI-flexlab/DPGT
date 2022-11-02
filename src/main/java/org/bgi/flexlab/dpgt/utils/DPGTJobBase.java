package org.bgi.flexlab.dpgt.utils;

import java.io.Serializable;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import com.alibaba.fastjson.JSON;

/**
 * a job that checks job state(success, not run, failed) before runing the job
 */
public abstract class DPGTJobBase<R> implements Serializable {
    private static final Logger logger = LoggerFactory.getLogger(DPGTJobBase.class);

    protected String stateFile = null;
    protected DPGTJobState jobState = null;

    public void writeStateFile() {
        final String jsonStr =  JSON.toJSONString(jobState);
        try {
            BufferedOutputStream outputStream = new BufferedOutputStream(new FileOutputStream(stateFile));
            outputStream.write(jsonStr.getBytes());
            outputStream.close();
        } catch (IOException e) {
            logger.error("{}", e.getMessage());
            System.exit(1);
        }
    }

    public void readStateFile() {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(stateFile));
            StringBuilder strBuilder = new StringBuilder();
            String line = null;
            while ((line = reader.readLine()) != null) {
                strBuilder.append(line);
            }
            jobState = JSON.parseObject(strBuilder.toString(), DPGTJobState.class);
            reader.close();
        } catch (IOException e) {
            jobState = new DPGTJobState();  // construct a default job state object, state is not run
        }
    }

    public abstract R load();

    public boolean isSuccess() {
        readStateFile();
        if (this.jobState.isSuccess()) {
            return true;
        } else {
            return false;
        }
    }
}
