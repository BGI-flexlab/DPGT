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

    public static DPGTJobState readStateFile(String inputStateFile) {
        DPGTJobState s = null;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(inputStateFile));
            StringBuilder strBuilder = new StringBuilder();
            String line = null;
            while ((line = reader.readLine()) != null) {
                strBuilder.append(line);
            }
            s = JSON.parseObject(strBuilder.toString(), DPGTJobState.class);
            reader.close();
        } catch (IOException e) {
            s = new DPGTJobState();  // construct a default job state object, state is not run
        }
        return s;
    }

    public abstract R load();

    public boolean isSuccess() {
        this.jobState = readStateFile(this.stateFile);
        if (this.jobState.isSuccess()) {
            return true;
        } else {
            return false;
        }
    }
}
