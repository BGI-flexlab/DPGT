/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.vcfstats.report;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfstats.VCFStatsOptions;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class VCFReport {

    private Map<String, PerSampleVCFReport> perSampleVCFReports;
    private VCFStatsOptions options;
    private ReferenceShare genomeShare;
    private VCFFileReader vcfReader;

    public VCFReport(){
        perSampleVCFReports = new HashMap<>();
    }

    public VCFReport(VCFStatsOptions options){
        perSampleVCFReports = new HashMap<>();
        this.options = options;
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReferenceSequencePath());
        if(options.getDbsnpFile() != null)
            vcfReader = new VCFFileReader(new File(options.getDbsnpFile()));
    }

    public void parseVariation(VariantContext vc){

        if(options.getDbsnpFile() != null)
            vc = setDbSNP(vc);

        PerSampleVCFReport aSampleVCFReport;
        for(String sample: vc.getSampleNames()){
            Genotype gt = vc.getGenotype(sample);
            if(!gt.isCalled())
                return;

            if(perSampleVCFReports.containsKey(sample))
                aSampleVCFReport = perSampleVCFReports.get(sample);
            else {
                aSampleVCFReport = new PerSampleVCFReport(genomeShare, options.isCountVarLength());
                perSampleVCFReports.put(sample, aSampleVCFReport);
            }

            aSampleVCFReport.add(vc, sample);
        }
    }

    private VariantContext setDbSNP(VariantContext vc) {
        CloseableIterator<VariantContext> vcfIter = vcfReader.query(vc.getContig(), vc.getStart(), vc.getEnd());

        String id = null;
        while (vcfIter.hasNext()) {
            VariantContext dbsnpvc = vcfIter.next();
            if (dbsnpvc.getStart() != vc.getStart() || dbsnpvc.getEnd() != vc.getEnd())
                continue;

            for(Allele allele: vc.getAlternateAlleles()){
                for(Allele dnsnpAllele: dbsnpvc.getAlternateAlleles()) {
                    if(allele.equals(dnsnpAllele)) {
                        id = dbsnpvc.getID();
                        break;
                    }
                }
                if(id != null)
                    break;
            }
            if(id != null)
                break;
        }

        if(id != null) {
            return new VariantContextBuilder(vc).id(id).make();
        }
        return vc;
    }


    public String parseReducerString(){
        return null;
    }

    public String toReducerString(){
        return null;
    }

    public void readFromHdfs(Path path, Configuration conf) throws IOException {
        FileSystem fs = path.getFileSystem(conf);
        FSDataInputStream FSinput = fs.open(path);

        LineReader lineReader = new LineReader(FSinput, conf);
        Text line = new Text();
        PerSampleVCFReport perSampleVCFReport;
        PerSampleVCFReport totalVCFReport;
        if(!perSampleVCFReports.containsKey("total")) {
            totalVCFReport = new PerSampleVCFReport(genomeShare, options.isCountVarLength());
            perSampleVCFReports.put("total", totalVCFReport);
        }else {
            totalVCFReport = perSampleVCFReports.get("total");
        }
        while ((lineReader.readLine(line)) != 0) {
            String reducerStr = line.toString();
            if(reducerStr.isEmpty()) continue;
            String[] fields = reducerStr.split("\t",2);
            String sample = fields[0];
            if(perSampleVCFReports.containsKey(sample))
                perSampleVCFReport = perSampleVCFReports.get(sample);
            else {
                perSampleVCFReport = new PerSampleVCFReport(genomeShare, options.isCountVarLength());
                perSampleVCFReport.setSampleName(sample);
                perSampleVCFReports.put(sample, perSampleVCFReport);
            }
            perSampleVCFReport.parseReducerString(fields[1]);
            totalVCFReport.parseReducerString(fields[1]);
        }
        lineReader.close();
    }

    private void write(String fileName, Configuration conf, String report)
            throws IOException {
        Path reportFilePath = new Path(fileName);
        FileSystem fs = reportFilePath.getFileSystem(conf);
        FSDataOutputStream writer = fs.create(reportFilePath);
        writer.write(report.getBytes());
        writer.close();
    }

    public void mergeReport(Path input, Configuration conf, Path outputDir) throws IOException {
        FileSystem fs = input.getFileSystem(conf);
        FileStatus filelist[] = fs.listStatus(input);
//        FileStatus filelist[] = fs.listStatus(input,new StaticPathFilter());

        for (FileStatus aFilelist : filelist) {
            if (!aFilelist.isDirectory()) {
                readFromHdfs(aFilelist.getPath(), conf);
                fs.delete(aFilelist.getPath(), false);
            }
        }

        fs.close();
        for(String sample: getPerSampleVCFReports().keySet()){
            PerSampleVCFReport PerSampleVCFReport = getPerSampleVCFReports().get(sample);
            String fileName = outputDir + "/" + sample + ".vcfstats.report.txt";
            write(fileName, conf, PerSampleVCFReport.getReport());
        }
    }

    public Map<String, PerSampleVCFReport> getPerSampleVCFReports() {
        return perSampleVCFReports;
    }

    static class StaticPathFilter implements PathFilter {
        @Override
        public boolean accept(Path path) {
            if (path.getName().startsWith("Statistic"))
                return true;
            return false;
        }
    }
}
