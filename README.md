# DPGT: A Distributed Population Genetics analysis Tool

DPGT is a distributed population genetics analysis tool which enabled joint calling on millions  of WGS(whole genome sequencing) samples.

## Getting started


### Usage

```
DPGT: Distributed Population Genetics analysis Tools
Version: 1.1.0.0

Options:
 -i,--input <FILE>                        input gvcf list in a file.
 -o,--output <DIR>                        output directory.
 -r,--reference <FILE>                    reference fasta file.
 -l,--target-regions <STRING>             target regions to process.
 -j,--jobs <INT>                          number of spark exculators. [10]
 -n,--num-combine-partitions <INT>        number of partitions for combining gvcfs. [10]
 -w,--window <INT>                        window size for each combine-genotype cycle. [300M]
 -d,--delete <Boolean>                    delete combine and genotype gvcf intermediate results. Possible values: true, false. [true]
 -s,--stand-call-conf <FLOAT>             the minimum phred-scaled confidence threshold at which variants should be called. [30.0]
    --dbsnp <FILE>                        dbsnp vcf file for annotation.
    --use-old-qual-calculator <Boolean>   use the old AF model. Possible values: true, false. [true]
    --heterozygosity <FLOAT>              heterozygosity value used to compute prior likelihoods for any locus. [0.001]
    --indel-heterozygosity <FLOAT>        heterozygosity for indel calling. [1.25E-4]
    --heterozygosity-stdev <FLOAT>        standard deviation of heterozygosity for SNP and indel calling. [0.01]
    --ploidy <INT>                        ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool *
                                          Sample Ploidy). [2]
    --local                               run spark in local mode, useful for debug.
 -h,--help                                print this message and exit.
 -V,--version                             show version.
```


### Run DPGT on a local computer

```sh
export LD_LIBRARY_PATH=/path_to_dpgt_dir/build/lib:${LD_LIBRARY_PATH}
export LD_PRELOAD=/path_to_jemalloc_5.3.0/lib/libjemalloc.so
spark-submit \
    --conf "spark.dynamicAllocation.enabled=false" \
    --conf spark.memory.fraction=0.01 \
    --conf spark.memory.storageFraction=0.01 \
    --master local[32] \
    --driver-memory 140g --executor-cores 1 \
    --conf spark.dynamicAllocation.enabled=false \
    --class org.bgi.flexlab.dpgt.jointcalling.JointCallingSpark \
    dpgt-<version>.jar \
    -i /path_to/vcfs.list \
    -r hg38.fasta \
    -o result -n 256 -j 32 -l chr1:1000000-2000000
```

`--master local[32]` means run spark application using 32-threads.


### Run DPGT on a yarn cluster

```sh
export LD_LIBRARY_PATH=/path_to_dpgt_dir/build/lib:${LD_LIBRARY_PATH}
export LD_PRELOAD=/path_to_jemalloc_5.3.0/lib/libjemalloc.so
spark-submit \
    --conf spark.dynamicAllocation.enabled=false \
    --conf spark.excutorEnv.LD_PRELOAD=${LD_PRELOAD} \
    --conf spark.driver.extraLibraryPath=${LD_LIBRARY_PATH} \
    --conf spark.executor.extraLibraryPath=${LD_LIBRARY_PATH} \
    --conf spark.driver.maxResultSize=10g \
    --conf spark.memory.fraction=0.01 \
    --conf spark.memory.storageFraction=0.01 \
    --master yarn \
    --driver-memory 40g --num-executors 256 \
    --executor-memory 10g --executor-cores 1 \
    --class org.bgi.flexlab.dpgt.jointcalling.JointCallingSpark \
    dpgt-<version>.jar \
    -i /path_to/vcfs.list \
    -r hg38.fasta \
    -o result -n 50 -j 256 -l chr1:1000000-2000000 
```

