# DPGT: A Distributed Population Genetics analysis Tool

DPGT is a distributed population genetics analysis tool which enabled joint calling on millions of WGS(whole genome sequencing) samples.

### Dependency

- JDK 1.8
- boost 1.74 or later
- dependencies of [htslib](https://github.com/samtools/htslib), please read htslib's INSTALL guide to install them.
- jemalloc

### Build

Get DGPT from github.
```sh
git clone --recursive git@github.com:BGI-flexlab/DPGT.git
# switch to the tag we need
git checkout v1.3.2.0
# update submodule, we need to do this because different version of DPGT may use different versions of submodules
git submodule update --init
```

DPGT is implemented using C++ and java, first we need to compile C++ libraries using following commands.
```sh
# in DPGT root directory
mkdir build
cd build && cmake -DCMAKE_PREFIX_PATH=/path/to/boost -DJAVA_INCLUDE_PATH=/path/to/jdk1.8/include ../src/main/native/
make -j 8
```
the compiled C++ libraries will be placed in `build/lib` directory.

DPGT uses [htslib](https://github.com/samtools/htslib) as a submodule, please read htslib's INSTALL guide
to install dependencies of it.

DPGT requires JDK 1.8, we have not conducted detailed tests to determine whether it can run well under other versions of JDK.

Build DPGT jar package using maven.
```sh
# in DPGT root directory
mvn package
```
the jar package will be placed in `target` directory.

### Usage

```
# print help message
export LD_LIBRARY_PATH=/path_to_dpgt_dir/build/lib:${LD_LIBRARY_PATH}
java -jar dpgt-<version>.jar -h

DPGT: Distributed Population Genetics analysis Tools
Version: 1.3.2.0

Options:
 -i,--input <FILE>                        input gvcf list in a file, one gvcf file per line.
 -o,--output <DIR>                        output directory.
 -r,--reference <FILE>                    reference fasta file.
 -l,--target-regions <STRING>             target regions to process.
 -j,--jobs <INT>                          number of spark exculators. [10]
 -n,--num-combine-partitions <INT>        number of partitions for combining gvcfs, default value is the square root of number of input gvcf files,
                                          rounded down. [-1]
 -w,--window <INT>                        window size for each combine-genotype cycle. [300M]
 -d,--delete <Boolean>                    delete combine and genotype gvcf intermediate results. Possible values: true, false. [true]
 -s,--stand-call-conf <FLOAT>             the minimum phred-scaled confidence threshold at which variants should be called. [30.0]
    --min-variant-sites <INT>             minimum number of variant sites of small partion of region. [1]
    --dbsnp <FILE>                        dbsnp vcf file for annotation.
    --use-old-qual-calculator <Boolean>   use the old AF model. Possible values: true, false. [true]
    --heterozygosity <FLOAT>              heterozygosity value used to compute prior likelihoods for any locus. [0.001]
    --indel-heterozygosity <FLOAT>        heterozygosity for indel calling. [1.25E-4]
    --heterozygosity-stdev <FLOAT>        standard deviation of heterozygosity for SNP and indel calling. [0.01]
    --max-alternate-alleles <INT>         Maximum number of alternate alleles to genotype. [6]
    --ploidy <INT>                        ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool *
                                          Sample Ploidy). [2]
    --local                               run spark in local mode, useful for debug.
 -h,--help                                print this message and exit.
 -V,--version                             show version.
```

### Run DPGT

DPGT uses [jemalloc](https://github.com/jemalloc/jemalloc) memory allocator for reducing memory fragmentation. Please read [INSTALL.md](https://github.com/jemalloc/jemalloc/blob/dev/INSTALL.md) to install jemalloc.

**On a local computer**

```sh
export LD_LIBRARY_PATH=/path_to_dpgt_dir/build/lib:${LD_LIBRARY_PATH}
export LD_PRELOAD=/path_to_jemalloc_5.3.0/lib/libjemalloc.so
spark-submit \
    --conf spark.dynamicAllocation.enabled=false \
    --conf spark.memory.fraction=0.01 \
    --conf spark.memory.storageFraction=0.01 \
    --master local[32] \
    --driver-memory 140g --executor-cores 1 \
    --class org.bgi.flexlab.dpgt.jointcalling.JointCallingSpark \
    dpgt-<version>.jar \
    -i /path_to/vcfs.list \
    -r hg38.fasta \
    -o result -j 32 -l chr1:1000000-2000000
```

`--master local[32]` means run spark application using 32-threads.


**On a yarn cluster**

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
    -o result -j 256 -l chr1:1000000-2000000 
```

### License

Licensed under the GPLv3 License. See the [LICENSE.txt](./LICENSE.txt) file for details.