# DPGT
DPGT: Distributed Population Genetics tools, a Spark based high-performance joint variant calling tool for million whole genome sequencing samples.
## Introduction
- Joint variant calling of multiple samples to obtain vcf files from a set of gVCFs. DPGT is able to make variant calls to millions of populations, which is efficient and accurate and extensible.


## Built
DPGT makes genomic data work with Spark. This project is built using maven and Java 8.

The user needs to clone the project locally, then switch to the root path of the project, and finally compile the project into a jar package. The command is shown below.
```
git clone git@github.com:BGI-flexlab/DPGT.git
cd DPGT
mvn package
```


## Get started

The main running commands are as follows:
```
spark-submit  \
  --conf "spark.network.timeout=10000000" \
  --class org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark.JointCallingSpark \
  --name DPGT DPGT.jar \
  -i file:///path/gvcf.list \
  -r file:///path/WindowIndex/reference/ref_bn.list \
  -k /path/dbsnp.vcf \
  -N 100 -n 100 -s 10.0 -S 10.0 -w 5000 \
  -o outdir \
  -c -l chr1:1000000-2000000
```
**Note: Details of the parameters are as follows.**

| opt      | longOpt     | hasArgs| required     |Description     |
| :-----: | :--------:  | :---------:| :---------: |:-----------: |
|  a | allSitePLs    | false  | false   |Annotate all sites with PLs.|
|  A | annotateNDA     | false | false |If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site").|
| b | hets    | true  | false   |Heterozygosity value used to compute prior likelihoods for any locus.|
|  B | indel_hets     | true | false |Heterozygosity for indel calling.|
|  c | mergeChrom     |  false | false|  output files of each with one chromosome data.|
|  C | sample_ploidy     | true | false |Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).|
|  G | gt_mode     | true | false |Specifies how to determine the alternate alleles to use for genotyping(DISCOVERY or GENOTYPE_GIVEN_ALLELES).|
|  i | input     | true | true |a gvcf list for input.|
|  I | include_non_variant   | false  | false|Include loci found to be non-variant after genotyping.|
|  j | heterozygosity_stdev     | true  | false|Standard deviation of eterozygosity for SNP and indel calling.|
|  k | knowSite     | true  | false|known snp/indel file,the format is VCF4.|
|  l | targetRegion     | true | false |target region to process.|
|  m | max_num_PL_values  | true  | false|  Maximum number of PL values to output.|
|  M | max_alternate_alleles   | true | false |  Maximum number of alternate alleles to genotype.|
|  N | combine     |  true | false |  core number used in combine step, default is 100.|
|  n | genotype     |  true | false|  core number used in genotype step, default is 100.|
|  o | output     |  true | true|  output directory.|
|  O | output_mode     | true  | false|  output mode(EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES).|
|  p | input_prior     |  true | false|  Input prior for calls(separation by Comma(,)).|
|  r | reference     |  true | true|  reference index(generation by GaeaIndex) file path, see note below for details.|
|  s | stand_emit_conf     |  true | false|  The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold").|
|  S | stand_call_conf     | true  | false|  The minimum phred-scaled confidence threshold at which variants should be called.|
|  u | uniquifySamples     |  false | false|  Assume duplicate samples are present and uniquify all names with '.variant' and file number index.|
|  U | useNewAFCalculator     |  false | false|  Use new AF model instead of the so-called exact model.|
|  w | keyWindow     | true  | false|  window size for key, default is 5000.|
|  W | regionSize     | true  | false|  region size per cycle process.|
|  z | keepCombine     |  false | false|  do not remove combine output.|

**Note: The reference can be obtained from the following command, where gaea-1.0.0.jar is in the resources directory. This will give you the WindowIndex directory.
```aidl
java -cp gaea-1.0.0.jar org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex ref.fa dbsnp.vcf
```


## License
This project is licensed under the GPLv3 License, see the [LICENSE](LICENSE) file for details.

## Citation
