### sc_mut_pipline

SeekOne FAST break the limitations of the conventional 3 'or 5' end transcriptome, and use random primers for random capture and detection of the full transcriptome, resulting in more abundant transcriptome information in the sequencing data. Therefore, the use of full sequence data for variation analysis of single-cell transcriptome data has become a new analysis content. At the same time, for some samples (such as tumor samples), insufficient depth of full sequence data may result in missing detection of some loci. In order to improve the accuracy of gene detection, the corresponding panel enrichment of the full sequence library can be carried out.
In view of the FAST sequence or data variation analysis, after the enrichment of this document reference souporcell software (https://github.com/wheaton5/souporcell) to construct the following pipline.

![图片.png](attachment:88965d55-4f99-4c42-bb2a-55ebcaa4104a.png)

### requirements

seeksoultools: http://seeksoul.seekgene.com/index.html 

STAR: https://github.com/alexdobin/STAR

hisat: https://github.com/DaehwanKimLab/hisat2

freebayes: https://github.com/freebayes/freebayes

bcftools: https://github.com/samtools/bcftools

VEP: http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html

VEP database: http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

vatrix: https://github.com/10XGenomics/vartrix

The above software can be downloaded and installed separately or integrated into one environment through conda installation.

```
ref="hisat_idx/hg38.hisat2.genome" （以human hg38为例）
star_index="star_index"
fq1="fastq.R1.gz"
fq2="fastq.R2.gz"
vepcache="software/vep-database"
barcode="barcodes.tsv"
out="outpupath"
```

### Usage

#### 1、extract cell barcode and umi by seeksoultools

```
seeksoultools fast step1 \
    --fq1 $fq1  \
    --fq2 $fq2  \
    --chemistry DD-Q \
    --samplename $sam  \
    --outdir  $out
```

output：sam_1.fq.gz sam_2.fq.gz are paired reads containing valid barcode, respectively.  sam_multi_1.fq.gz and sam_multi_2.fq.gz indicates the paired reads  that matches multiple barcode whitlists during barcode identification; sam_multi.json is the statistics file about the data composition.


#### 2、mapping by hisat2 and STAR

```
fq1=${out}/${sam}/step1/${sam}_1.fq.gz
fq2=${out}/${sam}/step1/${sam}_2.fq.gz

star --runThreadN 13 --limitOutSJcollapsed  5000000 --readMapNumber -1 --genomeDir ${star_index} --readFilesCommand zcat --outFileNamePrefix star --outSAMtype BAM Unsorted --readFilesIn $fq1 $fq2

samtools sort -O BAM -o sort.star.bam starAligned.out.bam

seeksoultools utils addtag --inbam sort.star.bam  --outbam star.tagged_sorted.bam

samtools index star.tagged_sorted.bam

# hisat
hisat2  -t 8 -1 $fq1 -2 $fq2 -x  $ref | samtools sort  -O BAM -o hisat.bam -

seeksoultools utils addtag --inbam hisat.bam  --outbam hisat.tagged_sorted.bam

samtools index hisat.tagged_sorted.bam
```

```
输出说明：starAligned.out.bam 是STAR软件直接比对输出的bam 
        sort.star.bam  是对STAR输出bam按照染色体位置进行排序后的bam
        star.tagged_sorted.bam 是在bam中增加barcode和UMI信息，tag中的“CB” “UB”分别表示barcode和UMI序列
        star.tagged_sorted.bam.bai 索引文件
        hisat.bam 是hisat2软件进行比对并按照染色体位置排序后的bam
        hisat.tagged_sorted.bam 在hisat比对完的bam中增加了barcode和UMI信息的bam
        hisat.tagged_sorted.bam.bai 索引文件
```

#### 3、calling candidate variants by freebayes and normlize was performed by bcftools

```
# star
freebayes -f ${ref}.fa -Xu -C 2 -q 20 -n 3 -E 40 -m 30 --min-coverage 6 --min-alternate-fraction 0.03  star.tagged_sorted.bam > star.freebayes.vcf

bcftools norm -m +both -Ov -f ${ref}.fa star.freebayes.vcf > norm.star.freebayes.vcf

# hisat

freebayes -f ${ref}.fa -Xu -C 2 -q 20 -n 3 -E 40 -m 30 --min-coverage 6 --min-alternate-fraction 0.03  hisat.tagged_sorted.bam > hisat.freebayes.vcf

bcftools norm -m +both -Ov -f ${ref}.fa hisat.freebayes.vcf > norm.hisat.freebayes.vcf
```

```
输出说明：star.freebayes.vcf 是STAR软件比对结果的变异信息文件
        norm.star.freebayes.vcf 是对变异位点进行变异标准化后的vcf文件
        hisat.freebayes.vcf 是hisat2软件比对结果的变异信息文件 
        norm.hisat.freebayes.vcf 是对变异位点进行变异标准化后的vcf文件
```

#### 4、predicting the functional effects of genomic variants by VEP

```
# star

vep -i norm.star.freebayes.vcf -o star.freebayes.vep.vcf --species homo_sapiens --cache    --force --dir_cache $vepcache --sift b --polyphen b --ccds  --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna  --merged  --exclude_predicted --tab --gencode_basic --per_gene --pick_allele --vcf

# hisat

vep -i norm.hisat.freebayes.vcf -o hisat.freebayes.vep.vcf --species homo_sapiens --cache    --force --dir_cache $vepcache --sift b --polyphen b --ccds  --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --max_af --pubmed --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna  --merged  --exclude_predicted --tab --gencode_basic --per_gene --pick_allele --vcf


```

```
输出说明：star.freebayes.vep.vcf 是STAR软件比对结果的变异信息注释结果文件
        star.freebayes.vep.vcf_summary.html 是VEP软件注释结果统计信息
        hisat.freebayes.vep.vcf 是hisat2软件比对结果的变异信息注释结果文件
        hisat.freebayes.vep.vcf_summary.html 是VEP软件注释结果统计信息
```

#### 5、intersect and filter variants 

```
python filter.py --minimapvcf star.freebayes.vep.vcf --hisatvcf hisat.freebayes.vep.vcf --finalvcf filter.vcf
```

Filter conditions include: Non-synonymous mutations (including "frameshift_variant","inframe_deletion","inframe_insertion","missense_variant","start_lost","stop_gained"); sor value less than 3; Overlay read greater than 10; The number of supporting mutation reads was greater than 3; EAS_AF does not exist in the EAS_AF database or the frequency is less than 0.01.

#### 6、cell allele counting by Vartrix

```
mkdir matrix

vartrix -v filter.vcf  -b star.tagged_sorted.bam -f ${ref}.fa -c  ${barcode}  -o matrix/alt_matrix.mtx   --out-variants matrix/variants.txt   --scoring-method coverage   --ref-matrix matrix/ref_matrix.mtx   --umi  --threads 4    --log-level debug --mapq 30

Rscript tomatrix.R  --refmatrix matrix/ref_matrix.mtx --altmatrix matrix/alt_matrix.mtx --annofile filter.vcf.txt --barcode ${barcode} --variants matrix/variants.txt --outdir ${outdir}
```

In this step, the matrix is further filtered to retain only sites where the number of mutant cells is greater than 2 and the total mutation UMI is greater than 2.

#### 7、fisher-test and plot

```
Rscipt snv_module_v2.R
    --rna_mat seurat.rds \
    --name ${sam} \
    --outdir  ${outdir}/mutation_umap \
    --snv_cover_mat  final.all.matrix \
    --snv_mut_mat final.alt.matrix
```

```
结果说明：sam_snv_markers.xls文件是位点在cluster中做的fisher检验
         SNV：突变位点信息
         p_val：fisher检验的p值
         ident1_cover：目标cluster/celltype中有多少细胞覆盖到了这个突变位点但是没发生突变
         ident1_mut：目标cluster/celltype中有多少细胞带有这个突变
         ident2_cover: 非目标cluster/celltype中有多少细胞覆盖到了这个突变位点但是没发生突变
         ident2_mut：非目标cluster/celltype中有多少细胞带有这个突变
         cluster：目标细胞类型
         
```

![图片.png](attachment:15498b1b-8392-437e-a359-a0a9405e024c.png)

```
Top left: Each dot represents a cell, and the colored dot represents the sequence of that cell covering the mutation site.
Bottom left: Each dot represents a cell, and the red dot represents the cell with the mutation.
Top right: Each dot represents a cell, and the shade of blue represents the expression level of the mutant gene at the RNA level. The expression of RNA is shown in order to rule out whether the cells in which the mutation was not detected were not able to detect/capture the mRNA sequence of the gene.
Bottom right: Each dot represents a cell, and different cell types are represented by different colors.
```

The rds file  adds mutation and coverage information on the basis of expressing rds, which can be directly used for subsequent analysis.

### Issues

#### 1、the software about align

For RNA sequence, the comparison software would be at the end of reads or the comparison of splice would be poor, and different software would cause different false positives, so two comparison software were used and their mutation results were intersected. Referring to souprocell software, minimap and hisat, minimap is not suitable for double-ended comparison of short reads, so it is replaced by the STAR software most commonly used in single-cell RNA data.（PMID: 32366989， 28680106 ）

#### 2、检出位点情况说明


利用我们的混合细胞系数据进行测试发现： 对于两个已知的位点：chr12:25245351-25245351（G12S）chr7:55174772-55174772(exon19:c.2236_2250del)，对于G12S这个snp位点可以正常检测到，对于后者indel位点，最终结果中不存在，其原因是hisat比对得到的bam中不存在该位点的deletion。

我们对比6个样本的seekoneFAST数据和WES数据：

| WES | ovelap |scFAST |
| :---: | :---: | :---: |
| 77 | 13（8.725%）| 59 |
| 95 | 9 （5.389%）| 63 |
|77|16（10.323%） |62
|112|12（6.704%）|55|
|68|15（11.194%）|51|
|76|7（4.828%）|62|


scFAST特有的位点，在WES的bam中可以看到位点覆盖度不高，可能是因为突变丰度较低，WES深度不够，也可能位点是是在转录层面的突变；
|summary|value|
|:---: | :---: |
|Min.  :|8|
|1st Qu.:|28|
|Median :|45|
| Mean   :|54.75|
|3rd Qu.:|80|
|Max.   :|194|

WES特有的位点在scFAST数据中可看到，一些存在链偏好性（SOR>3），一些不存在突变或alt reads极少，该类型位点在WES数据中vaf也不太高。
|summary|value|
|:---: | :---: |
|Min.  :|0.03200|
|1st Qu.:|0.05412|
|Median :|0.09375|
| Mean   :|0.21360|
|3rd Qu.:|0.40228|
|Max.   :|0.65395|

本文档旨在提供一个scFAST检测变异的思路及参考示例，然而因为数据多样性，最终得到的结果可能会存在某些假阳性或者假阴性，可根据自己数据及结果对软件参数或者过滤条件进行调整。

#### 3、关于输入文件

构建索引时，需要分别构建STAR索引及hisat索引，请保持所用fasta文件为同一个。

下载vep database时请下载indexed_vep_cache，并推荐采用ossline模式，这样会加快运行速度，且大大减少内存消耗。vep设置的参数请与下载的database一致。


