#建立文件夹
cd /local_data1/GUEST/3d_train/TrainDir
mkdir test
cd test

mkdir atac-seq
cd atac-seq

#添加环境变量
#export PATH=/local_data1/GUEST/3d_train/software/miniconda2/bin/:$PATH

#第一步：去接头
mkdir data
cd data

cp /local_data1/GUEST/3d_train/rawdata/ATAC-seq/data/Rep1_R?.fq.gz .  #这里的?代表任意一个字符
cp /local_data1/GUEST/3d_train/rawdata/ATAC-seq/data/Rep2_R?.fq.gz .

cutadapt -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -q 10 --trim-n -m 10 -o Rep1_trim.R1.fq -p Rep1_trim.R2.fq Rep1_R1.fq.gz Rep1_R2.fq.gz
cutadapt -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -q 10 --trim-n -m 10 -o Rep2_trim.R1.fq -p Rep2_trim.R2.fq Rep2_R1.fq.gz Rep2_R2.fq.gz

fastqc -o . Rep1_trim.R1.fq -p Rep1_trim.R2.fq Rep2_trim.R1.fq -p Rep2_trim.R2.fq

#第二步 建立index
cd ..
mkdir index
cd index

cp /local_data1/GUEST/3d_train/rawdata/ATAC-seq/index/ref.fa .
#bowtie2-build ref.fa ref
#samtools faidx ref.fa
cp /local_data1/GUEST/3d_train/rawdata/ATAC-seq/index/* .

#第三步 比对
cd ../data/
bowtie2 -x ../index/ref -x 1000 -1 Rep1_trim.R1.fq -2 Rep1_trim.R2.fq -s Rep1.sam
bowtie2 -x ../index/ref -x 1000 -1 Rep2_trim.R1.fq -2 Rep2_trim.R2.fq -s Rep2.sam

cp /local_data1/GUEST/3d_train/ATAC-seq/data/Rep1.sam .
cp /local_data1/GUEST/3d_train/ATAC-seq/data/Rep2.sam .

#第四步 挑选可靠的比对结果
samtools view -b -f 2 -q 30 -o Rep1.pairs.bam Rep1.sam
samtools view -b -f 2 -q 30 -o Rep2.pairs.bam Rep2.sam
#samtools view Rep1.pairs.sort.bam |less -S

#第五步 去除PCR重复
samtools sort -o Rep1.pairs.sort.bam Rep1.pairs.bam
samtools sort -o Rep2.pairs.sort.bam Rep2.pairs.bam

java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=Rep1.pairs.sort.bam O=Rep1.pairs.sort.dedup.bam M=Rep1.duplicates.log
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=Rep2.pairs.sort.bam O=Rep2.pairs.sort.dedup.bam M=Rep2.duplicates.log

#第六步 查看线粒体污染
samtools index Rep1.pairs.sort.dedup.bam
samtools index Rep2.pairs.sort.dedup.bam

#第七步 去除线粒体污染。假定线粒体序列名为chrM，叶绿体为Pt
samtools view -h Rep1.pairs.sort.dedup.bam |grep -v "chrM" |grep -v "Pt" |samtools view -bS -o Rep1.final.bam
samtools view -h Rep2.pairs.sort.dedup.bam |grep -v "chrM" |grep -v "Pt" |samtools view -bS -o Rep2.final.bam

#第八步 统计片段插入长度分布
plotinsertSize Rep1.final.bam Rep1
plotinsertSize Rep2.final.bam Rep2

#第九步 peak calling
bedtools bamtobed -i Rep1.final.bam >Rep1.final.bed
bedtools bamtobed -i Rep2.final.bam >Rep2.final.bed

macs2 callpeak -t Rep1.final.bed -n Rep1 --shift -100 --extsize 200 --nomodel -B SPMR -g 2.4e8
macs2 callpeak -t Rep2.final.bed -n Rep2 --shift -100 --extsize 200 --nomodel -B SPMR -g 2.4e8

#第十步 TSS富集
#replicate 1
bedSort Rep1_treat_pileup.bdg stdout |bedClip -truncate stdin ../index/ref.fa.fai stdout |perl -ane 'print if($F[1]<$F[2])' >Rep1_treat_pileup.bedGraph
bedGraphToBigWig Rep1_treat_pileup.bedGraph ../index/ref.fa.fai Rep1.bw
#replicate 2
bedSort Rep2_treat_pileup.bdg stdout |bedClip -truncate stdin ../index/ref.fa.fai stdout |perl -ane 'print if($F[1]<$F[2])' >Rep2_treat_pileup.bedGraph
bedGraphToBigWig Rep2_treat_pileup.bedGraph ../index/ref.fa.fai Rep2.bw

computeMatrix reference-point -S Rep1.bw Rep2.bw -R ../index/genes.gtf -a 3000 -b 3000 -p 1 -o tss3k.matrix.gz
plotHeatmap -m tss3k.matrix.gz -o tss3k.heatmap.png --colorMap Reds

#第十一步 peak的重复性
bedtools intersect -a Rep1_peaks.narrowPeak -b Rep2_peaks.narrowPeak |wc -l
bedtools intersect -a Rep1_peaks.narrowPeak -b Rep2_peaks.narrowPeak -f 0.5 -F 0.5 |wc -l

idr -s Rep1_peaks.narrowPeak Rep2_peaks.narrowPeak -o idr_result --plot

#第十二步 peak在基因上的分布
Rscript PeakAnnotation.R ../index/genes.gtf Rep1_peaks.narrowPeak Rep1
Rscript PeakAnnotation.R ../index/genes.gtf Rep2_peaks.narrowPeak Rep2





