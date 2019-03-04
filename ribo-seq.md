## remove ploy A
### a1-remove.sh
```sh
# 1.
[[ -d b-clean-reads ]] || mkdir -p b-clean-reads
#
for i in { "Col-0_FKDL171663936-1A_1.clean.fq.gz" };do
echo $i;
echo "cutadapt -a AAAAAAAA -q 30 -m 16 --trim-n -o $i.cutadapt.fastq $i";
cutadapt -a AAAAAAAA -q 30 -m 16 --trim-n -o $i.cutadapt.fastq $i;
done;
```
## remove head G+
### a2-removeG.pl
```perl
open DATA, "<$ARGV[0]";
open OUT, ">$ARGV[0].clean.fq";
while(<DATA>){
    $seq=<DATA>;
    $qua=<DATA>;
    $quality=<DATA>;
    chomp($seq);
    chomp($quality);
            if($seq=~/^(G+)(.*)/){
            $len_g=length($1);
            $len_f=length($seq);
            $len_q=$len_f-$len_g;
            $start=$len_g;
            $seq=$2;
            $quality=substr $quality, $start, $len_q;
            }
    print OUT "$_";
    print OUT "$seq\n";
    print OUT "$qua";
    print OUT "$quality\n";
}
```

## Quality of raw reads.
### 
```sh
fastqc -t 4 Col-0_FKDL171663936-1A_1.clean.fq.gz.cutadapt.fastq.clean.fq
```
## Generate the STAR mapping index
```sh
STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /home/xugang/index/tair.star \
--genomeFastaFiles /home/xugang/index/tair.annotation/tair10.fa \
--sjdbGTFfile /home/xugang/index/tair.annotation/tair10.gtf

```
## Mapping ribo-seq data into genome with STAR.
### a3-mapping.sh
```sh
[[ -d d-bam ]] || mkdir d-bam
## STAR index in /home/xugang/index/tair.star
clean=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/b-clean-reads/Col-0_FKDL171663936-1A_1.clean.fq.gz.cutadapt.fastq.clean.fq
out_dir=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/d-bam/map
STAR --runThreadN 4 --runMode alignReads --genomeDir /home/xugang/index/tair.star  --outFilterMismatchNmax 3 --outFilterMultimapNmax 8 --chimScoreSeparation 10 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --readFilesIn $clean --outFileNamePrefix $out_dir
```

## Create the annotation file for the subsequent analysis. [create_annotation.sh]
### a4-create-psite.sh
```sh
[[ -d /home/xugang/index/tair_ribowave ]] || mkdir /home/xugang/index/tair_ribowave
/home/xugang/RiboWave_v1.0/script/create_annotation.sh -G /home/xugang/index/tair.annotation/Arabidopsis_thaliana.TAIR10.34.refine.gtf -f /home/xugang/index/tair.annotation/tair10.fa  -o /home/xugang/index/tair_ribowave  -s /home/xugang/RiboWave_v1.0/script/
```


