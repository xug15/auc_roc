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
    $length=length($seq);
    if($length < 16){
	#print "$length\t$seq\n";
	next;
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

## Determine the P-site position of Ribo-seq. [P-site_determination.sh]
### a5-P-site_determination.sh

```sh
out=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/
bam=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/d-bam/mapAligned.sortedByCoord.out.bam
bed=/home/xugang/index/tair_ribowave/start_codon.bed
script=/home/xugang/RiboWave_v1.0/script/
name=Col

/home/xugang/RiboWave_v1.0/script/P-site_determination.sh -i $bam -S $bed -o $out -s $script -n $name
```

## Generate P-site file
### a6-P-site-create.py
```python
import sys
data=open(sys.argv[1])
out=open('P-site/psite1nt.txt','w')
for d in data:
    d=d.rstrip("\n")
    num=d.split("\t")
    name=num.pop(0)
    num=list(map(int,num))
    index=num.index(max(num))+1
    name_l=name.split(' ')[0]
    out.write(str(name_l)+'\t'+str(index)+"\n")
```

## Generating P-site track
### a7-p-site-track
```sh
out=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/f-P-site-track
bam=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/d-bam/mapAligned.sortedByCoord.out.bam
genome_size=/home/xugang/index/tair.annotation/tair-genome-size.txt
psite_position=/home/xugang/data_guoruixin-20190302/Col-0_FKDL171663936-1A/e-p-site/P-site/psite1nt.txt
name=Col
gtf=/home/xugang/index/tair_ribowave/exons.gtf
script=/home/xugang/RiboWave_v1.0/script/

[[ -d f-P-site-track ]] || mkdir f-P-site-track

/home/xugang/RiboWave_v1.0/script/create_track_Ribo.sh -i $bam -G $gtf -g $genome_size -P $psite_position -o $out -s $script -n $name

```

## Denoise the P-site track
### a8-Denoise-P-site-track.sh
```sh
psite=f-P-site-track/bedgraph/Col/final.psite
ORFs=/home/xugang/index/tair_ribowave/final.ORFs
script=/home/xugang/RiboWave_v1.0/script/

mkdir -p g-Ribowave;
/home/xugang/RiboWave_v1.0/script/Ribowave  -a $psite -b $ORFs -o g-Ribowave -n Col -s $script -p 8
```

## Identifying translated ORF
### a9-Identifying-translated-ORF.sh
```sh
psite=f-P-site-track/bedgraph/Col/final.psite
ORFs=/home/xugang/index/tair_ribowave/final.ORFs
script=/home/xugang/RiboWave_v1.0/script/
name=Col

mkdir -p h-Identifying-translated-ORF;
/home/xugang/RiboWave_v1.0/script/Ribowave -P -a $psite -b $ORFs -o h-Identifying-translated-ORF -n $name -s $script -p 4

```

## Estimating abundance
### a10-Estimating-abundance.sh
```sh
psite=f-P-site-track/bedgraph/Col/final.psite
ORFs=/home/xugang/index/tair_ribowave/final.ORFs
script=/home/xugang/RiboWave_v1.0/script/
name=Col

mkdir -p i-Estimating-abundance
/home/xugang/RiboWave_v1.0/script/Ribowave -D -a $psite -b $ORFs -o i-Estimating-abundance -n $name -s $script -p 4;
```

## Estimating TE
**IMPORTANT : when estimating TE, user should input the sequenced depth of Ribo-seq and the FPKM value from paired RNA-seq**
### a11-Estimating-TE.sh
```sh
psite=f-P-site-track/bedgraph/Col/final.psite
ORFs=/home/xugang/index/tair_ribowave/final.ORFs
script=/home/xugang/RiboWave_v1.0/script/
name=Col


mkdir -p j-Estimating-TE;
/home/xugang/RiboWave_v1.0/script/Ribowave -T 9012445  GSE52799/mRNA/SRR1039761.RPKM -a $psite -b $ORFs -o j-Estimating-TE -n $name -s $script -p 4

```





