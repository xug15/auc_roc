## remove ploy A
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

```sh
## STAR index in /home/xugang/index/tair.star
clean=
out_dir=
STAR --runThreadN 4 --runMode alignReads --genomeDir /home/xugang/index/tair.star  --outFilterMismatchNmax 3 --outFilterMultimapNmax 8 --chimScoreSeparation 10 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --readFilesIn $dir/c-remove-rRNAs/fq/$data_name.rmrRNA.rmsnoRNA.tRNA.fq --outFileNamePrefix $dir/d-star/$data_name

```



