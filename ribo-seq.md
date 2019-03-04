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
chmop;
$seq=<DATA>;
$qua=<DATA>;
$quality=<DATA>;
if($seq=~/^(G+)(.*)/){
$len=length($1);
print "$len";
}

}
```





