#!/bin/bash
# run_all.sh <prefix>  -> runs every tool on <prefix>.fa, writes <prefix>.<tool>.bed and <prefix>.times
x=$1; TRF=../ext/TRF/build/src/trf; U=../ext/ULTRA/ultra; M=../ext/mreps/mreps
tm() { python3 tm.py "$1" $x.times "$2"; }
: > $x.times
tm tandex "../tool/tandex $x.fa > $x.tandex.tsv"; cut -f1-4 $x.tandex.tsv > $x.tandex.bed
tm TRF "$TRF $x.fa 2 7 7 80 10 50 2000 -h -ngs > $x.trf.ngs"
awk '/^@/{s=substr($1,2);next}{print s"\t"$1-1"\t"$2"\t"$3}' $x.trf.ngs > $x.trf.bed
tm ULTRA-p500 "timeout 3000 $U -p 500 -t 2 --tsv --hs --disable_summary -o $x.ultra.tsv $x.fa >/dev/null 2>&1"
tm mreps-res3 "timeout 3000 $M -fasta -res 3 -minsize 25 -noprint $x.fa > $x.mreps.txt"
tm pytrf "timeout 3000 python3 run_pytrf.py $x.fa $x.pytrf.bed 500 2>/dev/null"

awk -F'\t' 'NR>1 && $2~/^[0-9]+$/{print $1"\t"$2"\t"$3"\t"$4}' $x.ultra.tsv > $x.ultra.bed
awk -v s=sim '$2=="->"{gsub(/[<>]/,"",$6); print s"\t"$1-1"\t"$3"\t"$6}' $x.mreps.txt > $x.mreps.bed
for t in tandex trf ultra mreps pytrf; do python3 evaluate.py $x.truth.tsv $x.$t.bed $t; done > $x.eval.txt
cat $x.times $x.eval.txt
