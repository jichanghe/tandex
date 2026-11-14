#!/bin/bash
# run_group.sh g1  : rebuild group from staged parts, run tandex per chromosome
g=$1; U=/mnt/user-data/uploads/c/chm13/chr; cd /home/claude/chm13
cat $U/$g.part.* > $g.fa.gz && gzip -t $g.fa.gz || { echo "$g gzip FAILED" >> wg.log; exit 1; }
zcat $g.fa.gz | awk -v g=$g '/^>/{f=substr($1,2)".fa"; print $1 > f; print substr($1,2) > (g".list"); next}{print > f}'
rm -f $g.fa.gz
for c in $(cat $g.list); do
  [ "$c" = chrM ] && { rm -f $c.fa; continue; }
  python3 ../bench/tm.py tandex-p50000-t2-$c wg.times "../tool/tandex -t 2 -p 50000 $c.fa > $c.wg.tsv"
  rm -f $c.fa
done
echo "$g done" >> wg.log
