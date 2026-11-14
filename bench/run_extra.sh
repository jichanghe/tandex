#!/bin/bash
# extra tools: tantan (default window and -w 2000), longdust, SRF
cd /home/claude/bench
T=tm.py
python3 $T tantan-w2000 sim2.times "tantan -f 4 -w 2000 sim2.fa > sim2.tantan2k.txt"
cd ../real
python3 ../bench/$T tantan chr21.times "tantan -f 4 chr21.fa > chr21.tantan.txt"
python3 ../bench/$T longdust chr21.times "../ext/longdust/longdust chr21.fa > chr21.longdust.bed"
cd ../ath
python3 ../bench/$T tantan-w200 cen1.times "tantan -f 4 -w 200 cen1.fa > cen1.tantan.txt"
python3 ../bench/$T longdust cen1.times "../ext/longdust/longdust cen1.fa > cen1.longdust.bed"
mkdir -p kmc_tmp
python3 ../bench/$T SRF-kmc cen1.times "kmc -fm -k151 -t1 -ci20 -cs100000 cen1.fa cen1.kmc kmc_tmp >/dev/null 2>&1 && kmc_dump cen1.kmc cen1.count.txt && ../ext/srf/srf -p cen1 cen1.count.txt > cen1.srf.fa"
echo done > /home/claude/bench/extra.done
