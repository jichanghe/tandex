#!/usr/bin/env python3
# tm.py <label> <timesfile> <shell command>: append wall time and peak RSS of the command
import sys, time, subprocess, resource
t = time.time(); rc = subprocess.call(sys.argv[3], shell=True)
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
open(sys.argv[2], 'a').write(f"{sys.argv[1]}\t{time.time()-t:.1f} s\t{ru.ru_maxrss/1024:.0f} MB\trc={rc}\n")
