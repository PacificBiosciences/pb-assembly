#Depends
You need the following tools in your path:
- RACON
- PBMM2
- SAMTOOLS

#Setup
Modify the config.json. Make sure the CPU settings agree between the config and cluster config. Total CPUs for PBMM2 are split by sorting and aligning


#submit
```snakemake -j 50 --cluster-config cluster.config.sge.json --cluster "qsub -S {cluster.S} -N {cluster.N} {cluster.P} -q {cluster.Q} {cluster.CPU} -e {cluster.E} -o {cluster.O} -V" -s ccs.racon.snakemake --verbose -p --latency-wait 60```
