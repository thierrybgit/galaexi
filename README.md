for compiling: 
in the galaexi cloned directory, type 

     ./compile_galaexi.sh

for running test case with rocprof 
go to Tests directory
 1) check the batch file (project number...) and
 2) type

    sbatch job_rocprof.sh

    4 runs are performed

    a) with --stats (outputs in stat-only.stats.csv are correct and similar to Nvidia outputs)

    b) with --stats --hsa-trace (outputs in stat+hsatrace.stats.csv are similar to in stat-only.stats.csv, i.e. they are correct)

    c) with --stats --hip-trace (outputs in stat+hiptrace.stats.csv are DIFFERENT to in stat-only.stats.csv, i.e. they are incorrect: see top five kernels Name and # of calls!)

    d) with --stats --sys-trace (outputs in stat+systrace.stats.csv are DIFFERENT to in stat-only.stats.csv, i.e. they are incorrect: see top five kernels Name and # of calls!)
