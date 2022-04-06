echo "Fichier graph"
read graph

result_perf="perf_Gorder.txt"


	perf stat -B -e cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-store-misses,L1-dcache-prefetch-misses,LLC-loads,LLC-stores,dTLB-load-misses,dTLB-store-misses,cycles -o $result_perf --append    taskset -c 0-0 ./Gorder $graph
