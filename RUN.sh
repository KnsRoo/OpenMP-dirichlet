max_threads=4
step=1

dim=(100 200)
for ((omp_threads=1;omp_threads<=$max_threads;omp_threads+=step))
do
	export OMP_NUM_THREADS=$omp_threads
	for d in ${dim[@]}
	do
		./main.exe d > "result_${omp_threads}_$d.txt"
	done
done
