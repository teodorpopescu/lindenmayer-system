CC = gcc
CFLAGS = -std=c99 -lm

build: build-seq build-omp build-mpi-sync build-mpi-batch build-pth build-hy

build-seq: lm_seq
build-omp: lm_omp
build-mpi-sync: lm_mpi_sync
build-mpi-batch: lm_mpi_batch
build-pth: lm_pth
build-hy: lm_hy

run-seq: lm_seq
	./lm_seq

run-omp: lm_omp
	./lm_omp

run-mpi: lm_mpi
	mpirun -np 4 ./lm_mpi

run-pth: lm_pth
	./lm_pth

run-hy: lm_hy
	mpirun -np 2 ./lm_hy1

lm_seq: lindenmayer_basic.c lindenmayer.c lindenmayer_dp.c pixmap.c
	$(CC) lindenmayer_basic.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -o lm_seq

lm_omp: lindenmayer_openmp.c lindenmayer.c lindenmayer_dp.c pixmap.c
	$(CC) lindenmayer_openmp.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -fopenmp -o lm_omp

lm_mpi_sync: lindenmayer_mpi_sync.c lindenmayer.c lindenmayer_dp.c pixmap.c
	mpicc lindenmayer_mpi_sync.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -o lm_mpi_sync

lm_mpi_batch: lindenmayer_mpi_batch.c lindenmayer.c lindenmayer_dp.c pixmap.c
	mpicc lindenmayer_mpi_batch.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -o lm_mpi_batch


lm_pth: lindenmayer_pthreads.c lindenmayer.c lindenmayer_dp.c pixmap.c
	$(CC) lindenmayer_pthreads.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -pthread -o lm_pth

lm_hy: lindenmayer_hybrid.c lindenmayer.c lindenmayer_dp.c pixmap.c
	mpicc lindenmayer_hybrid.c lindenmayer.c lindenmayer_dp.c pixmap.c $(CFLAGS) -fopenmp -o lm_hy

clean:
	rm -f lm_seq lm_omp lm_mpi_sync lm_mpi_batch lm_pth lm_hy
