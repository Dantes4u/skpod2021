heat-3d: heat-3d.out
	mpirun -n 8 --mca shmem posix --mca opal_event_include poll --map-by :OVERSUBSCRIBE --with-ft ulfm ./heat-3d.out 10 24
heat-3d.out: heat-3d.c
	mpicc heat-3d.c -o heat-3d.out
clean:
	rm -rf *.out
