time_bench: MurmurHash3.h MurmurHash3.cpp cm.h time_bench.cpp 
	g++ -O3 MurmurHash3.cpp time_bench.cpp -o time_bench -Wall

error_bench: MurmurHash3.h MurmurHash3.cpp cm.h error_bench.cpp 
	g++ -O3 MurmurHash3.cpp error_bench.cpp -o error_bench -Wall

clean:
	rm time_bench error_bench

run_time_bench: time_bench
	./time_bench 200000 200000 10000 && ./time_bench 400000 400000 10000 && ./time_bench 600000 600000 10000 && ./time_bench 800000 800000 10000 && ./time_bench 1000000 100000 10000

run_error_bench: error_bench
	./error_bench
	
experiments: run_time_bench run_error_bench
