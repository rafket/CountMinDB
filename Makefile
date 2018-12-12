CXX = g++

time_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h time_bench.cpp
	$(CXX) -g lz4.c MurmurHash3.cpp time_bench.cpp -lz -o time_bench -Wall

error_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h error_bench.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp error_bench.cpp -lz -o error_bench -Wall

ssd_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h ssd_bench.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp ssd_bench.cpp -lz -o ssd_bench -Wall

time_ip_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h time_ip_bench.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp time_ip_bench.cpp -lz -o time_ip_bench -Wall

time_btw_ip_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h time_btw_ip_bench.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp time_btw_ip_bench.cpp -lz -o time_btw_ip_bench -Wall

buffered_bench: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h buffered_bench.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp buffered_bench.cpp -lz -o buffered_bench -Wall

experiments: time_bench error_bench ssd_bench

clean:
	rm time_bench error_bench ssd_bench time_ip_bench time_btw_ip_bench pq_test totaltest spacetest

run_time_bench: time_bench
	./time_bench 200000 200000 10000 && ./time_bench 400000 400000 10000 && ./time_bench 600000 600000 10000 && ./time_bench 800000 800000 10000 && ./time_bench 1000000 100000 10000

run_error_bench: error_bench
	./error_bench 100000 10000 100000

run_ssd_bench: ssd_bench
	./ssd_bench 100000 10000 100000

run_buffered_bench: buffered_bench
	./buffered_bench 10000000 1000

run_experiments: run_time_bench run_error_bench run_ssd_bench

pq_test: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h pq_test.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp pq_test.cpp -lz -o pq_test -Wall

run_pq_test: pq_test
	./pq_test 100000 10000 100000

totaltest: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h zipf.h totaltest.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp totaltest.cpp -lz -o totaltest -Wall

run_totaltest: totaltest
	./totaltest

spacetest: lz4.h lz4.c MurmurHash3.h MurmurHash3.cpp cm.h zipf.h spacetest.cpp
	$(CXX) -O3 lz4.c MurmurHash3.cpp spacetest.cpp -lz -o spacetest -Wall

run_spacetest: spacetest
	./spacetest
