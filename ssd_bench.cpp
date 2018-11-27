#include <bits/stdc++.h>
#include "cm.h"

int main(int argc, char** argv) {
    if(argc!=4) {
        printf("Usage:\n    %s [number of elements in ssd] [number of updates in buffer] [number of queries]\n", argv[0]);
        exit(0);
    }

    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n = atoi(argv[1]),
        u = atoi(argv[2]),
        q = atoi(argv[3]);
    double epsilon = M_E/(10*n), delta = 1/pow(M_E, 3), epsilon_u = M_E/(10*u);
    CountMin cm_ssd(epsilon, delta, 1337, false, "file.cm", "count-min on SSD");
    CountMin cm_optimized(epsilon_u, delta, 1337, false, nullptr, "optimized count-min in RAM");
    CountMin cm_sparse(epsilon, delta, 1337, true, nullptr, "sparse count-min in RAM");

        vector<pair<uint64_t, uint64_t> > arr(u);
    clock_t start, finish;

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_ssd.update(key, value);
    }
    finish=clock();
    printf("building countMin in SSD took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("SSD countMin takes up %lfMB\n", (double)cm_ssd.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_optimized.update(key, value);
    }
    finish=clock();
    printf("building optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("optimized countMin takes up %lfMB\n", (double)cm_optimized.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_sparse.update(key, value);
    }
    finish=clock();
    printf("building sparse countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("sparse countMin takes up %lfMB\n", (double)cm_sparse.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        arr[i] = {key, value};
    }
    finish=clock();
    printf("building buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("buffer takes up %lfMB\n", (double)arr.size()*(sizeof(uint64_t)+sizeof(uint64_t))/1024/1024);

    uint64_t sum1=0, sum2=0, sum3=0, sum4=0;


    start=clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_ssd.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying ssd countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum4 += cm_optimized.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_sparse.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying sparse countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum3 += queryRawLog(arr, dist(mt), u);
    }
    // printf("pay no attention to these numbers: %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4);
    finish=clock();
    printf("querying buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    // merging with the sparse count min
    start=clock();
    cm_ssd.mergeCMs(cm_sparse);
    finish=clock();

    printf("Merging sparse countmin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    // merging with the raw log
    start=clock();
    cm_ssd.mergeRawLog(arr, u);
    finish=clock();

    printf("Merging raw log took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
}
