#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

int main(int argc, char** argv) {
    if(argc!=4) {
        printf("Usage:\n    %s [number of elements] [number of updates] [number of queries]\n", argv[0]);
        exit(0);
    }
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), u=atoi(argv[2]), q=atoi(argv[3]);
    printf("I will have %d elements in the original data, insert %d random elements and perform %d random queries\n", n, u, q);
    double epsilon=M_E/(1*n), delta=1/pow(M_E, 3);
    double epsilon_u=M_E/(10*u);
    printf("epsilon: %lf, epsilon_u: %lf, delta: %lf\n", epsilon, epsilon_u, delta);
    CountMin cm_normal(epsilon, delta, 1337, Uncompressed, u, "file.cm");
    CountMin cm_hashtable(epsilon, delta, 1337, HashTable, u);
    CountMin cm_tree(epsilon, delta, 1337, Tree, u);
    CountMin cm_zlib(epsilon, delta, 1337, ChunksZlib, u);

    // CountMin cm_optimized(epsilon_u, delta, 1337, Uncompressed);
    vector<pair<uint64_t, uint64_t> > arr(u);
    clock_t start, finish;

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_normal.update(key, value);
    }
    finish=clock();
    printf("building normal countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("normal countMin takes up %lfMB\n", (double)cm_normal.getMem()/1024/1024);

    // start=clock();
    // for(int i=0; i<u; ++i) {
    //     uint64_t key = dist(mt);
    //     uint64_t value = 10;
    //     cm_optimized.update(key, value);
    // }
    // finish=clock();
    // printf("building optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    // printf("optimized countMin takes up %lfMB\n", (double)cm_optimized.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_hashtable.update(key, value);
    }
    finish=clock();
    printf("building hash table countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("hash table countMin takes up %lfMB\n", (double)cm_hashtable.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_tree.update(key, value);
    }
    finish=clock();
    printf("building treecountMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    // printf("tree countMin takes up %lfMB\n", (double)cm_tree.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_zlib.update(key, value);
    }
    finish=clock();
    printf("building zlib countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("zlib countMin takes up %lfMB\n", (double)cm_zlib.getMem()/1024/1024);

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

    // start=clock();
    // for(int i=0; i<q; ++i) {
    //     sum4 += cm_optimized.pointQuery(dist(mt));
    // }
    // finish=clock();
    // printf("querying optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_normal.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying normal countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_hashtable.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying hash table countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);


    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_tree.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying tree countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);


    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_zlib.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying zlib countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);


    start=clock();
    for(int i=0; i<q; ++i) {
        sum3 += queryRawLog(arr, dist(mt), u);
    }
    // printf("pay no attention to these numbers: %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4);
    finish=clock();
    printf("querying buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    // // merging with the sparse count min
    // start=clock();
    // cm_normal.mergeCMs(cm_hashtable);
    // finish=clock();

    // printf("Merging hash table countmin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);


    // // merging with the sparse count min
    // start=clock();
    // cm_normal.mergeCMs(cm_tree);
    // finish=clock();

    // printf("Merging tree countmin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);


    // // // merging with the sparse count min
    // // start=clock();
    // // cm_normal.mergeCMs(cm_zlib);
    // // finish=clock();

    // // printf("Merging zlibcountmin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    // merging with the raw log
    start=clock();
    cm_normal.mergeRawLog(arr, u);
    finish=clock();

    printf("Merging raw log took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

}
