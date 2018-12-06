#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

// This is testing the inner producting between the sketch in SSD and the updates sketch
int main(int argc, char** argv) {
    if(argc!=4) {
        printf("Usage:\n    %s [number of elements] [number of updates] [number of queries]\n", argv[0]);
        exit(0);
    }
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), u=atoi(argv[2]), q=atoi(argv[3]);
    printf("I will have %d elements in the original data, insert %d random elements and perform %d random queries\n", n, u, q);
    double epsilon=M_E/(1*n*n), delta=1/pow(M_E, 3);
    printf("epsilon: %lf, delta: %lf\n", epsilon, delta);
    CountMin cm_normal(epsilon, delta, 1337, Uncompressed, u, "file.cm");
    CountMin cm_hashtable(epsilon, delta, 1337, HashTable, u);
    CountMin cm_tree(epsilon, delta, 1337, Tree, u);

    vector<pair<uint64_t, uint64_t> > arr(u);
    clock_t start, finish;

    for(int i=0; i<n / 2; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        cm_normal.update(key, value);
    }
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        cm_hashtable.update(key, value);
    }
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        cm_tree.update(key, value);
    }
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        arr[i] = {key, value};
    }

    uint64_t sum1=0, sum2=0, sum3=0, sum4=0;

    start=clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_normal.innerProductQuery(cm_hashtable);
    }
    finish=clock();
    printf("querying from hash table took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_normal.innerProductQuery(cm_tree);
    }
    finish=clock();
    printf("querying from tree took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum3 += cm_normal.innerProductQueryRawLog(arr, u);
    }
    finish=clock();
    printf("querying from buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    printf("pay no attention to these numbers: %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4);

}
