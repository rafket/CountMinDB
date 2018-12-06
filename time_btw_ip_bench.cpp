#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

// This is testing the inner producting between the two sketches in main memory (both updates)
int main(int argc, char** argv) {
    if(argc!=5) {
        printf("Usage:\n    %s [number of elements] [number of updates] [number of queries]\n", argv[0]);
        exit(0);
    }
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), u1=atoi(argv[2]), u2=atoi(argv[3]), q =atoi(argv[4]);
    printf("I will have %d elements in the original data, insert %d random elements and then %d random elements and perform %d random queries\n", n, u1, u2, q);
    double epsilon=M_E/(1*n*n), delta=1/pow(M_E, 3);
    double epsilon_u1=M_E/(1*u1*u1);
    double epsilon_u2=M_E/(1*u2*u1);
    printf("epsilon: %lf, delta: %lf\n", epsilon, delta);
    CountMin cm_hashtable1(epsilon, delta, 1337, HashTable, u1);
    CountMin cm_hashtable2(epsilon, delta, 1337, HashTable, u2);
    CountMin cm_tree1(epsilon, delta, 1337, Tree, u1);
    CountMin cm_tree2(epsilon, delta, 1337, Tree, u2);
    CountMin cm_optimized1(epsilon_u1, delta, 1337, Uncompressed, u1);
    CountMin cm_optimized2(epsilon_u2, delta, 1337, Uncompressed, u2);


    vector<pair<uint64_t, uint64_t> > arr1(u1);
    vector<pair<uint64_t, uint64_t> > arr2(u2);
    clock_t start, finish;

    for(int i=0; i<u1; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        cm_hashtable1.update(key, value);
        cm_tree1.update(key, value);
        arr1[i] = {key, value};
        cm_optimized1.update(key, value);

        // cm_hashtable2.update(key, value);
        // cm_tree2.update(key, value);
        // arr2[i] = {key, value};
        // cm_optimized2.update(key, value);
    }
    for(int i=0; i<u2; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 1;
        cm_hashtable2.update(key, value);
        cm_tree2.update(key, value);
        arr2[i] = {key, value};
        cm_optimized2.update(key, value);
    }


    uint64_t sum1=0, sum2=0, sum3=0, sum4=0, sum5=0, sum6=0;


    start=clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_hashtable1.innerProductQuery(cm_hashtable2);
    }
    finish=clock();
    printf("querying to hashtbl-hasbtl countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_tree1.innerProductQuery(cm_tree2);
    }
    finish=clock();
    printf("querying to tree-to-tree countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

      start=clock();
    for(int i=0; i<q; ++i) {
        sum3 += cm_optimized1.innerProductQuery(cm_optimized2);
    }
    finish=clock();
    printf("querying to countMin-to-countmin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum4 += cm_hashtable1.innerProductQueryRawLog(arr2, u2);
    }
    finish=clock();
    printf("querying hash table-buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum5 += cm_tree1.innerProductQueryRawLog(arr2, u2);
    }
    finish=clock();
    printf("querying tree-buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum6 += cm_optimized1.innerProductQueryRawLog(arr2, u2);
    }
    finish=clock();
    printf("querying count-min-tolog took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    printf("pay no attention to these numbers: %lu %lu %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4, sum5, sum6);

}
