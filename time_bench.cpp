#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

int main(int argc, char** argv) {
    if(argc!=3) {
        printf("Usage:\n    %s [number of elements] [number of queries]\n", argv[0]);
    }
    // cm should be a factor of 10 smaller
    // epsilon = e/(10*n), additive error of e/10
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), q=atoi(argv[2]);
    printf("I will insert %d random elements and perform %d random queries\n", n, q);
    double epsilon=M_E/(10*n), delta=1/pow(M_E, 3);
    printf("epsilon: %lf, delta: %lf\n", epsilon, delta);
    CountMin cm_normal(epsilon, delta, 1337, false);
    CountMin cm_sparse(epsilon, delta, 1337, true);
    vector<pair<uint64_t, uint64_t> > arr(n);
    clock_t start, finish;

    start=clock();
    for(int i=0; i<n; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_normal.update(key, value);
    }
    finish=clock();
    printf("building normal countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("normal countMin takes up %lfMB\n", (double)cm_normal.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<n; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        cm_sparse.update(key, value);
    }
    finish=clock();
    printf("building sparse countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("sparse countMin takes up %lfMB\n", (double)cm_sparse.getMem()/1024/1024);

    start=clock();
    for(int i=0; i<n; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        arr[i] = {key, value};
    }
    finish=clock();
    printf("building buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("buffer takes up %lfMB\n", (double)arr.size()*(sizeof(uint64_t)+sizeof(uint64_t))/1024/1024);

    uint64_t sum1=0, sum2=0, sum3=0;
    start=clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_normal.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying normal countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_sparse.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying sparse countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        uint64_t key = dist(mt);
        for(int j=0; j<n; ++j) {
            if(key==arr[j].first) {
                sum3 += arr[j].second;
            }
        }
    }
    printf("pay no attention to these numbers: %lu %lu %lu\n", sum1, sum2, sum3);
    finish=clock();
    printf("querying buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    cm_normal.mergeCMs(cm_sparse);
    finish=clock();

    printf("Merging took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
}
