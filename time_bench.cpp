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
    double epsilon=M_E/(10*n), delta=1/pow(M_E, 3);
    double epsilon_u=M_E/(10*u);
    printf("epsilon: %lf, epsilon_u: %lf, delta: %lf\n", epsilon, epsilon_u, delta);
    CountMin cm_normal(epsilon, delta, 1337, false, "file.cm");
    CountMin cm_sparse(epsilon, delta, 1337, true);
    CountMin cm_optimized(epsilon_u, delta, 1337, false);
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
        cm_optimized.update(key, value);
    }
    finish=clock();
    printf("building optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("optimized countMin takes up %lfMB\n", (double)cm_optimized.getMem()/1024/1024);

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

    uint64_t sum1=0, sum2=0, sum3=0, sum4=0;

    start=clock();
    for(int i=0; i<q; ++i) {
        sum4 += cm_optimized.pointQuery(dist(mt));
    }
    finish=clock();
    printf("querying optimized countMin took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

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
    printf("pay no attention to these numbers: %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4);
    finish=clock();
    printf("querying buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start=clock();
    cm_normal.mergeCMs(cm_sparse);
    finish=clock();

    printf("Merging took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
}
