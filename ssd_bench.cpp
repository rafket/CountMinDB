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
        b = atoi(argv[2]),
        q = atoi(argv[3]);
    double epsilon = M_E/(10*n), delta = 1/pow(M_E, 3), epsilon_small = M_E/(10*b);
    CountMin cm_ssd(epsilon, delta, 1337, false, "file.cm");
    CountMin cm_sparse(epsilon, delta, 1337, true);
    CountMin cm_optimized(epsilon_small, delta, 1337, false);
    vector<pair<uint64_t, uint64_t> > arr(n);
    clock_t start, finish;

    start = clock();
    for(int i=0; i<n; ++i) {
        uint32_t key = dist(mt);
        uint32_t value = dist(mt);
        cm_ssd.update(key, value);
    }
    finish = clock();
    printf("building countMin on SSD took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<b; ++i) {
        uint32_t key = dist(mt);
        uint32_t value = dist(mt);
        cm_optimized.update(key, value);
    }
    finish = clock();
    printf("building optimized countMin took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<b; ++i) {
        uint32_t key = dist(mt);
        uint32_t value = dist(mt);
        cm_sparse.update(key, value);
    }
    finish = clock();
    printf("building sparse countMin took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<b; ++i) {
        uint32_t key = dist(mt);
        uint32_t value = dist(mt);
        arr[i] = {key, value};
    }
    finish = clock();
    printf("building buffer took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    uint64_t sum1=0, sum2=0, sum3=0, sum4=0;

    start = clock();
    for(int i=0; i<q; ++i) {
        sum1 += cm_ssd.pointQuery(dist(mt));
    }
    finish = clock();
    printf("querying countMin on SSD took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<q; ++i) {
        sum2 += cm_optimized.pointQuery(dist(mt));
    }
    finish = clock();
    printf("querying optimized CountMin took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<q; ++i) {
        sum3 += cm_sparse.pointQuery(dist(mt));
    }
    finish = clock();
    printf("querying sparse countMin took %lfms\n",
            (double)(finish-start)*1000/CLOCKS_PER_SEC);

    start = clock();
    for(int i=0; i<q; ++i) {
        uint64_t key = dist(mt);
        for(int j=0; j<n; ++j) {
            if(key==arr[j].first) {
                sum3 += arr[j].second;
            }
        }
    }
    finish = clock();
    printf("querying buffer took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
    printf("do not uncomment this line: %lu %lu %lu %lu\n", sum1, sum2, sum3, sum4);

    start = clock();
    cm_ssd.mergeCMs(cm_sparse);
    finish = clock();

    printf("Merging took %lfms\n", (double)(finish-start)*1000/CLOCKS_PER_SEC);
}
