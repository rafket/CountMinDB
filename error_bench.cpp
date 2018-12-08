#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

void print_histogram(vector<int> errors) {
    for(size_t i=0; i<errors.size(); ++i) {
        if(i==0) {
            printf("0 ");
        }
        else {
            printf("%lu--%lu ", 1+(i-1)*((uint32_t)-2)/(errors.size()-1)/2, 1+(i)*((uint32_t)-2)/(errors.size()-1)/2);
        }
    }
    printf("\n");
    for(auto &v: errors) {
        printf("%d ", v);
    }
    printf("\n");
}

int main(int argc, char** argv) {
    // cm should be a factor of 10 smaller
    // epsilon = e/(10*n), additive error of e/10
    assert(argc == 4);
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), u=atoi(argv[2]), q=atoi(argv[3]);
    printf("I will insert %d random elements\n", n);
    double epsilon=M_E/(10*n), delta=1/pow(M_E, 3);
    double epsilon_u=M_E/(10*u);
    printf("epsilon: %lf, epsilon_u: %lf, delta: %lf\n", epsilon, epsilon_u, delta);
    CountMin cm_normal(epsilon, delta, 1337, Uncompressed, u);
    CountMin cm_optimized(epsilon_u, delta, 1337, Uncompressed, u);
    CountMin cm_compressed(epsilon, delta, 1337, ChunksZlib, u);
    map<uint64_t, int> arr;

    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = dist(mt);
        cm_normal.update(key, value);
        cm_optimized.update(key, value);
        cm_compressed.update(key, value);
        // printf("update: %d, %d\n", key, value);
        arr[key] += value;
    }

    vector<int> cnt_error_normal(20), cnt_error_optimized(20);
    for(int i=0; i<q; ++i) {
        uint64_t key = dist(mt);
        uint32_t err_normal = abs(cm_normal.pointQuery(key)-arr[key]);
        if(err_normal>0) {
            err_normal/=((uint32_t)-1)/cnt_error_normal.size()/2;
        }
        ++cnt_error_normal[err_normal];

        // printf("key: %d, %d\n", cm_normal.pointQuery(key), cm_compressed.pointQuery(key));


        assert(cm_normal.pointQuery(key) == cm_compressed.pointQuery(key));

        uint32_t err_optimized = abs(cm_optimized.pointQuery(key)-arr[key]);
        if(err_optimized>0) {
            err_optimized/=((uint32_t)-1)/cnt_error_optimized.size()/2;
        }
        ++cnt_error_optimized[err_optimized];
    }
    printf("Histogram of errors for normal count min (last value is sum of errors after the max displayed error)\n");
    print_histogram(cnt_error_normal);

    printf("Histogram of errors for optimized count min (last value is sum of errors after the max displayed error)\n");
    print_histogram(cnt_error_optimized);
}
