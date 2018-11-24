#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

int main(int argc, char** argv) {
    assert(argc==3);
    // cm should be a factor of 10 smaller
    // epsilon = e/(10*n), additive error of e/10
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), q=atoi(argv[2]);
    printf("I will insert %d random elements\n", n);
    double epsilon=M_E/(4*n), delta=1/pow(M_E, 3);
    printf("epsilon: %lf, delta: %lf\n", epsilon, delta);
    CountMin cm_normal(epsilon, delta, 1337, false);
//    CountMin cm_sparse(epsilon, delta, 1337, true);
    map<uint64_t, int> arr;

    for(int i=0; i<n; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = dist(mt);
        cm_normal.update(key, value);
//        cm_sparse.update(key, value);
        arr[key] += value;
    }

    vector<int> cnt_error_normal(20), cnt_error_sparse(20);
/*    for(auto &v: arr) {
        uint32_t err = abs(cm_normal.pointQuery(v.first)-v.second);
        if(err>0) {
            err/=((uint32_t)-1)/cnt_error_normal.size()/2;
        }
        ++cnt_error_normal[err];
    }*/
    for(int i=0; i<q; ++i) {
        uint64_t key = dist(mt);
        uint32_t err = abs(cm_normal.pointQuery(key)-arr[key]);
        if(err>0) {
            err/=((uint32_t)-1)/cnt_error_normal.size()/2;
        }
        ++cnt_error_normal[err];
    }

    printf("Histogram of errors (last value is sum of errors after the max displayed error)\n");
    for(size_t i=0; i<cnt_error_normal.size(); ++i) {
        if(i==0) {
            printf("0 ");
        }
        else {
            printf("%lu--%lu ", 1+(i-1)*((uint32_t)-2)/(cnt_error_normal.size()-1)/2, 1+(i)*((uint32_t)-2)/(cnt_error_normal.size()-1)/2);
        }
    }
    printf("\n");
    for(auto &v: cnt_error_normal) {
        printf("%d ", v);
    }
    printf("\n");
    /*
    printf("Histogram of errors (last value is sum of errors after the max displayed error)\n");
    for(auto &v: cnt_error_sparse) {
        printf("%d ", v);
    }
    printf("\n");*/
}
