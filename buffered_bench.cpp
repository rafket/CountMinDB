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
    assert(argc==3);
    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n=atoi(argv[1]), q=atoi(argv[2]);
    printf("I will insert %d random elements\n", n);
    double epsilon=M_E/(10*n), delta=1/pow(M_E, 3);
    printf("epsilon: %lf, delta: %lf\n", epsilon, delta);
    CountMin *cm_buffered = new CountMin(epsilon, delta, 1337, BufferedHash, n);
    map<uint64_t, int> arr;

    for(int i=0; i<n; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = dist(mt);
        cm_buffered->update(key, value);
        arr[key] += value;
    }

    vector<int> cnt_error_buffered(20);
    for(int i=0; i<q; ++i) {
        uint64_t key = dist(mt);
        uint32_t err_buffered = abs(cm_buffered->pointQuery(key)-arr[key]);
        if(err_buffered>0) {
            err_buffered/=((uint32_t)-1)/cnt_error_buffered.size()/2;
        }
        ++cnt_error_buffered[err_buffered];
    }
    print_histogram(cnt_error_buffered);
}
