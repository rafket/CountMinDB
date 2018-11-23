#include <bits/stdc++.h>
#include "MurmurHash3.h"
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

using namespace std;

class CountMin {
public:
    CountMin(double eps, double delta, uint64_t seed, bool is_sparse, char* filename);
    void update(uint64_t i, int c);                              // increase a_i by c
    int pointQuery(uint64_t i) const;                            // return a_i
    int innerProductQuery(const CountMin &other) const;     // return a * b (inner product)
    void mergeCMs(const CountMin &other);

    uint64_t getWidth() const {
        return w;
    }

    uint64_t getDepth() const {
        return d;
    }

    const int** getCounts() const {
        return (const int**)counts;
    }

    const vector<map<uint64_t, int>> &getSparseCounts() const {
        return sparse_counts;
    }

    bool isSparse() const {
        return is_sparse;
    }

    uint64_t getMem() {
        if(!is_sparse) {
            return w*d*sizeof(int)+d*sizeof(uint32_t);
        }
        else {
            int out=0;
            for(size_t i=0; i<d; ++i) {
                out += sparse_counts[i].size()*(sizeof(uint64_t)+sizeof(int));
            }
            out += d*sizeof(uint32_t);
            return out;
        }
    }

private:
    int* flatcounts;
    int** counts;
    //vector<vector<int> > counts;
    vector<map<uint64_t, int> > sparse_counts;
    vector<uint32_t> hash_seed;
    uint64_t w, d;
    bool is_sparse;
    char* filename;

    uint64_t hash(uint64_t key, uint32_t seed) const {
        uint64_t out[2];
        MurmurHash3_x64_128(&key, sizeof(uint64_t), seed, out);
        return out[1];
    }
};

CountMin::CountMin(double eps, double delta, uint64_t seed, bool is_sparse = false, char* filename = NULL) {
    this->is_sparse = is_sparse;
    w = ceil(exp(1) / eps);
    d = ceil(log(1 / delta));
    printf("w: %lu, d: %lu\n", w, d);
    if(!is_sparse) {
        if(filename == NULL) {
            this->filename = NULL;
            flatcounts = (int*)calloc(d*w, sizeof(int));
        }
        else {
            this->filename = strdup(filename);
            int fd = open(filename, O_RDWR | O_CREAT, 0666);
            truncate(filename, d*w*sizeof(int));
            flatcounts = (int*)mmap(0, d*w*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        }
        counts = (int**)malloc(d*sizeof(int*));
        for(size_t i=0; i<d; ++i) {
            counts[i] = flatcounts + i*w;
        }
    }
    else {
        assert(filename == NULL);
        sparse_counts.assign(d, map<uint64_t, int>());
    }
    hash_seed.resize(d);
    mt19937_64 mt(seed);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    for(size_t i = 0; i < d; i++) {
        hash_seed[i] = dist(mt);      // generate random seeds for hash function
    }
}

void CountMin::update(uint64_t i, int c) {
    if(!is_sparse) {
        for(size_t j = 0; j < d; j++) {
            counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    }
    else {
        for(size_t j=0; j<d; ++j) {
            sparse_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    }
}

int CountMin::pointQuery(uint64_t i) const {
    int res = 0;
    if(!is_sparse) {
        res = counts[0][hash(i, hash_seed[0]) % w];
        for(size_t j = 1; j < d; j++) {
            res = min(res, counts[j][hash(i, hash_seed[j]) % w]);
        }
    }
    else {
        auto it = sparse_counts[0].find(hash(i, hash_seed[0]) % w);
        if(it == sparse_counts[0].end()) {
            res = 0;
        }
        else {
            res = it->second;
        }
        for(size_t j = 1; j < d; j++) {
            auto it = sparse_counts[j].find(hash(i, hash_seed[j]) % w);
            if(it == sparse_counts[j].end()) {
                res = min(res, 0);
            }
            else {
                res = min(res, it->second);
            }
        }
    }
    return res;
}

int CountMin::innerProductQuery(const CountMin &other) const { // sparse unimplemented
    const int** otherCounts = other.getCounts();
    int res = -1;
    if(is_sparse) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    for(size_t j = 0; j < d; j++) {
        int current_res = 0;
        for(size_t i = 0; i < w; i++) {
            current_res += counts[j][i] * otherCounts[j][i];
        }
        if(res == -1 || res > current_res) {
            res = current_res;
        }
    }
    return res;
}

void CountMin::mergeCMs(const CountMin& other) {
    if(is_sparse) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }

    if(!other.is_sparse) {
        const int** otherCounts = other.getCounts();
        for(size_t i=0; i<d; ++i) {
            for(size_t j=0; j<w; ++j) {
                counts[i][j] += otherCounts[i][j];
            }
        }
    }
    else {
        const auto &otherSparseCounts = other.getSparseCounts();
        for(size_t i=0; i<d; ++i) {
            for(auto &it: otherSparseCounts[i]) {
                counts[i][it.first] += it.second;
            }
        }
    }
}
