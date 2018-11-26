#ifndef CM_H
#define CM_H

#include <bits/stdc++.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include "MurmurHash3.h"

using namespace std;

struct Hashtable {
    uint32_t seed;
    vector<pair<uint64_t, int>> arr;

    Hashtable(uint32_t n, uint32_t seed_) : seed(seed_) {
        arr.resize(n, make_pair((uint64_t)-1, 0));
    }

    uint32_t findIndex(uint64_t key) const {
        uint32_t n = arr.size(), idx = hash(key) % n, cnt = 0;
        while (arr[idx].first != (uint64_t)-1 && arr[idx].first != key && cnt < n) {
            cnt++;
            idx++;
            idx -= (idx == n) * n;
        }
        return cnt == n ? (uint32_t)-1 : idx;
    }

    // Return a reference to the value at key. If key is new then assert that there is a slot for the new key.
    int &operator[](uint64_t key) {
        auto idx = findIndex(key);
        assert(idx != (uint32_t)-1); 
        if (arr[idx].first == (uint64_t)-1) {
            arr[idx].first = key;
        }
        return arr[idx].second;
    }

    // Return value at key if key exists otherwise -1. Assumes values are non-negative.
    int find(uint64_t key) const {
        auto idx = findIndex(key);
        if (idx == (uint32_t)-1 || arr[idx].first == (uint64_t)-1) {
            return -1;
        }
        return arr[idx].second;
    }
    
    uint64_t hash(uint64_t key) const {
        uint64_t out[2];
        MurmurHash3_x64_128(&key, sizeof(uint64_t), seed, out);
        return out[1];
    }

    uint64_t getMem() const {
        return sizeof(seed) + arr.size() * sizeof(arr.front());
    }
};

class CountMin {
public:
    CountMin(double eps, double delta, uint64_t seed, bool is_sparse, const char* filename, string name);
    void update(uint64_t i, int c);                         // increase a_i by c
    int pointQuery(uint64_t i) const;                       // return a_i
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

    const vector<Hashtable> &getSparseCounts() const {
        return sparse_counts;
    }

    bool isSparse() const {
        return is_sparse;
    }

    uint64_t getMem() {
        if (!is_sparse) {
            return w * d * sizeof(int)  /* size of the CM table */
                + d * sizeof(uint32_t)  /* size of seed for hash function */;
        }
        uint64_t out = 0;
        for (size_t i = 0; i < d; ++i) {
            out += sparse_counts[i].getMem();  // size of each hashtable
        }
        out += d * sizeof(uint32_t);
        return out;
    }

private:
    int* flatcounts;
    int** counts;
    vector<Hashtable> sparse_counts;
    vector<uint32_t> hash_seed;
    uint64_t w, d;
    bool is_sparse;
    char* filename;
    string name;

    uint64_t hash(uint64_t key, uint32_t seed) const {
        uint64_t out[2];
        MurmurHash3_x64_128(&key, sizeof(uint64_t), seed, out);
        return out[1];
    }
};

CountMin::CountMin(double eps, double delta, uint64_t seed, bool is_sparse = false, const char* filename = NULL, string name = "noname") {
    this->is_sparse = is_sparse;
    this->name = name;
    w = ceil(M_E / eps);
    d = ceil(log(1 / delta));
    printf("%s : w: %llu, d: %llu\n", name.c_str(), w, d);
    if (!is_sparse) {
        if (filename == NULL) {
            this->filename = NULL;
            flatcounts = (int*)calloc(d * w, sizeof(int));
        } else {
            // TODO: implement this
            this->filename = strdup(filename);
            int fd = open(filename, O_RDWR | O_CREAT, 0666);
            ftruncate(fd, d * w * sizeof(int));
            flatcounts = (int*)mmap(0, d * w * sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_NOCACHE, fd, 0);
        }
        counts = (int**)malloc(d * sizeof(int*));
        for (size_t i = 0; i < d; ++i) {
            counts[i] = flatcounts + i * w;
        }
    } else {
        assert(filename == NULL);
        // TODO: check the correctness: since the keys are from 0 to w-1, the size should be 2*w?
        sparse_counts.resize(d, Hashtable(2 * w, 1337));
    }
    hash_seed.resize(d);
    mt19937_64 mt(seed);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    for (size_t i = 0; i < d; i++) {
        hash_seed[i] = dist(mt);        // generate random seeds for hash function
    }
}

void CountMin::update(uint64_t i, int c) {
    if (!is_sparse) {
        for (size_t j = 0; j < d; j++) {
            counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    } else {
        for (size_t j = 0; j < d; ++j) {
            sparse_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    }
}

int CountMin::pointQuery(uint64_t i) const {
    int res = 0;
    if (!is_sparse) {
        res = counts[0][hash(i, hash_seed[0]) % w];
        for (size_t j = 1; j < d; j++) {
            res = min(res, counts[j][hash(i, hash_seed[j]) % w]);
        }
    } else {
        auto it = sparse_counts[0].find(hash(i, hash_seed[0]) % w);
        res = max(0, it);
        for (size_t j = 1; j < d; j++) {
            auto it = sparse_counts[j].find(hash(i, hash_seed[j]) % w);
            res = min(res, max(0, it));
        }
    }
    return res;
}

int CountMin::innerProductQuery(const CountMin &other) const { // sparse unimplemented
    const int** otherCounts = other.getCounts();
    int res = -1;
    if (is_sparse) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    for (size_t j = 0; j < d; j++) {
        int current_res = 0;
        for (size_t i = 0; i < w; i++) {
            current_res += counts[j][i] * otherCounts[j][i];
        }
        if (res == -1 || res > current_res) {
            res = current_res;
        }
    }
    return res;
}

void CountMin::mergeCMs(const CountMin& other) {
    if (is_sparse) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }

    if (!other.is_sparse) {
        const int** otherCounts = other.getCounts();
        for (size_t i = 0; i < d; ++i) {
            for (size_t j = 0; j < w; ++j) {
                counts[i][j] += otherCounts[i][j];
            }
        }
    } else {
        const auto &otherSparseCounts = other.getSparseCounts();
        for (size_t i = 0; i < d; ++i) {
            for (const auto &it: otherSparseCounts[i].arr) {
                if (it.first == (uint64_t)-1) {
                    continue;
                }
                counts[i][it.first] += it.second;
            }
        }
    }
}

#endif
