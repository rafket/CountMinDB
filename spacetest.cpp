/*
 * Total testing of the count-min structures
 *
 * Run all the tests on uniform and Zipf distributions. For each distribution,
 * run the tests to record the following: update, point query, and inner
 * product query times. Obtain times values for the following structures: naked
 * count-min on disk, count-min on disk with a count-min buffer, hashtable
 * buffer, tree buffer, etc. For point and inner product queries, obtain times
 * for intermediate steps as well.
 */

#include <bits/stdc++.h>
#include "cm.h"
#include "zipf.h"

using namespace std;

mt19937_64 mt(1337);
uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
zipf_distribution dist_zipf((uint32_t)-1, 0.99);

double clock_to_ms(clock_t t) {
    return t * 1000.0 / CLOCKS_PER_SEC;
}

void update_cm(
    CountMin &cm, const vector<pair<uint64_t, int>> &kvs
) {
    for (auto &&p : kvs) {
        auto key = p.first;
        auto val = p.second;
        cm.update(key, val);
    }
}

vector<pair<uint64_t, int>> generate_uniform_kvs(int u) {
    vector<pair<uint64_t, int>> res(u);
    for (int i = 0; i < u; i++) {
        res[i].first = dist(mt);
        res[i].second = 10;
    }
    return res;
}

vector<pair<uint64_t, int>> generate_zipf_kvs(int u) {
    vector<pair<uint64_t, int>> res(u);
    for (int i=0; i < u; i++) {
        uint32_t tmp = dist_zipf(mt), out[4];
        MurmurHash3_x64_128(&tmp, sizeof(uint32_t), 1337, out);
        res[i] = {out[0], 10};
    }
    return res;
}

int pq_full_cm(CountMin &cm, const vector<uint64_t> &keys) {
    int res = 0;
    for (auto key : keys) {
        res += cm.pointQuery(key);
    }
    return res;
}

vector<uint64_t> generate_uniform_keys(int u) {
    vector<uint64_t> res(u);
    for (int i = 0; i < u; i++) {
        res[i] = dist(mt);
    }
    return res;
}

void foo(int n, int u, double delta) {
    double epsilon = M_E / (10 * n),
        //    delta = 1 / pow(M_E, 3),
           epsilon_u = M_E / (10 * u);

    // naked count-min on disk
    CountMin cm_ssd(
        epsilon, delta, 1337, Uncompressed, u, "naked.cm",
        "naked count-min on disk"
    );

    // count-min on disk with a count-min buffer
    CountMin cm_cmbuf(
        epsilon_u, delta, 1337, Uncompressed, u, nullptr,
        "in memory count-min buffer"
    );

    // count-min on disk with a hashtable buffer
    CountMin cm_hashtable(
        epsilon, delta, 1337, HashTable, u, nullptr,
        "in memory hashtable buffer"
    );

    // count-min on disk with a tree buffer
    CountMin cm_tree(
        epsilon, delta, 1337, Tree, u, nullptr,
        "in memory tree buffer"
    );

    // count-min on disk with a raw log buffer
    CountMin cm_rawlog(
        0, 0, 1337, RawLog, u, nullptr,
        "in memory raw log buffer"
    );

    CountMin cm_zlib(
        epsilon, delta, 1337, ChunksZlib, u, nullptr,
        "in memory zlib buffer"
    );

    CountMin cm_lz4(
        epsilon, delta, 1337, LZ4, u, nullptr,
        "in memory lz4 buffer"
    );

    // update time
    auto uniform_kvs = generate_uniform_kvs(u);
    auto zipf_kvs = generate_zipf_kvs(u);

    // Need for some to figure out the size
    update_cm(cm_hashtable, uniform_kvs);
    update_cm(cm_tree, uniform_kvs);
    update_cm(cm_rawlog, uniform_kvs);
    update_cm(cm_zlib, uniform_kvs);
    update_cm(cm_lz4, uniform_kvs);

    printf("naked count-min on disk size:    %lld\n", cm_ssd.getMem());
    printf("in memory count-min buffer size: %lld\n", cm_cmbuf.getMem());
    printf("in memory hashtable buffer size: %lld\n", cm_hashtable.getMem());
    printf("in memory tree buffer size:      %lld\n", cm_tree.getMem());
    printf("in memory raw log buffer size:   %lld\n", cm_rawlog.getMem());
    printf("in memory zlib buffer size:      %lld\n", cm_zlib.getMem());
    printf("in memory lz4 buffer size:       %lld\n", cm_lz4.getMem());
}

int main(int argc, char** argv) {
    int n = atoi(argv[1]),
        u = atoi(argv[2]),
        q = atoi(argv[3]);

    for (int u = 10000; u <= 50000; u += 10000) {
        foo(100000, u, 0.05);
    }
    // for (int i = 0; i < 10; i++) {
    //     int n = 100000 + i * 100000;
    // }

    // for (int i = 0; i < 10; i++) {
    // }

    // Print sizes of normal count-min on ssd for different n and delta.
    // 
}
