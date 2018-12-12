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

void experiment_space(long long n, int u, double delta, double c) {
    double epsilon = M_E / (c * n),
           epsilon_u = M_E / (c * u);

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

    CountMin cm_bufhash(
        epsilon, delta, 1337, BufferedHash, u, "bufhash.cm",
        "buffer hash countmin"
    );

    CountMin cm_buftree(
        epsilon, delta, 1337, BufferedTree, u, "buftree.cm",
        "buffer tree countmin"
    );
    
    CountMin cm_bufraw(
        epsilon, delta, 1337, BufferedRaw, u, "bufraw.cm",
        "buffer raw countmin"
    );

    auto uniform_kvs = generate_uniform_kvs(u);

    update_cm(cm_hashtable, uniform_kvs);
    update_cm(cm_tree, uniform_kvs);
    update_cm(cm_rawlog, uniform_kvs);
    update_cm(cm_zlib, uniform_kvs);
    update_cm(cm_lz4, uniform_kvs);
    update_cm(cm_bufhash, uniform_kvs);
    update_cm(cm_buftree, uniform_kvs);
    update_cm(cm_bufraw, uniform_kvs);

    printf("naked count-min on disk size:\t%lld\n", cm_ssd.getMem());
    printf("in memory count-min buffer size:\t%lld\n", cm_cmbuf.getMem());
    printf("in memory hashtable buffer size:\t%lld\n", cm_hashtable.getMem());
    printf("in memory tree buffer size:\t%lld\n", cm_tree.getMem());
    printf("in memory raw log buffer size:\t%lld\n", cm_rawlog.getMem());
    printf("in memory zlib buffer size:\t%lld\n", cm_zlib.getMem());
    printf("in memory lz4 buffer size:\t%lld\n", cm_lz4.getMem());
    printf("buffer hash countmin size:\t%lld\n", cm_bufhash.getMem());
    printf("buffer tree countmin size:\t%lld\n", cm_buftree.getMem());
    printf("buffer raw countmin size:\t%lld\n", cm_bufraw.getMem());
}

void full_experiment_space() {
    long long ns[] = {
        500000000
    };
    int us[] = {
        5000,
        15000,
        25000,
        35000,
        45000
    };
    double deltas[] = {
        0.02,
        0.05,
        0.08,
        0.1
    };
    double cs[] = {
        10.0,
        5.0,
        2.5,
        0.7,
        // 0.15,
        // 0.05,
        // 0.01
    };
    constexpr long long N = 500000000;
    constexpr int U = 25000;
    constexpr double DELTA = 0.05;
    constexpr double C = 5.0;
    for (auto u : us) {
        printf("n = %lld\nu = %d\ndelta = %f\neps = %e\n",
               N, u, DELTA, M_E / (C * N));
        experiment_space(N, u, DELTA, C);
        printf("\n");
        fprintf(stderr, "yay\n");
    }
    for (auto delta : deltas) {
        printf("n = %lld\nu = %d\ndelta = %f\neps = %e\n",
               N, U, delta, M_E / (C * N));
        experiment_space(N, U, delta, C);
        printf("\n");
        fprintf(stderr, "yay\n");
    }
    for (auto c : cs) {
        printf("n = %lld\nu = %d\ndelta = %f\neps = %e\n",
               N, U, DELTA, M_E / (c * N));
        experiment_space(N, U, DELTA, c);
        printf("\n");
        fprintf(stderr, "yay\n");
    }
}

int main(int argc, char** argv) {
    full_experiment_space();
}
