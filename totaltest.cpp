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

int main(int argc, char** argv) {
    int n = atoi(argv[1]),
        u = atoi(argv[2]),
        q = atoi(argv[3]);
    double epsilon = M_E / (10 * n),
           delta = 1 / pow(M_E, 3),
           epsilon_u = M_E / (10 * u);
    clock_t start, finish;

    // Test uniform distribution

    // naked count-min on disk
    CountMin cm_ssd(
        epsilon, delta, 1337, Uncompressed, u, "naked.cm",
        "naked count-min on disk"
    );

    // count-min on disk with a count-min buffer
    CountMin cm_ssd_cmbuf(
        epsilon, delta, 1337, Uncompressed, u, "cmbuf.cm",
        "count-min on disk for count-min buffer"
    );
    CountMin cm_cmbuf(
        epsilon_u, delta, 1337, Uncompressed, u, nullptr,
        "in memory count-min buffer"
    );

    // count-min on disk with a hashtable buffer
    CountMin cm_ssd_ht(
        epsilon, delta, 1337, Uncompressed, u, "hashtable.cm",
        "count-min on disk for hashtable buffer"
    );
    CountMin cm_hashtable(
        epsilon, delta, 1337, HashTable, u, nullptr,
        "in memory hashtable buffer"
    );

    // count-min on disk with a tree buffer
    CountMin cm_ssd_tree(
        epsilon, delta, 1337, Uncompressed, u, "tree.cm",
        "count-min on disk for tree buffer"
    );
    CountMin cm_tree(
        epsilon, delta, 1337, Tree, u, nullptr,
        "in memory tree buffer"
    );

    // count-min on disk with a raw log buffer
    CountMin cm_ssd_rawlog(
        epsilon, delta, 1337, Uncompressed, u, "rawlog.cm",
        "count-min on disk for raw log buffer"
    );
    CountMin cm_rawlog(
        0, 0, 1337, RawLog, u, nullptr,
        "in memory raw log buffer"
    );

    // update time
    auto uniform_kvs = generate_uniform_kvs(u);
    auto zipf_kvs = generate_zipf_kvs(u);

    start = clock();
    update_cm(cm_ssd, uniform_kvs);
    finish = clock();
    printf("update time for naked count-min on disk:    %lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_cmbuf, uniform_kvs);
    finish = clock();
    printf("update time for in memory count-min buffer: %lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_hashtable, uniform_kvs);
    finish = clock();
    printf("update time for in memory hashtable buffer: %lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_tree, uniform_kvs);
    finish = clock();
    printf("update time for in memory tree buffer:      %lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_rawlog, uniform_kvs);
    finish = clock();
    printf("update time for raw log buffer:             %lfms\n",
           clock_to_ms(finish - start));

    printf("naked count-min on disk size:    %lld\n", cm_ssd.getMem());
    printf("in memory count-min buffer size: %lld\n", cm_cmbuf.getMem());
    printf("in memory hashtable buffer size: %lld\n", cm_hashtable.getMem());
    printf("in memory tree buffer size:      %lld\n", cm_tree.getMem());
    printf("in memory raw log buffer size:   %lld\n", cm_rawlog.getMem());

    // point query time
    auto uniform_keys = generate_uniform_keys(q);
    int tmp;

    start = clock();
    tmp = pq_full_cm(cm_ssd, uniform_keys);
    finish = clock();
    printf("point query time for naked count-min on disk:    %lfms\n",
           clock_to_ms(finish - start));
    printf("### ignore this value: %d\n", tmp);

    start = clock();
    tmp = pq_full_cm(cm_cmbuf, uniform_keys);
    finish = clock();
    printf("point query time for in memory count-min buffer: %lfms\n",
           clock_to_ms(finish - start));
    printf("### ignore this value: %d\n", tmp);

    start = clock();
    tmp = pq_full_cm(cm_hashtable, uniform_keys);
    finish = clock();
    printf("point query time for in memory hashtable buffer: %lfms\n",
           clock_to_ms(finish - start));
    printf("### ignore this value: %d\n", tmp);

    start = clock();
    tmp = pq_full_cm(cm_tree, uniform_keys);
    finish = clock();
    printf("point query time for in memory tree buffer:      %lfms\n",
           clock_to_ms(finish - start));
    printf("### ignore this value: %d\n", tmp);

    start = clock();
    tmp = pq_full_cm(cm_rawlog, uniform_keys);
    finish = clock();
    printf("point query time for in memory raw log buffer:   %lfms\n",
           clock_to_ms(finish - start));
    printf("### ignore this value: %d\n", tmp);

    // inner product query time
    // TODO

    // Test Zipf distribution
    // TODO
}
