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

void populate_cm(CountMin &cm, int n) {
    auto uniform_kvs = generate_uniform_kvs(n);
    update_cm(cm, uniform_kvs);
}

void experiment_update_time(long long n, int u, double delta, double c, bool need_merge = true) {
    double epsilon = M_E / (c * n),
           epsilon_u = M_E / (c * u);
    clock_t start, finish;

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

    printf("updating with u elements...\n");
    // update time
    auto uniform_kvs = generate_uniform_kvs(u);
    // auto zipf_kvs = generate_zipf_kvs(u);

    start = clock();
    update_cm(cm_ssd, uniform_kvs);
    finish = clock();
    printf("update time for naked count-min on disk:\t%lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_cmbuf, uniform_kvs);
    finish = clock();
    printf("update time for in memory count-min buffer:\t%lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_hashtable, uniform_kvs);
    finish = clock();
    printf("update time for in memory hashtable buffer:\t%lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_tree, uniform_kvs);
    finish = clock();
    printf("update time for in memory tree buffer:\t%lfms\n",
           clock_to_ms(finish - start));

    start = clock();
    update_cm(cm_rawlog, uniform_kvs);
    finish = clock();
    printf("update time for in memory raw log buffer:\t%lfms\n",
           clock_to_ms(finish - start));

    printf("done\n");

    // time to merge
    if (need_merge) {
        printf("merging with ssd countmins...\n");

        start = clock();
        cm_ssd_ht.mergeCMs(cm_hashtable);
        finish = clock();
        printf("merge time for hashtable:\t%lfms\n", clock_to_ms(finish - start));

        start = clock();
        cm_ssd_tree.mergeCMs(cm_tree);
        finish = clock();
        printf("merge time for tree:\t%lfms\n", clock_to_ms(finish - start));

        start = clock();
        cm_ssd_rawlog.mergeCMs(cm_rawlog);
        finish = clock();
        printf("merge time for raw log:\t%lfms\n", clock_to_ms(finish - start));

        printf("done\n");
    }
}

void full_experiment_update_time() {
    long long ns[] = {
        500000000
    };
    int us[] = {
        5000,
        // 10000,
        // 15000,
        // 20000,
        // 25000,
        // 30000,
        // 35000,
        // 40000,
        // 45000,
        // 50000,
        // 75000,
        // 100000
    };
    double deltas[] = {
        0.0001,
        0.001,
        0.01,
        0.05,
        0.08,
        0.1
    };
    double cs[] = {
        0.01, 0.03, 0.05, 0.07, 0.1
        // 10.0,
        // 5.0,
        // 2.5,
        // 0.7,
        // 0.15,
        // 0.05,
        // 0.01
    };
    constexpr long long N = 500000000;
    constexpr int U = 40000;
    constexpr double DELTA = 0.05;
    constexpr double C = 2.5;
    // for (auto u : us) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, u, DELTA, C);
    //     experiment_update_time(N, u, DELTA, C);
    //     fprintf(stderr, "yay\n");
    // }
    // for (auto delta : deltas) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, U, delta, C);
    //     experiment_update_time(N, U, delta, C);
    //     fprintf(stderr, "yay\n");
    // }
    for (auto c : cs) {
        printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
               N, U, DELTA, c);
        experiment_update_time(N, U, DELTA, c);
        fprintf(stderr, "yay\n");
    }
}

int query_cm(
    CountMin &cm, CountMin &cm_ssd,
    const vector<uint64_t> &keys, bool need_full
) {
    int res = 0;
    if (need_full) {
        for (auto key : keys) {
            res += cm.pointQuery(key) + cm_ssd.pointQuery(key);
        }
    } else {
        for (auto key : keys) {
            res += cm.pointQuery(key);
        }
    }
    return res;
}

void experiment_point_query_time(long long n, int u, double delta, double c, int q) {
    double epsilon = M_E / (c * n),
           epsilon_u = M_E / (c * u);
    clock_t start, finish;
    int sum;

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
    // auto zipf_kvs = generate_zipf_kvs(u);

    printf("updating cms...\n");
    update_cm(cm_ssd, uniform_kvs);
    update_cm(cm_cmbuf, uniform_kvs);
    update_cm(cm_hashtable, uniform_kvs);
    update_cm(cm_tree, uniform_kvs);
    update_cm(cm_rawlog, uniform_kvs);
    update_cm(cm_zlib, uniform_kvs);
    update_cm(cm_lz4, uniform_kvs);
    update_cm(cm_bufhash, uniform_kvs);
    update_cm(cm_buftree, uniform_kvs);
    update_cm(cm_bufraw, uniform_kvs);
    printf("done\n");

    auto uniform_keys = generate_uniform_keys(q);
    int totalsum = 0;

    // do this to cache
    sum = query_cm(cm_ssd, cm_ssd, uniform_keys, false);
    totalsum += sum;

    start = clock();
    sum = query_cm(cm_ssd, cm_ssd, uniform_keys, false);
    finish = clock();
    totalsum += sum;
    printf("query throughput for naked count-min on disk:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_cmbuf, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory count-min buffer:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_hashtable, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory hashtable buffer:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_tree, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory tree buffer:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_rawlog, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory raw log buffe:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_zlib, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory zlib buffer:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    start = clock();
    sum = query_cm(cm_lz4, cm_ssd, uniform_keys, true);
    finish = clock();
    totalsum += sum;
    printf("query throughput for in memory lz4 buffer:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    // do this to cache
    sum = query_cm(cm_bufhash, cm_bufhash, uniform_keys, false);
    totalsum += sum;

    start = clock();
    sum = query_cm(cm_bufhash, cm_bufhash, uniform_keys, false);
    finish = clock();
    totalsum += sum;
    printf("query throughput for buffer hash countmin:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    // do this to cache
    sum = query_cm(cm_buftree, cm_buftree, uniform_keys, false);
    totalsum += sum;

    start = clock();
    sum = query_cm(cm_buftree, cm_buftree, uniform_keys, false);
    finish = clock();
    totalsum += sum;
    printf("query throughput for buffer tree countmin:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    // do this to cache
    sum = query_cm(cm_bufraw, cm_bufraw, uniform_keys, false);
    totalsum += sum;

    start = clock();
    sum = query_cm(cm_bufraw, cm_bufraw, uniform_keys, false);
    finish = clock();
    totalsum += sum;
    printf("query throughput for buffer raw countmin:\t%lf q/sec\n",
           q / (clock_to_ms(finish - start) / 1000));

    printf("### %d\n", totalsum);
}

void full_experiment_point_query_time() {
    long long ns[] = {
        500000000
    };
    int us[] = {
        5000,
        20000,
        35000,
        50000,
    };
    double deltas[] = {
       0.0001  , 0.002595, 0.00509 , 0.007585, 0.01008 , 0.012575,
       0.01507 , 0.017565, 0.02006 , 0.022555, 0.02505 , 0.027545,
       0.03004 , 0.032535, 0.03503 , 0.037525, 0.04002 , 0.042515,
       0.04501 , 0.047505, 0.05
    };
    double cs[] = {
        10.0,
        5.0,
        2.5,
        0.7,
        0.15,
        0.05,
        0.01
    };
    constexpr long long N = 500000000;
    constexpr int U = 5000;
    constexpr double DELTA = 0.0001;
    constexpr double C = 0.25;
    constexpr int Q = 30000;
    // experiment_point_query_time(N, U, DELTA, C, Q);
    // for (auto u : us) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, u, DELTA, C);
    //     experiment_update_time(N, u, DELTA, C);
    //     fprintf(stderr, "yay\n");
    // }
    for (auto delta : deltas) {
        printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
               N, U, delta, C);
        experiment_point_query_time(N, U, delta, C, Q);
        fprintf(stderr, "yay\n");
    }
    // for (auto c : cs) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, U, DELTA, c);
    //     experiment_update_time(N, U, DELTA, c);
    //     fprintf(stderr, "yay\n");
    // }
}

void experiment_inner_product_query_time(long long n, int u, double delta, double c) {
    double epsilon = M_E / (c * n * n),
           epsilon_u = M_E / (c * u * u);
    clock_t start, finish;
    int sum;

    // naked count-min on disk
    CountMin cm_ssd(
        epsilon, delta, 1337, Uncompressed, u, "naked.cm",
        "naked count-min on disk"
    );
    CountMin cm_ssd1(
        epsilon, delta, 1337, Uncompressed, u, "naked1.cm",
        "naked count-min on disk"
    );

    // count-min on disk with a count-min buffer
    CountMin cm_cmbuf(
        epsilon_u, delta, 1337, Uncompressed, u, nullptr,
        "in memory count-min buffer"
    );
    CountMin cm_cmbuf1(
        epsilon_u, delta, 1337, Uncompressed, u, nullptr,
        "in memory count-min buffer"
    );

    // count-min on disk with a hashtable buffer
    CountMin cm_hashtable(
        epsilon, delta, 1337, HashTable, u, nullptr,
        "in memory hashtable buffer"
    );
    CountMin cm_hashtable1(
        epsilon, delta, 1337, HashTable, u, nullptr,
        "in memory hashtable buffer"
    );

    // count-min on disk with a tree buffer
    CountMin cm_tree(
        epsilon, delta, 1337, Tree, u, nullptr,
        "in memory tree buffer"
    );
    CountMin cm_tree1(
        epsilon, delta, 1337, Tree, u, nullptr,
        "in memory tree buffer"
    );

    // count-min on disk with a raw log buffer
    CountMin cm_rawlog(
        0, 0, 1337, RawLog, u, nullptr,
        "in memory raw log buffer"
    );
    CountMin cm_rawlog1(
        0, 0, 1337, RawLog, u, nullptr,
        "in memory raw log buffer"
    );

    CountMin cm_zlib(
        epsilon, delta, 1337, ChunksZlib, u, nullptr,
        "in memory zlib buffer"
    );
    CountMin cm_zlib1(
        epsilon, delta, 1337, ChunksZlib, u, nullptr,
        "in memory zlib buffer"
    );

    CountMin cm_lz4(
        epsilon, delta, 1337, LZ4, u, nullptr,
        "in memory lz4 buffer"
    );
    CountMin cm_lz41(
        epsilon, delta, 1337, LZ4, u, nullptr,
        "in memory lz4 buffer"
    );

    CountMin cm_bufhash(
        epsilon, delta, 1337, BufferedHash, u, "bufhash.cm",
        "buffer hash countmin"
    );
    CountMin cm_bufhash1(
        epsilon, delta, 1337, BufferedHash, u, "bufhash1.cm",
        "buffer hash countmin"
    );

    CountMin cm_buftree(
        epsilon, delta, 1337, BufferedTree, u, "buftree.cm",
        "buffer tree countmin"
    );
    CountMin cm_buftree1(
        epsilon, delta, 1337, BufferedTree, u, "buftree1.cm",
        "buffer tree countmin"
    );

    CountMin cm_bufraw(
        epsilon, delta, 1337, BufferedRaw, u, "bufraw.cm",
        "buffer raw countmin"
    );
    CountMin cm_bufraw1(
        epsilon, delta, 1337, BufferedRaw, u, "bufraw1.cm",
        "buffer raw countmin"
    );

    auto uniform_kvs = generate_uniform_kvs(u);
    // auto zipf_kvs = generate_zipf_kvs(u);

    // printf("updating cms...\n");
    update_cm(cm_ssd, uniform_kvs);
    update_cm(cm_cmbuf, uniform_kvs);
    update_cm(cm_hashtable, uniform_kvs);
    update_cm(cm_tree, uniform_kvs);
    update_cm(cm_rawlog, uniform_kvs);
    update_cm(cm_zlib, uniform_kvs);
    update_cm(cm_lz4, uniform_kvs);
    update_cm(cm_bufhash, uniform_kvs);
    update_cm(cm_buftree, uniform_kvs);
    update_cm(cm_bufraw, uniform_kvs);
    update_cm(cm_ssd1, uniform_kvs);
    update_cm(cm_cmbuf1, uniform_kvs);
    update_cm(cm_hashtable1, uniform_kvs);
    update_cm(cm_tree1, uniform_kvs);
    update_cm(cm_rawlog1, uniform_kvs);
    update_cm(cm_zlib1, uniform_kvs);
    update_cm(cm_lz41, uniform_kvs);
    update_cm(cm_bufhash1, uniform_kvs);
    update_cm(cm_buftree1, uniform_kvs);
    update_cm(cm_bufraw1, uniform_kvs);
    // printf("done\n");

    int totalsum = 0;
    double ssd_ssd, ssd_ht, ssd_tree, ssd_rawlog, ht_rawlog, tree_rawlog, rawlog_rawlog, ht_ht, tree_tree, bufhash_x2, buftree_x2, bufraw_x2;
    // warm up the cache
    totalsum += cm_ssd.innerProductQuery(cm_ssd);

    start = clock();
    totalsum += cm_ssd.innerProductQuery(cm_ssd);
    finish = clock();
    // printf("ssd countmin x ssd countmin:\t%lfms\n",
    //        ssd_ssd = clock_to_ms(finish - start));
    ssd_ssd = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_ssd.innerProductQuery(cm_hashtable);
    finish = clock();
    // printf("in memory hashtable x ssd countmin:\t%lfms\n",
    //        ssd_ht = clock_to_ms(finish - start));
    ssd_ht = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_ssd.innerProductQuery(cm_tree);
    finish = clock();
    // printf("in memory tree x ssd countmin:\t%lfms\n",
    //        ssd_tree = clock_to_ms(finish - start));
    ssd_tree = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_ssd.innerProductQuery(cm_rawlog);
    finish = clock();
    // printf("in memory raw log x ssd countmin:\t%lfms\n",
    //        ssd_rawlog = clock_to_ms(finish - start));
    ssd_rawlog = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_hashtable.innerProductQuery(cm_rawlog);
    finish = clock();
    // printf("in memory raw log x in memory hashtable:\t%lfms\n",
        //    ht_rawlog = clock_to_ms(finish - start));
    ht_rawlog = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_tree.innerProductQuery(cm_rawlog);
    finish = clock();
    // printf("in memory raw log x in memory tree:\t%lfms\n",
        //    tree_rawlog = clock_to_ms(finish - start));
    tree_rawlog = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_rawlog.innerProductQuery(cm_rawlog1);
    finish = clock();
    // printf("in memory raw log x in memory raw log:\t%lfms\n",
        //    rawlog_rawlog = clock_to_ms(finish - start));
    rawlog_rawlog = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_hashtable.innerProductQuery(cm_hashtable1);
    finish = clock();
    // printf("in memory hashtable x in memory hashtable:\t%lfms\n",
        //    ht_ht = clock_to_ms(finish - start));
    ht_ht = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_tree.innerProductQuery(cm_tree1);
    finish = clock();
    // printf("in memory tree x in memory tree:\t%lfms\n",
        //    tree_tree = clock_to_ms(finish - start));
    tree_tree = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_bufhash.innerProductQuery(cm_bufhash1);
    finish = clock();
    // printf("buf hash x buf hash:\t%lfms\n",
        //    bufhash_x2 = clock_to_ms(finish - start));
    bufhash_x2 = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_buftree.innerProductQuery(cm_buftree1);
    finish = clock();
    // printf("buf tree x buf tree:\t%lfms\n",
        //    buftree_x2 = clock_to_ms(finish - start));
    buftree_x2 = clock_to_ms(finish - start);

    start = clock();
    totalsum += cm_bufraw.innerProductQuery(cm_bufraw1);
    finish = clock();
    // printf("buf raw x buf raw:\t%lfms\n",
        //    bufraw_x2 = clock_to_ms(finish - start));
    bufraw_x2 = clock_to_ms(finish - start);

    printf("### %d\n", totalsum);

    // printf("ignore everything above\n");

    printf("ssd x ssd:\t%lfms\n",
           ssd_ssd);
    printf("ssd x ht:\t%lfms\n",
           ssd_ssd + ssd_ht);
    printf("ssd x tree:\t%lfms\n",
           ssd_ssd + ssd_tree);
    printf("ssd x raw log:\t%lfms\n",
           ssd_ssd + ssd_rawlog);
    printf("raw log x ht:\t%lfms\n",
           ssd_ssd + ssd_ht + ssd_rawlog + ht_rawlog);
    printf("raw log x tree:\t%lfms\n",
           ssd_ssd + ssd_tree + ssd_rawlog + tree_rawlog);
    printf("raw log x raw log:\t%lfms\n",
           ssd_ssd + ssd_rawlog + ssd_rawlog + rawlog_rawlog);
    printf("ht x ht:\t%lfms\n",
           ssd_ssd + ssd_ht + ssd_ht + ht_ht);
    printf("tree x tree:\t%lfms\n",
           ssd_ssd + ssd_tree + ssd_tree + tree_tree);
    printf("bufht x bufht:\t%lfms\n", bufhash_x2);
    printf("buftree x buftree:\t%lfms\n", buftree_x2);
    printf("bufraw x bufraw:\t%lfms\n", bufraw_x2);
}

void full_experiment_inner_product_query_time() {
    long long ns[] = {
        500000000
    };
    int us[] = {
        5000,
        20000,
        35000,
        50000,
    };
    double deltas[] = {
       0.0001 , 0.00509, 0.01008, 0.01507, 0.02006, 0.02505, 0.03004,
       0.03503, 0.04002, 0.04501, 0.05
    };
    double cs[] = {
        10.0,
        5.0,
        2.5,
        0.7,
        0.15,
        0.05,
        0.01
    };
    constexpr long long N = 10000;
    constexpr int U = 5000;
    constexpr double DELTA = 0.001;
    constexpr double C = 0.25;
    // experiment_inner_product_query_time(N, U, DELTA, C);
    // for (auto u : us) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, u, DELTA, C);
    //     experiment_inner_product_query_time(N, u, DELTA, C);
    //     fprintf(stderr, "yay\n");
    // }
    for (auto delta : deltas) {
        printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
               N, U, delta, C);
        experiment_inner_product_query_time(N, U, delta, C);
        printf("\n");
        fprintf(stderr, "yay\n");
    }
    // for (auto c : cs) {
    //     printf("n = %lld\nu = %d\ndelta = %f\nc = %f\n",
    //            N, U, DELTA, c);
    //     experiment_inner_product_query_time(N, U, DELTA, c);
    //     fprintf(stderr, "yay\n");
    // }
}

int main(int argc, char** argv) {
    // full_experiment_update_time();
    // full_experiment_point_query_time();
    full_experiment_inner_product_query_time();
}
