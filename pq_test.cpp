#include <bits/stdc++.h>
#include "cm.h"

using namespace std;

struct TimeData {
    int cnt_upd_buf, cnt_upd_cm, cnt_upd_full;
    int cnt_pq_buf, cnt_pq_cm, cnt_pq_full;
    int cnt_ipq_buf_buf, cnt_ipq_buf_cm, cnt_ipq_cm_cm, cnt_ipq_full;

    clock_t t_upd_buf, t_upd_cm, t_upd_full;
    clock_t t_pq_buf, t_pq_cm, t_pq_full;
    clock_t t_ipq_buf_buf, t_ipq_buf_cm, t_ipq_cm_cm, t_ipq_full;
};

struct TimeCM {
    CountMin *buffer, *cm;
    TimeData td;
    bool record_full;

    TimeCM() : buffer(nullptr), cm(nullptr), record_full(true) {
        memset(&td, 0, sizeof(td));
    }

    ~TimeCM() {
        delete buffer;
        delete cm;
    }

    CountMin* bufferCreate(double eps, double delta, uint64_t seed,
                           storage_type type, uint64_t num_updates,
                           const char* name);

    CountMin* cmCreate(double eps, double delta, uint64_t seed,
                       const char* filename, const char* name);

    void update(uint64_t i, int c);
    void updateBuffer(uint64_t i, int c);
    void updateCM(uint64_t i, int c);
    int pointQuery(uint64_t i);
    int pointQueryBuffer(uint64_t i);
    int pointQueryCM(uint64_t i);
    int innerProductQuery(TimeCM &other);
    int innerProductQueryBufferBuffer(TimeCM &other);
    int innerProductQueryBufferCM(TimeCM &other);
    int innerProductQueryCMBuffer(TimeCM &other);
    int innerProductQueryCMCM(TimeCM &other);
};

CountMin* TimeCM::bufferCreate(
    double eps, double delta, uint64_t seed, storage_type type,
    uint64_t num_updates, const char* name = "noname"
) {
    buffer = new CountMin(eps, delta, seed, type, num_updates, nullptr, name);
    return buffer;
}

CountMin* TimeCM::cmCreate(
    double eps, double delta, uint64_t seed, const char* filename,
    const char* name = "noname"
) {
    cm = new CountMin(eps, delta, seed, Uncompressed, 0, filename, name);
    return cm;
}

void TimeCM::update(uint64_t i, int c) {
    clock_t start = 0;
    if (record_full) start = clock();
    if (buffer) {
        updateBuffer(i, c);
    } else if (cm) {
        updateCM(i, c);
    }
    if (record_full) {
        td.t_upd_full += clock() - start;
        td.cnt_upd_full++;
    }
}

void TimeCM::updateBuffer(uint64_t i, int c) {
    clock_t start = 0;
    if (buffer) {
        if (!record_full) start = clock();
        buffer->update(i, c);
        if (!record_full) {
            td.t_upd_buf += clock() - start;
            td.cnt_upd_buf++;
        }
    }
}

void TimeCM::updateCM(uint64_t i, int c) {
    clock_t start = 0;
    if (cm) {
        if (!record_full) start = clock();
        cm->update(i, c);
        auto finish = clock();
        if (!record_full) {
            td.t_upd_cm += clock() - start;
            td.cnt_upd_cm++;
        }
    }
}

int TimeCM::pointQuery(uint64_t i) {
    clock_t start = 0;
    if (record_full) start = clock();
    int res = pointQueryCM(i) + pointQueryBuffer(i);
    if (record_full) {
        td.t_pq_full += clock() - start;
        td.cnt_pq_full++;
    }
    return res;
}

int TimeCM::pointQueryBuffer(uint64_t i) {
    clock_t start = 0;
    if (buffer) {
        if (!record_full) start = clock();
        int res = buffer->pointQuery(i);
        if (!record_full) {
            td.t_pq_buf += clock() - start;
            td.cnt_pq_buf++;
        }
        return res;
    }
    return 0;
}

int TimeCM::pointQueryCM(uint64_t i) {
    clock_t start = 0;
    if (cm) {
        if (!record_full) start = clock();
        int res = cm->pointQuery(i);
        if (!record_full) {
            td.t_pq_cm += clock() - start;
            td.cnt_pq_cm++;
        }
        return res;
    }
    return 0;
}

int TimeCM::innerProductQuery(TimeCM &other) {
    return innerProductQueryBufferBuffer(other)
        + innerProductQueryBufferCM(other)
        + innerProductQueryCMBuffer(other)
        + innerProductQueryCMCM(other);
}

int TimeCM::innerProductQueryBufferBuffer(TimeCM &other) {
    if (buffer == nullptr || other.buffer == nullptr) {
        return 0;
    }
    return buffer->innerProductQuery(*other.buffer);
}

int TimeCM::innerProductQueryBufferCM(TimeCM &other) {
    if (buffer == nullptr || other.cm == nullptr) {
        return 0;
    }
    return buffer->innerProductQuery(*other.cm);
}

int TimeCM::innerProductQueryCMBuffer(TimeCM &other) {
    return other.innerProductQueryBufferCM(*this);
}

int TimeCM::innerProductQueryCMCM(TimeCM &other) {
    if (cm == nullptr || other.cm == nullptr) {
        return 0;
    }
    return cm->innerProductQuery(*other.cm);
}


int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Usage:\n    %s [number of elements in ssd] [number of updates in buffer] [number of queries]\n", argv[0]);
        exit(0);
    }

    mt19937_64 mt(1337);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    int n = atoi(argv[1]),
        u = atoi(argv[2]),
        q = atoi(argv[3]);
    double epsilon = M_E / (10 * n), delta = 1 /pow(M_E, 3),
           epsilon_u = M_E / (10 * u);

    TimeCM tcm_ssd, tcm_ram, tcm_hashtable;
    tcm_ssd.cmCreate(epsilon, delta, 1337, "file.cm", "count-min on SSD");
    tcm_ram.bufferCreate(epsilon_u, delta, 1337, Uncompressed, u, "count-min buffer in RAM");
    tcm_hashtable.bufferCreate(epsilon, delta, 1337, HashTable, u, "hashtable buffer in RAM");
    // CountMin cm_ssd(epsilon, delta, 1337, Uncompressed, u, "file.cm", "count-min on SSD");
    // CountMin cm_optimized(epsilon_u, delta, 1337, Uncompressed, u, nullptr, "optimized count-min in RAM");
    // CountMin cm_sparse(epsilon, delta, 1337, HashTable, u, nullptr, "sparse count-min in RAM");
    vector<pair<uint64_t, uint64_t>> arr(u);
    clock_t start, finish;

    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        tcm_ssd.update(key, value);
    }
    printf("building countMin in SSD took %lfms\n",
           (double)tcm_ssd.td.t_upd_full * 1000 / CLOCKS_PER_SEC);
    printf("SSD countMin takes up %lfMB\n",
           (double)tcm_ssd.cm->getMem() / 1024 / 1024);

    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        tcm_ram.update(key, value);
    }
    printf("building optimized countMin took %lfms\n",
           (double)tcm_ram.td.t_upd_full * 1000 / CLOCKS_PER_SEC);
    printf("optimized countMin takes up %lfMB\n",
           (double)tcm_ram.buffer->getMem() / 1024 / 1024);

    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        tcm_hashtable.update(key, value);
    }
    printf("building sparse countMin took %lfms\n",
           (double)tcm_hashtable.td.t_upd_full * 1000 / CLOCKS_PER_SEC);
    printf("sparse countMin takes up %lfMB\n",
           (double)tcm_hashtable.buffer->getMem() / 1024 / 1024);

    start=clock();
    for(int i=0; i<u; ++i) {
        uint64_t key = dist(mt);
        uint64_t value = 10;
        arr[i] = {key, value};
    }
    finish=clock();
    printf("building buffer took %lfms\n",
           (double)(finish - start) * 1000 / CLOCKS_PER_SEC);
    printf("buffer takes up %lfMB\n",
           (double)arr.size() * (sizeof(uint64_t) + sizeof(uint64_t))
           / 1024 / 1024);

    uint64_t sum1=0, sum2=0, sum3=0, sum4=0;

    for(int i=0; i<q; ++i) {
        sum1 += tcm_ssd.pointQuery(dist(mt));
    }
    printf("querying ssd countMin took %lfms\n",
           (double)tcm_ssd.td.t_pq_full * 1000 / CLOCKS_PER_SEC);

    for(int i=0; i<q; ++i) {
        sum4 += tcm_ram.pointQuery(dist(mt));
    }
    printf("querying optimized countMin took %lfms\n",
           (double)tcm_ram.td.t_pq_full * 1000 / CLOCKS_PER_SEC);

    for(int i=0; i<q; ++i) {
        sum2 += tcm_hashtable.pointQuery(dist(mt));
    }
    printf("querying sparse countMin took %lfms\n",
           (double)tcm_hashtable.td.t_pq_full * 1000 / CLOCKS_PER_SEC);

    start=clock();
    for(int i=0; i<q; ++i) {
        sum3 += queryRawLog(arr, dist(mt), u);
    }
    finish=clock();
    printf("pay no attention to these numbers: %lu %lu %lu %lu\n",
           sum1, sum2, sum3, sum4);
    printf("querying buffer took %lfms\n",
           (double)(finish - start) * 1000 / CLOCKS_PER_SEC);
}
