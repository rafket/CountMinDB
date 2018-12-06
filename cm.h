#ifndef CM_H
#define CM_H
#define CHUNKSIZE 1024

#include <bits/stdc++.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include "MurmurHash3.h"
#include "zlib.h"

using namespace std;

// decompress function from zlib (which uses adaptive Huffman coding apparently)
int* zlib_decompress(char* chunk, size_t compressed_size, size_t uncompressed_size) {
    int* output = new int[uncompressed_size]();

    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = (uInt) compressed_size;
    infstream.next_in = (Bytef*) chunk;
    infstream.avail_out = (uInt) (uncompressed_size);
    infstream.next_out = (Bytef*) output;
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream); 

    return output;
}

// compress function from zlib (which uses adaptive Huffman coding apparently)
size_t zlib_compress(int* count_array, size_t uncompressed_size, char** array) {
    // start with 0.01 * the size for now
    size_t compressed_size = max(1.0, uncompressed_size * 0.01);
    char* chunk = NULL;
    z_stream defstream;

    int r = Z_BUF_ERROR;
    // keep looping until we allocate the right amount for the compressed string
    while (r != Z_STREAM_END) {
        compressed_size *= 2;
        if (chunk) delete[] chunk;
        chunk = new char[compressed_size]();
        defstream.zalloc = Z_NULL;
        defstream.zfree = Z_NULL;
        defstream.opaque = Z_NULL;
        defstream.avail_in = (uInt) uncompressed_size;
        defstream.next_in = (Bytef*) count_array;
        defstream.avail_out = (uInt) compressed_size;
        defstream.next_out = (Bytef*) chunk;
        // can choose different parameter settings
       // deflateInit(&defstream, Z_BEST_COMPRESSION);
        deflateInit(&defstream, Z_BEST_SPEED);
        // deflateInit(&defstream, Z_NO_COMPRESSION);
        r = deflate(&defstream, Z_FINISH);
        deflateEnd(&defstream);
    }
    // printf("compressed size: %d, uncompressed size: %d\n", compressed_size, uncompressed_size);
    *array = chunk;
    return compressed_size;
}


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

enum storage_type {
    Uncompressed,
    HashTable,
    Tree,
    BufferedVersion,
    ChunksZlib
};

class CountMin {
public:
    CountMin(double eps, double delta, uint64_t seed, storage_type type, uint64_t num_updates, const char* filename, string name);
    ~CountMin();
    void update(uint64_t i, int c);                         // increase a_i by c
    int pointQuery(uint64_t i) const;                       // return a_i
    int innerProductQuery(const CountMin &other) const;     // return a * b (inner product)
    void mergeCMs(const CountMin &other);
    void mergeRawLog(const vector<pair<uint64_t, uint64_t>>& other, size_t u);
    int innerProductQueryRawLog(const vector<pair<uint64_t, uint64_t>>& other, size_t u);


    uint64_t getWidth() const {
        return w;
    }

    uint64_t getDepth() const {
        return d;
    }

    const int* getCounts() const {
        return (const int*)flatcounts;
    }

    const vector<Hashtable> &getHashTableCounts() const {
        return hash_table_counts;
    }

    const  vector<map<uint64_t,int>> &getTreeCounts() const {
        return tree_counts;
    }

    char*** getChunkZlibCompressions() const {
        return chunks_zlib;
    }

    size_t** getCompressedSizes() const {
        return compressed_sizes;
    }

    size_t getNumChunks() const {
        return num_chunks;
    }


    bool storageType() const {
        return type;
    }

    uint64_t getMem() {
        if (type == Uncompressed) {
            return w * d * sizeof(int)  /* size of the CM table */
                + d * sizeof(uint32_t)  /* size of seed for hash function */;
        }
        if (type == HashTable) {
            uint64_t out = 0;
            for (size_t i = 0; i < d; ++i) {
                out += hash_table_counts[i].getMem();  // size of each hashtable
            }
            out += d * sizeof(uint32_t);
            return out;
        }
        if (type == ChunksZlib) {
            uint64_t out = 0;
            for (size_t i = 0; i < d; i++) {
                for (size_t j = 0; j < num_chunks; j++) {
                    out += compressed_sizes[i][j];
                }
            }
            return out;
        }
        fprintf(stderr, "NOT IMPLEMENTED\n");
        return -1;   
    }

private:
    vector<uint32_t> hash_seed;
    uint64_t w, d;
    storage_type type;
    char* filename;
    int fd;
    string name;

    // storage formats
    int* flatcounts; // uncompressed
    vector<Hashtable> hash_table_counts; //hash table
    vector<map<uint64_t,int>> tree_counts; // tree
    char*** chunks_zlib; // zlib chunk
    size_t num_chunks; // number of chunks in each row
    size_t** compressed_sizes; // sizes of compressed chunks


    uint64_t hash(uint64_t key, uint32_t seed) const {
        uint64_t out[2];
        MurmurHash3_x64_128(&key, sizeof(uint64_t), seed, out);
        return out[1];
    }
};

CountMin::CountMin(double eps, double delta, uint64_t seed, storage_type type = Uncompressed, uint64_t num_updates = 0, const char* filename = NULL, string name = "noname") {
    this->type = type;
    this->name = name;
    w = ceil(M_E / eps);
    d = ceil(log(1 / delta));
    printf("%s : w: %" PRIu64 ", d: %" PRIu64 "\n", name.c_str(), w, d);
    if (type == Uncompressed) {
        if (filename == NULL) {
                this->filename = NULL;
                flatcounts = new int[d*w]();
        } else {
            // TODO: implement this
            this->filename = strdup(filename);
            fd = open(filename, O_RDWR | O_CREAT, 0666);
            ftruncate(fd, d * w * sizeof(int));
            flatcounts = (int*)mmap(0, d * w * sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        }
    }
    else if (type == HashTable) {
        assert(filename == NULL);
        hash_table_counts.resize(d, Hashtable(max(num_updates, (uint64_t) 100), 1337));
    }
    else if (type == Tree) {
        assert(filename == NULL);
        tree_counts.assign(d, map<uint64_t,int>());
    }
    else if (type == ChunksZlib) {
        assert(filename == NULL);
        num_chunks = (w / CHUNKSIZE);
        if (w % CHUNKSIZE != 0) {
            num_chunks++;
        }
        // printf("num chunks:%d\n", num_chunks);
        chunks_zlib = new char**[d]();
        compressed_sizes = new size_t*[d]();
        for (size_t i = 0; i < d; i++) {
            chunks_zlib[i] = new char*[num_chunks]();
            compressed_sizes[i] = new size_t[num_chunks](); 
        }
    }
    else {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    hash_seed.resize(d);
    mt19937_64 mt(seed);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    for (size_t i = 0; i < d; i++) {
        hash_seed[i] = dist(mt);        // generate random seeds for hash function
    }
}

CountMin::~CountMin() {
    if (filename != NULL) {
        munmap(flatcounts, w*d*sizeof(int));
        close(fd);
        unlink(filename); // Files should not be persistent in experiments
    } else if (type == Uncompressed) {
        delete[] flatcounts;
    }
}

void CountMin::update(uint64_t i, int c) {
    if (type == Uncompressed) {
        for (size_t j = 0; j < d; j++) {
            flatcounts[j*w + hash(i, hash_seed[j]) % w] += c;
        }
        return;
    }
    if (type == HashTable) {
        for (size_t j = 0; j < d; ++j) {
            hash_table_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
        return;
    }
    if (type == Tree) {
        for (size_t j = 0; j < d; ++j) {
            tree_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
        return;
    }
    if (type == ChunksZlib) {
        for (size_t j = 0; j < d; j++) {
            size_t index =  hash(i, hash_seed[j]) % w; // find out where the hash lives
            size_t chunk_block = (index / CHUNKSIZE); // find out which chunk block it is in
            size_t chunk_index = index % CHUNKSIZE; // find out which chunk index it is in the block
            size_t chunk_size = CHUNKSIZE; // size of uncompressed block
            // if it's the last block it may be smaller
            if (chunk_block == num_chunks - 1) {
                chunk_size = w % CHUNKSIZE;
            }
            int* inflated; 
            if (chunks_zlib[j][chunk_block]) {
                inflated = zlib_decompress(chunks_zlib[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int));
                delete[] chunks_zlib[j][chunk_block];
            }
            else {
                inflated = new int[CHUNKSIZE]();
            }
            inflated[chunk_index] += c;
            compressed_sizes[j][chunk_block] = zlib_compress(inflated, chunk_size * sizeof(int), &chunks_zlib[j][chunk_block]);
            delete[] inflated;
        }
        return;
    }
    fprintf(stderr, "NOT IMPLEMENTED\n");
    exit(0);
}

int CountMin::pointQuery(uint64_t i) const {
    int res = 0;
    if (type == Uncompressed) {
        res = flatcounts[hash(i, hash_seed[0]) % w];
        for (size_t j = 1; j < d; j++) {
            res = min(res, flatcounts[j*w + hash(i, hash_seed[j]) % w]);
        }
        return res;
    }
    if (type == HashTable) {
        auto it = hash_table_counts[0].find(hash(i, hash_seed[0]) % w);
        res = max(0, it);
        for (size_t j = 1; j < d; j++) {
            auto it = hash_table_counts[j].find(hash(i, hash_seed[j]) % w);
            res = min(res, max(0, it));
        }
        return res;
    }
    if (type == Tree) {
        auto it = tree_counts[0].find(hash(i, hash_seed[0]) % w);
        if (it == tree_counts[0].end()) {
            res = 0;
        }
        else {
            res = it->second;
        }
        for (size_t j = 1; j < d; j++) {
            auto it = tree_counts[j].find(hash(i, hash_seed[j]) % w);
            if (it == tree_counts[j].end()) {
                res = min(res, 0);
            }
            else {
                res = min(res, it->second);
            }
        }
        return res;
    }
    if (type == ChunksZlib) {
        for (size_t j = 0; j < d; j++) {
            size_t index =  hash(i, hash_seed[j]) % w; // find out where the hash lives
            size_t chunk_block = (index / CHUNKSIZE); // find out which chunk block it is in
            size_t chunk_index = index % CHUNKSIZE; // find out which chunk index it is in the block
            size_t chunk_size = CHUNKSIZE; // size of uncompressed block
            if (!chunks_zlib[j][chunk_block]) return 0;
            // if it's the last block it may be smaller
            if (chunk_block == num_chunks - 1) {
                chunk_size = w % CHUNKSIZE;
            }
            int* inflated; 
            inflated = zlib_decompress(chunks_zlib[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int));
            if (j == 0) res = inflated[chunk_index];
            res = min(res, inflated[chunk_index]);
            delete[] inflated;
        }
        return res;
    }
    fprintf(stderr, "NOT IMPLEMENTED\n");
    return -1;
}

int CountMin::innerProductQuery(const CountMin &other) const {
    int* totals = new int[d]();
    if (type == Uncompressed) {
        if (other.type == Uncompressed) {
            const int* otherCounts = other.getCounts();
            for (size_t i = 0; i < w; i++) {
                for (size_t j = 0; j < d; j++) {
                    totals[j] += flatcounts[j*w + i] * otherCounts[j*w + i];
                }
            }
        }
        if (other.type == HashTable) {
            const auto &otherSparseCounts = other.getHashTableCounts();
            for (size_t j = 0; j < d; j++) {
                for (const auto &it: otherSparseCounts[j].arr) {
                    if (it.first == (uint64_t)-1) {
                        continue;
                    }
                    totals[j] += flatcounts[j*w + it.first] * it.second;
                }
            }
        }
        if (other.type == Tree) {
            const auto &otherSparseCounts = other.getTreeCounts();
            for (size_t j = 0; j < d; j++) {
                for (const auto &it: otherSparseCounts[j]) {
                    totals[j] += flatcounts[j*w + it.first] * it.second;  
                }
            }
        }
        if (other.type == ChunksZlib) {
            char*** chunks = other.getChunkZlibCompressions();
            size_t** sizes = other.getCompressedSizes();
            size_t number_chunks = other.getNumChunks();
            size_t chunk_size = CHUNKSIZE;
            for (size_t j = 0; j < number_chunks; j++) {
                for (size_t i = 0; i < d; i++) {
                    chunk_size = CHUNKSIZE;
                    if (chunks[i][j]) {
                        if (j == number_chunks - 1) {
                            chunk_size = w % CHUNKSIZE;
                        }
                        int* inflated; 
                        inflated = zlib_decompress(chunks[i][j], sizes[i][j], chunk_size * sizeof(int));
                        for (size_t k = 0; k < chunk_size; k++) {
                            totals[j] += flatcounts[i * w + CHUNKSIZE * j + k] * inflated[k];
                        }
                        delete[] inflated;
                    }
                }
            }
        }
    }
    else if (type == HashTable && other.type == HashTable) {
        const auto &otherSparseCounts = other.getHashTableCounts();
        for (size_t j = 0; j < d; j++) {
            for (const auto it: otherSparseCounts[j].arr) {
                if (it.first == (uint64_t)-1) {
                    continue;
                }
                auto val = hash_table_counts[j].find(it.first);
                totals[j] += max(0, val) * it.second;
            }
        }
    }
    else if (type == Tree && other.type == Tree) {
        const auto &otherSparseCounts = other.getTreeCounts();
            for (size_t j = 0; j < d; j++) {
                for (const auto it: otherSparseCounts[j]) {
                    auto val = tree_counts[j].find(it.first);
                    if (val != tree_counts[j].end()) {
                        totals[j] += val->second * it.second;
                    }
                }
            }
    }
    else if (type == ChunksZlib && other.type == ChunksZlib) {
        char*** chunks = other.getChunkZlibCompressions();
        size_t** sizes = other.getCompressedSizes();
        size_t chunk_size = CHUNKSIZE;
        for (size_t j = 0; j < num_chunks; j++) {
            for (size_t i = 0; i < d; i++) {
                chunk_size = CHUNKSIZE;
                if (chunks_zlib[i][j] && chunks[i][j]) {
                    if (j == num_chunks - 1) {
                        chunk_size = w % CHUNKSIZE;
                    }
                    int* inflated1; 
                    int* inflated2;
                    inflated1 = zlib_decompress(chunks[i][j], sizes[i][j], chunk_size * sizeof(int));
                    inflated2 = zlib_decompress(chunks_zlib[i][j], compressed_sizes[i][j], chunk_size * sizeof(int));
                    for (size_t k = 0; k < chunk_size; k++) {
                        totals[j] += inflated1[k] * inflated2[k];
                    }
                    delete[] inflated1;
                    delete[] inflated2;
                }
            }
        }
    }
    else {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    int res = totals[0];
    for (size_t j=0; j<d; j++) {
        res = min(res, totals[j]);
    }
    delete[] totals;
    return res;
}

void CountMin::mergeCMs(const CountMin& other) {
    if (type != Uncompressed) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    if (other.type == Uncompressed) {
        const int* otherCounts = other.getCounts();
        for (size_t i = 0; i < d; ++i) {
            for (size_t j = 0; j < w; ++j) {
                flatcounts[i*w + j] += otherCounts[i*w + j];
            }
        }
        return;
    }
    if (other.type == HashTable) {
        const auto &otherSparseCounts = other.getHashTableCounts();
        for (size_t i = 0; i < d; ++i) {
            for (const auto &it: otherSparseCounts[i].arr) {
                if (it.first == (uint64_t)-1) {
                    continue;
                }
                flatcounts[i*w + it.first] += it.second;
            }
        }
        return;
    }
    if (other.type == Tree) {
        const auto &otherSparseCounts = other.getTreeCounts();
        for (size_t i = 0; i < d; i++) {
            for (auto &it: otherSparseCounts[i]) {
                flatcounts[i*w + it.first] += it.second;
            }
        }
        return;
    }
    if (other.type == ChunksZlib) {
        char*** chunks = other.getChunkZlibCompressions();
        size_t** sizes = other.getCompressedSizes();
        size_t number_chunks = other.getNumChunks();
        size_t chunk_size = CHUNKSIZE;
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < number_chunks; j++) {
                chunk_size = CHUNKSIZE;
                if (j == number_chunks - 1) {
                    chunk_size = w % CHUNKSIZE;
                }
                if (chunks[i][j]) {
                    int* inflated; 
                    inflated = zlib_decompress(chunks[i][j], sizes[i][j], chunk_size * sizeof(int));
                    for (size_t k = 0; k < chunk_size; k++) {
                        flatcounts[i * w + CHUNKSIZE * j + k] += inflated[k];
                    }
                    delete[] inflated;
                }
            }
        }
        return;
    }

    fprintf(stderr, "NOT IMPLEMENTED\n");
    exit(0);
}

// merge the raw log into the count-min directly with a for loop
void CountMin::mergeRawLog(const vector<pair<uint64_t, uint64_t>>& other, size_t u) {
    if (type != Uncompressed) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }

    for (size_t j = 0; j < u; j++) {
        update(other[j].first, other[j].second);
    }
}

// merge the raw log into the count-min directly with a for loop
int CountMin::innerProductQueryRawLog(const vector<pair<uint64_t, uint64_t>>& arr, size_t u) {
    int* totals = new int[d]();
    if (type == Uncompressed) {
        for(size_t i=0; i<u; ++i) {
            for (size_t j = 0; j < d; j++) {
                size_t index = hash(arr[i].first, hash_seed[j]) % w;
                totals[j] += flatcounts[j * w + index] * arr[i].second;
            }
        }
    }
    else if (type == HashTable) {
        for(size_t i=0; i<u; ++i) {
            for (size_t j = 0; j < d; j++) {
                size_t index = hash(arr[i].first, hash_seed[j]) % w;
                auto it = hash_table_counts[j].find(index);
                totals[j] += max(0, it) * arr[i].second;
            }
        }
    }
    else if (type == Tree) {
        for(size_t i=0; i<u; ++i) {
            for (size_t j = 0; j < d; j++) {
                size_t index = hash(arr[i].first, hash_seed[j]) % w;
                auto it = tree_counts[j].find(index);
                if (it != tree_counts[j].end()) {
                    totals[j] += it->second * arr[i].second;
                }
            }
        }
    }
    else if (type == ChunksZlib) {
        for(size_t i=0; i<u; ++i) {
            for (size_t j = 0; j < d; j++) {
                size_t index = hash(arr[i].first, hash_seed[j]) % w;
                size_t chunk_block = (index / CHUNKSIZE); // find out which chunk block it is in
                size_t chunk_index = index % CHUNKSIZE; // find out which chunk index it is in the block
                size_t chunk_size = CHUNKSIZE; // size of uncompressed block
                if (chunks_zlib[j][chunk_block]) {
                    if (chunk_block == num_chunks - 1) {
                        chunk_size = w % CHUNKSIZE;
                    }
                    int* inflated; 
                    inflated = zlib_decompress(chunks_zlib[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int));
                    totals += inflated[chunk_index] * arr[i].second;
                    delete[] inflated;
                }
            }
        }
    }
    else { 
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0); 
    }
    int res = totals[0];
    for (size_t j=0; j<d; j++) {
        res = min(res, totals[j]);
    }
    delete[] totals;
    return res;
}

// function for querying the raw log
uint64_t queryRawLog(const vector<pair<uint64_t, uint64_t>>& arr, uint64_t key, size_t u) {
    for(size_t j=0; j<u; ++j) {
        if(key==arr[j].first) {
            return arr[j].second;
        }
    }
    return 0;
}

// // function for inner product querying two raw logs
// uint64_t innerProductQueryRawLogs(const vector<pair<uint64_t, uint64_t>>& arr, const vector<pair<uint64_t, uint64_t>>& arr, uint64_t key, size_t u) {
//     // need to implement
// }

#endif
