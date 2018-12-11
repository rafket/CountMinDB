#ifndef CM_H
#define CM_H
#define CHUNKSIZE 4096
#define BLOCKSIZE 268435456 // 256 KB -- One SSD block
#define BUFFERSIZE 1024

#include <bits/stdc++.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include "lz4.h"
#include "MurmurHash3.h"
#include "zlib.h"

using namespace std;

enum storage_type {
    Uncompressed,
    HashTable,
    Tree,
    BufferedHash,
    BufferedTree,
    BufferedRaw,
    ChunksZlib,
    RawLog,
    LZ4
};

size_t lz4_compress(int* source, size_t source_size, char** dest) {
    size_t dest_size = 16;
    *dest = nullptr;
    do {
        dest_size *= 2;
        if(*dest != nullptr) {
            delete[] *dest;
        }
        *dest = new char[dest_size];
    } while( !LZ4_compress_fast((char*) source, *dest, source_size, dest_size, 10) );
    return dest_size;
}

int* lz4_decompress(char* source, size_t source_size, size_t dest_size) {
    int* dest = new int[dest_size];
    LZ4_decompress_safe(source, (char*) dest, source_size, dest_size);
    return dest;
}


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


size_t compression(int* count_array, size_t uncompressed_size, char** array, storage_type type) {
    if (type == ChunksZlib) {
        return zlib_compress(count_array, uncompressed_size, array);
    }
    if (type == LZ4) {
        return lz4_compress(count_array, uncompressed_size, array);
    }
    return 0;
}

int* decompression(char* chunk, size_t compressed_size, size_t uncompressed_size, storage_type type)  {
    if (type == ChunksZlib) {
        return zlib_decompress(chunk, compressed_size, uncompressed_size);
    }
    if (type == LZ4) {
        return lz4_decompress(chunk, compressed_size, uncompressed_size);
    }
    return nullptr;
}



struct Array {
    char* filename;
    int* data;
    size_t length;
    short int fd;

    Array() {
        data = NULL;
        filename = NULL;
        length = 0;
    }

    Array(size_t length) {
        this->length = length;
        filename = NULL;
        data = new int[length]();
    }

    Array(const char* filename, size_t length) {
        this->filename = strdup(filename);
        this->length = length;
        fd = open(this->filename, O_RDWR | O_CREAT, 0666);
        ftruncate(fd, length * sizeof(int));
        data = (int*)mmap(0, length * sizeof(int), PROT_READ | PROT_WRITE,
                MAP_SHARED, fd, 0);
    }

    /* This is so that flatcounts is not a pointer, because it can't be
     * initialized in the initializer list. It doesn't free existing data.
     */
    void reset(const char* filename, size_t length) {
        if(filename == NULL) {
            this->length = length;
            this->filename = NULL;
            data = new int[length]();
        }
        else {
            this->filename = strdup(filename);
            this->length = length;
            fd = open(this->filename, O_RDWR | O_CREAT, 0666);
            ftruncate(fd, length * sizeof(int));
            data = (int*)mmap(0, length * sizeof(int), PROT_READ | PROT_WRITE,
                    MAP_SHARED, fd, 0);
        }
    }

    void clear() {
        memset(data, 0, length*sizeof(int));
    }

    ~Array() {
        if(data == NULL) {
            return;
        }
        if(this->filename != NULL) {
            munmap(data, length * sizeof(int));
            close(fd);
            unlink(filename);
            free(filename);
        }
        else {
            delete[] data;
        }
    }

    int operator[](size_t idx) const {
        assert(idx<length);
        return data[idx];
    }

    int& operator[](size_t idx) {
        assert(idx<length);
        return data[idx];
    }
};

struct Hashtable {
    uint32_t seed;
    vector<pair<uint64_t, int>> arr;

    Hashtable(uint32_t n, uint32_t seed_) : seed(seed_) {
        arr.resize(n, make_pair((uint64_t)-1, 0));
    }

    size_t size() const {
        return arr.size();
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

    void clear() {
        fill(arr.begin(), arr.end(), make_pair((uint64_t)-1, 0));
    }
};


class CountMin {
public:
    CountMin(double eps, double delta, uint64_t seed, storage_type type, uint64_t num_updates, const char* filename, const char* name);
    ~CountMin();
    void update(uint64_t i, int c);                         // increase a_i by c
    int pointQuery(uint64_t i) const;                       // return a_i
    int innerProductQuery(const CountMin &other) const;     // return a * b (inner product)
    void mergeCMs(const CountMin &other);
    void mergeRawLog(const vector<pair<uint64_t, uint64_t>>& other, size_t u);
    int innerProductQueryRawLog(const vector<pair<uint64_t, uint64_t>>& other, size_t u);
    void clear();


    uint64_t getWidth() const {
        return w;
    }

    uint64_t getDepth() const {
        return d;
    }

    const Array& getCounts() const {
        return (const Array&)flatcounts;
    }

    const vector<Hashtable> &getHashTableCounts() const {
        return hash_table_counts;
    }

    const vector<map<uint64_t,int>> &getTreeCounts() const {
        return tree_counts;
    }

    const vector<pair<uint64_t, uint64_t>> &getRawLog() const {
        return arr;
    }

    const vector<CountMin>* getBchunks() const {
        return (const vector<CountMin>*)bchunks;
    }

    char*** getChunkZlibCompressions() const {
        return chunks_compression;
    }

    size_t** getCompressedSizes() const {
        return compressed_sizes;
    }

    size_t getNumChunks() const {
        return num_chunks;
    }

    size_t getNumUpdates() const {
        return num_updates;
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
        if (type == ChunksZlib || type == LZ4) {
            uint64_t out = 0;
            for (size_t i = 0; i < d; i++) {
                for (size_t j = 0; j < num_chunks; j++) {
                    out += compressed_sizes[i][j];
                }
            }
            return out;
        }
        if (type == Tree) {
            uint64_t out = 0;
            for (size_t i = 0; i < d; ++i) {
                out += tree_counts[i].size() * (16 + sizeof(pair<uint64_t,int>));  // size of each hashtable
            }
            return out;
        }
        if (type == RawLog) {
            uint64_t out = arr.size() * sizeof(arr[0]);
            return out;
        }
        if (type == BufferedHash || type == BufferedTree || type == BufferedRaw) {
            uint64_t out = 0;
            for(auto &i: bchunks[1]) {
                out += i.getMem();
            }
            return out;
        }
        fprintf(stderr, "NOT IMPLEMENTED\n");
        return -1;
    }

    size_t getNumContents() const {
        return num_contents;
    }

private:
    vector<uint32_t> hash_seed;
    uint64_t w, d;
    storage_type type;
    char* filename;
    int fd;
//    string name;
    uint32_t bchunk_hash_seed; // Buffered CountMin
    size_t num_contents;

    // storage formats
    Array flatcounts; // uncompressed
    vector<Hashtable> hash_table_counts; // hash table
    vector<map<uint64_t,int>> tree_counts; // tree
    vector<CountMin> bchunks[2]; // Buffered CountMin
    vector<pair<uint64_t, uint64_t>> arr; // raw log
    size_t num_updates;
    char*** chunks_compression; // zlib chunk
    size_t num_chunks; // number of chunks in each row
    size_t** compressed_sizes; // sizes of compressed chunks

    uint64_t hash(uint64_t key, uint32_t seed) const {
        uint64_t out[2];
        MurmurHash3_x64_128(&key, sizeof(uint64_t), seed, out);
        return out[1];
    }
};

CountMin::CountMin(double eps, double delta, uint64_t seed, storage_type type = Uncompressed, uint64_t number_updates = 100, const char* filename = NULL, const char* name = "noname") {
    this->type = type;
//    this->name = name;
    w = ceil(M_E / eps);
    d = ceil(log(1 / delta));
    num_contents = 0;
    mt19937_64 mt(seed);
    uniform_int_distribution<uint32_t> dist(0, (uint32_t)-1);
    //printf("%s-> w: %" PRIu64 ", d: %" PRIu64 "\n", name, w, d);
    if (type != BufferedHash and type != BufferedTree and type != BufferedRaw) {
        this->hash_seed.resize(d);
        for (size_t i = 0; i < d; i++) {
            this->hash_seed[i] = dist(mt);        // generate random seeds for hash function
        }
    }
    if (type == Uncompressed) {
        flatcounts.reset(filename, d*w);
    }
    else if (type == HashTable) {
        assert(filename == NULL);
        hash_table_counts.resize(d, Hashtable(2*number_updates, 1337));
    }
    else if (type == Tree) {
        assert(filename == NULL);
        tree_counts.assign(d, map<uint64_t,int>());
    }
    else if (type == ChunksZlib || type == LZ4) {
        assert(filename == NULL);
        num_chunks = (w / CHUNKSIZE);
        if (w % CHUNKSIZE != 0) {
            num_chunks++;
        }
        // printf("num chunks:%d\n", num_chunks);
        chunks_compression = new char**[d]();
        compressed_sizes = new size_t*[d]();
        for (size_t i = 0; i < d; i++) {
            chunks_compression[i] = new char*[num_chunks]();
            compressed_sizes[i] = new size_t[num_chunks]();
        }
    }
    // else if (type == LZ4) {
    //     assert(filename == NULL);
    //     char* tmp = new char[d*w*sizeof(int)]();
    //     lz4_size = lz4_compress(tmp, lz4_flatcounts, d*w*sizeof(int));
    //     delete[] tmp;
    // }
    else if (type == BufferedHash || type == BufferedTree || type == BufferedRaw) {
        assert(filename != NULL);
        if (w*d < BLOCKSIZE) {
            fprintf(stderr, "Filter is too small for buffered version\n");
            exit(0);
        }

        bchunk_hash_seed = dist(mt);
        num_chunks = ((w*d+BLOCKSIZE-1) / BLOCKSIZE);
        char tmp_filename[strlen(filename)+20];
        bchunks[0].reserve(num_chunks);
        for (size_t i=0; i<num_chunks; ++i) {
            sprintf(tmp_filename, "%s%lu", filename, i);
            uint32_t tmp_seed = dist(mt);
            // TODO how should eps increase?
            bchunks[0].emplace_back(d*M_E/BLOCKSIZE, delta, tmp_seed, Uncompressed,
                    number_updates, (const char*)tmp_filename);
            if (type == BufferedHash) {
                bchunks[1].emplace_back(d*M_E/BLOCKSIZE, delta, tmp_seed, HashTable, BUFFERSIZE);
            }
            else if (type == BufferedTree) {
                bchunks[1].emplace_back(d*M_E/BLOCKSIZE, delta, tmp_seed, Tree, BUFFERSIZE);
            }
            else {
                bchunks[1].emplace_back(d*M_E/BLOCKSIZE, delta, tmp_seed, RawLog, BUFFERSIZE);
            }
        }
    }
    else if (type == RawLog) {
        assert(filename == NULL);
        arr.resize(number_updates);
        num_updates = 0;
    }
    else {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
}

CountMin::~CountMin() {
}

void CountMin::update(uint64_t i, int c) {
    ++num_contents;
    if (type == Uncompressed) {
        for (size_t j = 0; j < d; j++) {
            flatcounts[j*w + hash(i, hash_seed[j]) % w] += c;
        }
    }
    else if (type == HashTable) {
        for (size_t j = 0; j < d; ++j) {
            hash_table_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    }
    else if (type == Tree) {
        for (size_t j = 0; j < d; ++j) {
            tree_counts[j][hash(i, hash_seed[j]) % w] += c;
        }
    }
    else if (type == ChunksZlib || type == LZ4) {
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
            if (chunks_compression[j][chunk_block]) {
                inflated = decompression(chunks_compression[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int), type);
                delete[] chunks_compression[j][chunk_block];
            }
            else {
                inflated = new int[CHUNKSIZE]();
            }
            inflated[chunk_index] += c;
            compressed_sizes[j][chunk_block] = compression(inflated, chunk_size * sizeof(int), &chunks_compression[j][chunk_block], type);
            delete[] inflated;
        }
    }
    else if (type == BufferedHash || type == BufferedTree || type == BufferedRaw) {
        size_t chunk = hash(i, bchunk_hash_seed)%num_chunks;
        if(bchunks[1][chunk].getNumContents() == BUFFERSIZE) {
            bchunks[0][chunk].mergeCMs(bchunks[1][chunk]);
            bchunks[1][chunk].clear();
        }
        bchunks[1][chunk].update(i, c);
    }
    else if (type == RawLog) {
        arr[num_updates] = {i, c};
        num_updates++;
    }
    else {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
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
    if (type == ChunksZlib || type == LZ4) {
        for (size_t j = 0; j < d; j++) {
            size_t index =  hash(i, hash_seed[j]) % w; // find out where the hash lives
            size_t chunk_block = (index / CHUNKSIZE); // find out which chunk block it is in
            size_t chunk_index = index % CHUNKSIZE; // find out which chunk index it is in the block
            size_t chunk_size = CHUNKSIZE; // size of uncompressed block
            if (!chunks_compression[j][chunk_block]) return 0;
            // if it's the last block it may be smaller
            if (chunk_block == num_chunks - 1) {
                chunk_size = w % CHUNKSIZE;
            }
            int* inflated;
            inflated = decompression(chunks_compression[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int), type);
            if (j == 0) res = inflated[chunk_index];
            res = min(res, inflated[chunk_index]);
            delete[] inflated;
        }
        return res;
    }
    // if (type == LZ4) {
    //     char* tmp = new char[d*w*sizeof(int)]();
    //     lz4_decompress(lz4_flatcounts, tmp, lz4_size, d*w*sizeof(int));
    //     res = ((int*)tmp)[hash(i, hash_seed[0]) % w];
    //     for (size_t j = 1; j < d; j++) {
    //         res = min(res, ((int*)tmp)[j*w + hash(i, hash_seed[j]) % w]);
    //     }
    //     delete[] tmp;
    //     return res;
    // }
    if (type == BufferedHash || type == BufferedTree || type == BufferedRaw) {
        size_t chunk = hash(i, bchunk_hash_seed)%num_chunks;
        return bchunks[0][chunk].pointQuery(i) + bchunks[1][chunk].pointQuery(i);
    }
    if (type == RawLog) {
        for(size_t j=0; j<num_updates; ++j) {
            if(i==arr[j].first) {
                res += arr[j].second;
            }
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
            const Array &otherCounts = other.getCounts();
            for (size_t j = 0; j < d; j++) {
                for (size_t i = 0; i < w; i++) {
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
        if (other.type == ChunksZlib || other.type == LZ4) {
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
                        inflated = decompression(chunks[i][j], sizes[i][j], chunk_size * sizeof(int), other.type);
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
    else if ((type == ChunksZlib && other.type == ChunksZlib)|| (type == LZ4 && other.type == LZ4)) {
        char*** chunks = other.getChunkZlibCompressions();
        size_t** sizes = other.getCompressedSizes();
        size_t chunk_size = CHUNKSIZE;
        for (size_t j = 0; j < num_chunks; j++) {
            for (size_t i = 0; i < d; i++) {
                chunk_size = CHUNKSIZE;
                if (chunks_compression[i][j] && chunks[i][j]) {
                    if (j == num_chunks - 1) {
                        chunk_size = w % CHUNKSIZE;
                    }
                    int* inflated1;
                    int* inflated2;
                    inflated1 = decompression(chunks[i][j], sizes[i][j], chunk_size * sizeof(int), type);
                    inflated2 = decompression(chunks_compression[i][j], compressed_sizes[i][j], chunk_size * sizeof(int), type);
                    for (size_t k = 0; k < chunk_size; k++) {
                        totals[j] += inflated1[k] * inflated2[k];
                    }
                    delete[] inflated1;
                    delete[] inflated2;
                }
            }
        }
    }
    else if (other.type == RawLog) {
        const auto &rawlog = other.getRawLog();
        size_t updates = other.getNumUpdates();
        if (type == Uncompressed) {
            for(size_t i=0; i<updates; ++i) {
                for (size_t j = 0; j < d; j++) {
                    size_t index = hash(rawlog[i].first, hash_seed[j]) % w;
                    totals[j] += flatcounts[j * w + index] * rawlog[i].second;
                }
            }
        }
        else if (type == HashTable) {
            for(size_t i=0; i<updates; ++i) {
                for (size_t j = 0; j < d; j++) {
                    size_t index = hash(rawlog[i].first, hash_seed[j]) % w;
                    auto it = hash_table_counts[j].find(index);
                    totals[j] += max(0, it) * rawlog[i].second;
                }
            }
        }
        else if (type == Tree) {
            for(size_t i=0; i<updates; ++i) {
                for (size_t j = 0; j < d; j++) {
                    size_t index = hash(rawlog[i].first, hash_seed[j]) % w;
                    auto it = tree_counts[j].find(index);
                    if (it != tree_counts[j].end()) {
                        totals[j] += it->second * rawlog[i].second;
                    }
                }
            }
        }
        else if (type == ChunksZlib || type == LZ4) {
            for(size_t i=0; i<updates; ++i) {
                for (size_t j = 0; j < d; j++) {
                    size_t index = hash(rawlog[i].first, hash_seed[j]) % w;
                    size_t chunk_block = (index / CHUNKSIZE); // find out which chunk block it is in
                    size_t chunk_index = index % CHUNKSIZE; // find out which chunk index it is in the block
                    size_t chunk_size = CHUNKSIZE; // size of uncompressed block
                    if (chunks_compression[j][chunk_block]) {
                        if (chunk_block == num_chunks - 1) {
                            chunk_size = w % CHUNKSIZE;
                        }
                        int* inflated;
                        inflated = decompression(chunks_compression[j][chunk_block], compressed_sizes[j][chunk_block], chunk_size * sizeof(int), type);
                        totals += inflated[chunk_index] * rawlog[i].second;
                        delete[] inflated;
                    }
                }
            }
        }
        else if ( (type == BufferedHash || type == BufferedTree || type == BufferedRaw)
                && (other.type == BufferedHash || other.type == BufferedTree || other.type == BufferedRaw) ) {
            assert(num_chunks == other.getNumChunks());
            const vector<CountMin>* otherBchunks = other.getBchunks();
            int res = 0;
            for(size_t k=0; k<num_chunks; ++k) {
                res += bchunks[0][k].innerProductQuery(otherBchunks[0][k]);
                res += bchunks[0][k].innerProductQuery(otherBchunks[1][k]);
                res += bchunks[1][k].innerProductQuery(otherBchunks[0][k]);
                res += bchunks[1][k].innerProductQuery(otherBchunks[1][k]);
            }
            delete[] totals;
            return res;
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

void CountMin::clear() {
    if (type == Uncompressed) {
        flatcounts.clear();
        num_contents = 0;
    }
    else if (type == HashTable) {
        for (auto &i: hash_table_counts) {
            i.clear();
        }
        num_contents = 0;
    }
    else if (type == Tree) {
        for (auto &i: tree_counts) {
            i.clear();
        }
        num_contents = 0;
    }
    else {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
}

void CountMin::mergeCMs(const CountMin& other) {
    if (type != Uncompressed) {
        fprintf(stderr, "NOT IMPLEMENTED\n");
        exit(0);
    }
    if (other.type == Uncompressed) {
        const Array &otherCounts = other.getCounts();
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
    if (other.type == ChunksZlib || other.type == LZ4) {
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
                    inflated = decompression(chunks[i][j], sizes[i][j], chunk_size * sizeof(int), other.type);
                    for (size_t k = 0; k < chunk_size; k++) {
                        flatcounts[i * w + CHUNKSIZE * j + k] += inflated[k];
                    }
                    delete[] inflated;
                }
            }
        }
        return;
    }
    if (other.type == RawLog) {
        const auto &rawlog = other.getRawLog();
        for (size_t j = 0; j < other.getNumUpdates(); j++) {
            update(rawlog[j].first, rawlog[j].second);
        }
        return;
    }

    fprintf(stderr, "NOT IMPLEMENTED\n");
    exit(0);
}

#endif
