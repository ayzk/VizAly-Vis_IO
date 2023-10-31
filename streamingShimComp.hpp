#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <fstream>
#include <cmath>


#include "blosc.h"
#include "sz.h"
#include "sz_stats.h"
#include "SZ3/api/sz.hpp"

//
// Generic memory deallocator
inline int deAllocMem(std::string dataType, void *&data) {
    if (data == NULL) // already deallocated!
        return 1;

    if (dataType == "float")
        delete[](float *) data;
    else if (dataType == "double")
        delete[](double *) data;
    else if (dataType == "int")
        delete[](int *) data;
    else if (dataType == "int8_t")
        delete[](int8_t *) data;
    else if (dataType == "int16_t")
        delete[](int16_t *) data;
    else if (dataType == "int32_t")
        delete[](int32_t *) data;
    else if (dataType == "int64_t")
        delete[](int64_t *) data;
    else if (dataType == "uint8_t")
        delete[](uint8_t *) data;
    else if (dataType == "uint16_t")
        delete[](uint16_t *) data;
    else if (dataType == "uint32_t")
        delete[](uint32_t *) data;
    else if (dataType == "uint64_t")
        delete[](uint64_t *) data;
    else
        return 0;

    data = NULL;

    return 1;
}


class StreamingShimComp {
    u_int32_t headerSize;
    u_int32_t numWrites;
    size_t blockIndex;

    std::stringstream log;
    size_t huffmanTreeSize;


    std::vector<u_int32_t> blkOffset;    // in bytes
    std::vector<u_int32_t> blkCmpSize;   // in bytes
    std::vector<u_int32_t> blkNumRows;

    std::map<std::string, std::string> params; // Key: <value, data type>

public:
    StreamingShimComp() {
        blockIndex = 0;
        huffmanTreeSize = 0;
    };

    ~StreamingShimComp();

    void setParam(std::string key, std::string val);

    size_t stmCompress(float *input, std::string dataType, size_t dataTypeSize, size_t n, int blockSize);

    size_t compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t numToCompress);

    int stmDecompressHeader();

    int stmDecompressNextBlock(void *&data);

    int stmDecompressData();

    size_t getHuffSize() { return huffmanTreeSize; }

private:
    size_t stmCompress_SZ3(float *input, std::string dataType, size_t dataTypeSize, size_t n, int blockSize);
};


inline ::StreamingShimComp::~StreamingShimComp() {
    std::ofstream logFile;
    logFile.open("output.log");
    logFile << log.str();
    logFile.close();

    std::cout << "huffman tree size: " << huffmanTreeSize << std::endl;
};


inline void StreamingShimComp::setParam(std::string key, std::string val) {
    params[key] = val;
}


inline size_t StreamingShimComp::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t numToCompress) {
    if (params["compressor"] == "SZ") {
        int dataType = SZ_FLOAT;
        int mode = ABS;
        double pwRel, abs, rel;
        pwRel = abs = rel = 0.0;

        if (params.count("abs") != 0) {
            mode = ABS;
            abs = std::stof(params["abs"]);
        } else if (params.count("rel") != 0) {
            mode = REL;
            rel = std::stof(params["rel"]);
        } else if (params.count("pwRel") != 0) {
            mode = PW_REL;
            pwRel = std::stof(params["pwRel"]);
        }

        if (params.count("type") != 0)
            if (params["type"] == "double")
                dataType = SZ_DOUBLE;

        std::size_t csize;
        size_t n[5] = {0, 0, 0, 0, 0};

        std::uint8_t *cdata = SZ_compress_args(dataType, (float *) input, &csize, mode, abs, rel, pwRel, n[4], n[3], n[2], n[1], numToCompress);
        output = cdata;
        return csize;
    } else if (params["compressor"] == "BLOSC") {
        size_t isize = dataTypeSize * numToCompress;
        size_t osize = isize + BLOSC_MAX_OVERHEAD;

        output = std::malloc(isize); //byte array;
        osize = blosc_compress(9, 1, dataTypeSize, isize, input, output, osize);

        if (osize < 0)
            throw std::runtime_error("Compression error. Error code: " + std::to_string(osize));

        if (osize > 0)
            output = std::realloc(output, osize);

        return osize;
    }
    return -1;
}


// Assume all the data is float
inline size_t StreamingShimComp::stmCompress(float *input, std::string dataType, size_t dataTypeSize, size_t numElements, int blockSize) {
    if (params["compressor"] == "SZ3-small-tree" || params["compressor"] == "SZ3-shared-tree" || params["compressor"] == "SZ3-first-tree") {
        return stmCompress_SZ3(input, dataType, dataTypeSize, numElements, blockSize);
    }

    u_int32_t numWrites = std::ceil((float) numElements / blockSize);
    std::vector<u_int32_t> cmpBlkOffset(numWrites);      // in bytes
    std::vector<u_int32_t> cmpBlkCmpSize(numWrites);     // in bytes
    std::vector<u_int32_t> cmpBlkNumRows(numWrites);
    u_int32_t cmpHeaderSize = sizeof(u_int32_t) * 2 + (sizeof(u_int32_t) * 3 * numWrites); // headerSize, numRows, header (offset, size, #rows)


    // Open file for output
    std::ofstream outputFile("output.gio", std::ios::binary);

    // Header
    outputFile.write(reinterpret_cast<const char *>(&cmpHeaderSize), sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(&numWrites), sizeof(u_int32_t));

    outputFile.write(reinterpret_cast<const char *>(cmpBlkOffset.data()), cmpBlkOffset.size() * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));



    // Compress Blocks
    std::size_t csize = 0;
    std::size_t pos = 0;
    int iteration = 0;
    std::size_t totalCompSize = 0;

    while (pos < numElements) {
        // Find how many particles this block has
        size_t numToCompress = blockSize;
        if (numToCompress + pos > numElements)
            numToCompress = numElements - pos;

        // Compress
        void *output;
        std::size_t csize = compress(&input[pos], output, dataType, dataTypeSize, numToCompress);

        //log << "huffman tree size = " << sz_stat.huffmanTreeSize << " bytes " << std::endl;
//        huffmanTreeSize += sz_stat.huffmanTreeSize;
//        log << sz_stat.huffmanTreeSize << std::endl;

        // Write
        outputFile.write(reinterpret_cast<const char *>(output), csize);
        deAllocMem(dataType, output);

        // Store info
        if (iteration == 0)
            cmpBlkOffset[iteration] = 0;
        else
            cmpBlkOffset[iteration] = cmpBlkOffset[iteration - 1] + cmpBlkCmpSize[iteration - 1];
        cmpBlkCmpSize[iteration] = csize;
        cmpBlkNumRows[iteration] = numToCompress;


        pos += numToCompress;
        totalCompSize += csize;
        iteration++;
    }
    std::cout << ", # blocks: " << numWrites;
    outputFile.close();


    //
    // Reopen file to modify offset and size
    std::fstream of2("output.gio", std::fstream::in | std::fstream::out | std::fstream::binary);
    of2.seekp(2 * sizeof(u_int32_t), std::ios::beg);
    of2.write(reinterpret_cast<const char *>(cmpBlkOffset.data()), cmpBlkOffset.size() * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char *>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char *>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));
    of2.close();


    return totalCompSize;
}

inline size_t StreamingShimComp::stmCompress_SZ3(float *input, std::string dataType, size_t dataTypeSize, size_t numElements, int blockSize) {
    if (params.count("type") != 0 && params["type"] == "double") {
        throw "SZ3 double not added yet";
    }
    u_int32_t numWrites = std::ceil((float) numElements / blockSize);
    std::vector<u_int32_t> cmpBlkOffset(numWrites);      // in bytes
    std::vector<u_int32_t> cmpBlkCmpSize(numWrites);     // in bytes
    std::vector<u_int32_t> cmpBlkNumRows(numWrites);
    u_int32_t cmpHeaderSize = sizeof(u_int32_t) * 2 + (sizeof(u_int32_t) * 3 * numWrites); // headerSize, numRows, header (offset, size, #rows)


    // Open file for output
    std::ofstream outputFile("output.gio", std::ios::binary);

    // Header
    outputFile.write(reinterpret_cast<const char *>(&cmpHeaderSize), sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(&numWrites), sizeof(u_int32_t));

    outputFile.write(reinterpret_cast<const char *>(cmpBlkOffset.data()), cmpBlkOffset.size() * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char *>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));

    SZ::HuffmanEncoder<int> encoder;
    bool isFirstComprIter = true;

    size_t bufferSize = 2 * sizeof(double) * numElements;
    auto *buffer = new SZ::uchar[bufferSize];
    SZ::uchar *buffer_pos = buffer;
    std::vector<int> quant_inds(numElements);

    // Compress Blocks
    std::size_t pos = 0;
    std::size_t csize = 0;
    int iteration = 0;
    std::size_t totalCompSize = 0;

    while (pos < numElements) {
        // Find how many particles this block has
        size_t numToCompress = blockSize;
        if (numToCompress + pos > numElements)
            numToCompress = numElements - pos;

        SZ::Config conf(numToCompress);
        conf.cmprAlgo = SZ::ALGO_LORENZO_REG;
        conf.lorenzo = true;
        conf.regression = false;
        if (params.count("abs") != 0) {
            conf.errorBoundMode = SZ::EB_ABS;
            conf.absErrorBound = std::stof(params["abs"]);
        } else if (params.count("rel") != 0) {
            conf.errorBoundMode = SZ::EB_ABS;
            conf.absErrorBound = std::stof(params["abs"]);
        } else {
            throw "error bound not supported";
        }
        void *output;
        std::vector<float> data_cpy(&input[pos], &input[pos + numToCompress]);
        SZ::calAbsErrorBound(conf, data_cpy.data());
        if (params["compressor"] == "SZ3-small-tree") {
            output = SZ_compress<float>(conf, &input[pos], csize);
        } else {
            std::vector<int> quant_inds_block;
            auto predictor = SZ::LorenzoPredictor<float, 1, 1>(conf.absErrorBound);
            auto quantizer = SZ::LinearQuantizer<float>(conf.absErrorBound, conf.quantbinCnt / 2);
            auto frontend = SZ::make_sz_general_frontend<float, 1>(conf, predictor, quantizer);
            quant_inds_block = frontend.compress(data_cpy.data());

            if (params["compressor"] == "SZ3-first-tree") {
                buffer_pos = buffer;
                if (isFirstComprIter) {
                    std::vector<bool> quant_ind_map(conf.quantbinCnt, false);
                    for (auto &quant: quant_inds_block) {
                        quant_ind_map[quant] = true;
                    }
                    std::vector<int> quant_inds_extended(quant_inds_block);
                    quant_inds_extended.reserve(conf.quantbinCnt);
                    for (int i = 0; i < conf.quantbinCnt; i++) {
                        if (!quant_ind_map[i]) {
                            quant_inds_extended.push_back(i);
                        }
                    }
                    encoder.preprocess_encode(quant_inds_extended, 0);

                    encoder.save(buffer_pos);
                }

                encoder.encode(quant_inds_block, buffer_pos);

                SZ::Lossless_zstd zstd;
                output = zstd.compress(buffer, buffer_pos - buffer, csize);
                isFirstComprIter = false;

            } else if (params["compressor"] == "SZ3-shared-tree") {
                memcpy(&quant_inds[pos], quant_inds_block.data(), sizeof(int) * numToCompress);
            }
        }
        if (params["compressor"] == "SZ3-small-tree" || params["compressor"] == "SZ3-first-tree") {
            // Write
            outputFile.write(reinterpret_cast<const char *>(output), csize);
            deAllocMem(dataType, output);

            // Store info
            if (iteration == 0)
                cmpBlkOffset[iteration] = 0;
            else
                cmpBlkOffset[iteration] = cmpBlkOffset[iteration - 1] + cmpBlkCmpSize[iteration - 1];
            cmpBlkCmpSize[iteration] = csize;
            cmpBlkNumRows[iteration] = numToCompress;
        }
        totalCompSize += csize;
        iteration++;
        pos += numToCompress;
    }


    if (params["compressor"] == "SZ3-shared-tree") {
        // Compress Blocks
        pos = 0;
        iteration = 0;
        buffer_pos = buffer;

        //build huffman tree
        SZ::HuffmanEncoder<int> huffman;
        huffman.preprocess_encode(quant_inds, 0);

        //write huffman tree to file
        huffman.save(buffer_pos);
        outputFile.write(reinterpret_cast<const char *>(buffer), buffer_pos - buffer);

        while (pos < numElements) {
            // Find how many particles this block has
            size_t numToCompress = blockSize;
            if (numToCompress + pos > numElements)
                numToCompress = numElements - pos;


            //encode each block and zstd it
            SZ::Lossless_zstd zstd;
            buffer_pos = buffer;
            huffman.encode(&quant_inds[pos], numToCompress, buffer_pos);
            void *output = zstd.compress(buffer, buffer_pos - buffer, csize);

            // Write
            outputFile.write(reinterpret_cast<const char *>(output), csize);
            deAllocMem(dataType, output);

            // Store info
            if (iteration == 0)
                cmpBlkOffset[iteration] = 0;
            else
                cmpBlkOffset[iteration] = cmpBlkOffset[iteration - 1] + cmpBlkCmpSize[iteration - 1];
            cmpBlkCmpSize[iteration] = csize;
            cmpBlkNumRows[iteration] = numToCompress;


            pos += numToCompress;
            totalCompSize += csize;
            iteration++;
        }
    }


    delete[] buffer;
    std::cout << ", # blocks: " << numWrites;
    outputFile.close();
    //
    // Reopen file to modify offset and size
    std::fstream of2("output.gio", std::fstream::in | std::fstream::out | std::fstream::binary);
    of2.seekp(2 * sizeof(u_int32_t), std::ios::beg);
    of2.write(reinterpret_cast<const char *>(cmpBlkOffset.data()), cmpBlkOffset.size() * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char *>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char *>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));
    of2.close();


    return totalCompSize;
}


/*
inline size_t StreamingShimComp::stmCompress(float *input, std::string dataType, size_t dataTypeSize, size_t * n, int blockSize)
{
    size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

    u_int32_t numWrites = std::ceil((float)numel/blockSize);
    std::vector<u_int32_t> cmpBlkOffset(numWrites);      // in bytes
    std::vector<u_int32_t> cmpBlkCmpSize(numWrites);     // in bytes
    std::vector<u_int32_t> cmpBlkNumRows(numWrites);     
    u_int32_t cmpHeaderSize = sizeof(u_int32_t)*2 + (sizeof(u_int32_t)*3 * numWrites); // headerSize, numRows, header (offset, size, #rows)


    // Open file for output
    std::ofstream outputFile("output.gio", std::ios::binary);

    // Header
    outputFile.write(reinterpret_cast<const char*>(&cmpHeaderSize), sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char*>(&numWrites),  sizeof(u_int32_t));

    outputFile.write(reinterpret_cast<const char*>(cmpBlkOffset.data()),  cmpBlkOffset.size()  * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char*>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    outputFile.write(reinterpret_cast<const char*>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));



    // Compress Blocks
    int mode = ABS;
    double pwRel = 0.0;
    double abs = 0.003;
    double rel = 0.0;

    std::size_t csize = 0;
    std::size_t pos = 0;
    int iteration = 0;
    std::size_t totalCompSize = 0;

    while (pos < numel)
    {   
        // Find how many particles this block has
        size_t numToCompress = blockSize;
        if (numToCompress + pos > numel)
            numToCompress = numel - pos;

        // Compress
	    std::uint8_t *cdata = SZ_compress_args(SZ_FLOAT, &input[pos], &csize, mode, abs, rel, pwRel, n[4], n[3], n[2], n[1], numToCompress);

        // Write
        outputFile.write(reinterpret_cast<const char*>(cdata), csize);
        delete []cdata;

        // Store info
        if (iteration == 0)
            cmpBlkOffset[iteration] = 0;
        else
            cmpBlkOffset[iteration] = cmpBlkOffset[iteration-1] + cmpBlkCmpSize[iteration-1];
        cmpBlkCmpSize[iteration] = csize;
        cmpBlkNumRows[iteration] = numToCompress;

 
        pos += numToCompress;
        totalCompSize += csize;
        iteration++;        
    }
    std::cout << ", # blocks: " << numWrites ;
    outputFile.close();


    //
    // Reopen file to modify offset and size
    std::fstream of2("output.gio", std::fstream::in | std::fstream::out | std::fstream::binary);
    of2.seekp(2*sizeof(u_int32_t), std::ios::beg);
    of2.write(reinterpret_cast<const char*>(cmpBlkOffset.data()),  cmpBlkOffset.size()  * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char*>(cmpBlkCmpSize.data()), cmpBlkCmpSize.size() * sizeof(u_int32_t));
    of2.write(reinterpret_cast<const char*>(cmpBlkNumRows.data()), cmpBlkNumRows.size() * sizeof(u_int32_t));
    of2.close();
    
    
    return totalCompSize;
}
*/


inline int StreamingShimComp::stmDecompressHeader() {
    // Open file for reading
    std::ifstream inputFile("output.gio", std::ios::binary);

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        return -1; // or handle the error in an appropriate way
    }


    //
    // Get header information
    inputFile.read(reinterpret_cast<char *>(&headerSize), sizeof(u_int32_t));
    inputFile.read(reinterpret_cast<char *>(&numWrites), sizeof(u_int32_t));

    blkOffset.resize(numWrites);
    blkCmpSize.resize(numWrites);
    blkNumRows.resize(numWrites);

    inputFile.read(reinterpret_cast<char *>(blkOffset.data()), blkOffset.size() * sizeof(u_int32_t));
    inputFile.read(reinterpret_cast<char *>(blkCmpSize.data()), blkCmpSize.size() * sizeof(u_int32_t));
    inputFile.read(reinterpret_cast<char *>(blkNumRows.data()), blkNumRows.size() * sizeof(u_int32_t));

    inputFile.close();


    return 1;
}


inline int StreamingShimComp::stmDecompressData() {
    // Open file for reading
    std::ifstream inputFile("output.gio", std::ios::binary);

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        return -1; // or handle the error in an appropriate way
    }

    inputFile.seekg(headerSize, std::ios::beg);

    //
    // Decompress each block
    size_t n[5] = {0};
    void *output;
    for (int i = 0; i < numWrites; i++) {
        // Read from file
        std::size_t csize = blkCmpSize[i];
        uint8_t *buffer = new uint8_t[csize];
        inputFile.read(reinterpret_cast<char *>(buffer), csize);


        // Decompress
        n[0] = blkNumRows[i];
        output = SZ_decompress(SZ_FLOAT, buffer, csize, n[4], n[3], n[2], n[1], n[0]);
        std::cout << i << ", " << ((float *) output)[0] << std::endl;

        delete[]buffer;
        delete[](float *) output;
    }

    inputFile.close();


    return 1;
}


inline int StreamingShimComp::stmDecompressNextBlock(void *&data) {
    // Open file for reading
    std::ifstream inputFile("output.gio", std::ios::binary);

    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        return -1; // or handle the error in an appropriate way
    }

    inputFile.seekg(headerSize + blkOffset[blockIndex], std::ios::beg);


    // Read from file
    std::size_t csize = blkCmpSize[blockIndex];
    uint8_t *buffer = new uint8_t[csize];
    inputFile.read(reinterpret_cast<char *>(buffer), csize);
    inputFile.close();


    // Decompress
    size_t n[5] = {0};
    n[0] = blkNumRows[blockIndex];
    data = SZ_decompress(SZ_FLOAT, buffer, csize, n[4], n[3], n[2], n[1], n[0]);

    delete[] buffer;
    blockIndex++;

    return n[0];
}

