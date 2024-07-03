/******************************************************************************
 * evaluator_ep.cpp
 * *
 * Adil Chhabra
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "data_structure/matrix/online_precalc_matrix.h"
#include "data_structure/matrix/online_binary_matrix.h"
#include "data_structure/matrix/full_matrix.h"
#include "graph_io_stream.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "tools/random_functions.h"
#include "timer.h"

#include "partition/onepass_partitioning/vertex_partitioning.h"
#include "partition/onepass_partitioning/fennel.h"
#include "partition/onepass_partitioning/fennel_approx_sqrt.h"
#include "partition/onepass_partitioning/ldg.h"

#include "Freight_Info_generated.h"

#define MIN(A,B) (((A)<(B))?(A):(B))
#define MAX(A,B) (((A)>(B))?(A):(B))

std::string extractBaseFilename(const std::string& fullPath);

int main(int argn, char **argv) {
    std::cout << R"(

    ╭━━━╮╱╱╱╱╱╱╱╱╭╮╱╭╮╭━━━╮
    ┃╭━━╯╱╱╱╱╱╱╱╱┃┃╭╯╰┫╭━━╯
    ┃╰━━┳━┳━━┳┳━━┫╰┻╮╭┫╰━━╮
    ┃╭━━┫╭┫┃━╋┫╭╮┃╭╮┃┃┃╭━━╯
    ┃┃╱╱┃┃┃┃━┫┃╰╯┃┃┃┃╰┫╰━━╮
    ╰╯╱╱╰╯╰━━┻┻━╮┣╯╰┻━┻━━━╯
    ╱╱╱╱╱╱╱╱╱╱╭━╯┃
    ╱╱╱╱╱╱╱╱╱╱╰━━╯

    )" << std::endl;
    PartitionConfig config;
    std::string graph_filename;
    quality_metrics qm;
    balance_configuration bc;

    bool is_graph_weighted = false;
    bool suppress_output   = false;
    bool recursive         = false;

    int ret_code = parse_parameters(argn, argv,
                                    config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if(ret_code) {
        return 0;
    }

    std::streambuf* backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if(suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }
    srand(config.seed);
    random_functions::setSeed(config.seed);

    config.LogDump(stdout);
    config.stream_input = true;

    bool already_fully_partitioned;

    std::cout << "Running evaluator..." << std::endl;

    // Open the binary FlatBuffer file for reading
    std::string baseFilename = extractBaseFilename(graph_filename);
    std::string outputFileNameStream;
    outputFileNameStream = config.output_path + baseFilename + "_" + std::to_string(config.k) + ".bin";

    std::fstream file(outputFileNameStream, std::ios::binary | std::ios::in | std::ios::out);

    if (!file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Read the contents of the file into a vector<char>
    file.seekg(0, std::ios::end);
    size_t length = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(length);
    file.read(buffer.data(), length);

    // Deserialize the FlatBuffer from the buffer
    flatbuffers::Verifier verifier(reinterpret_cast<const uint8_t*>(buffer.data()), length);
    if (!FreightInfo::VerifyEdgePartitionBuffer(verifier)) {
        std::cerr << "Error verifying FlatBuffer." << std::endl;
        return 1;
    }

    // Access the root table
    const FreightInfo::EdgePartition* edgePartition = FreightInfo::GetEdgePartition(buffer.data());

    // Access individual fields of the EdgePartition table
    const FreightInfo::GraphMetadata* metadata = edgePartition->graph_metadata();
    const std::string& filename = metadata->filename()->str();
    uint64_t numNodes = metadata->num_nodes();
    uint64_t numEdges = metadata->num_edges();
    std::cout << "(Hyper)Graph: " << filename << std::endl;
    std::cout << "Nodes: " << numNodes << std::endl;
    std::cout << "Edges: " << numEdges << std::endl;

    if (config.num_streams_passes >
        1 + config.restream_number) {
        config.stream_total_upperbound = ceil(
                ((100 + 1.5 * config.imbalance) / 100.) *
                (numEdges / (double)config.k));
    } else {
        config.stream_total_upperbound = ceil(
                ((100 + config.imbalance) / 100.) *
                (numEdges / (double)config.k));
    }

    // Access PartitionConfiguration
    const FreightInfo::PartitionConfiguration* config_data = edgePartition->partition_configuration();
    uint32_t k = config_data->k();
    int seed = config_data->seed();
    double imbalance = config_data->input_balance();
    std::cout << "Blocks (k): " << k << std::endl;
    std::cout << "Seed: " << seed << std::endl;
    std::cout << "Imbalance: " << imbalance << std::endl;

    // Access RunTime
    const FreightInfo::RunTime* runtime = edgePartition->runtime();
    double ioTime = runtime->io_time();
    double mappingTime = runtime->mapping_time();
    double totalTime = runtime->total_time();
    std::cout << "IO time: " << ioTime << std::endl;
    std::cout << "Mapping time: " << mappingTime << std::endl;
    std::cout << "Total time: " << totalTime << std::endl;

    // Access MemoryConsumption
    const FreightInfo::MemoryConsumption* memory_consumption = edgePartition->memory_consumption();
    long maxRSS = memory_consumption->max_rss();
    std::cout << "Maximum Resident Set Size (KB): " << maxRSS << std::endl;

    // Access PartitionMetrics
    const FreightInfo::PartitionMetrics* metrics = edgePartition->metrics();
    LongNodeID vertexCut = metrics->vertex_cut();
    LongNodeID replicas = metrics->replicas();
    double balance = metrics->balance();
    double replicationFactor = metrics->replication_factor();

    if(!config.filename_output.compare("")) {
        config.filename_output = "part_" + baseFilename + "_" + std::to_string(config.k) + ".txt";
    }

    if(replicationFactor == 0) {
        vertexCut = 0;
        replicas = 0;

        if(config.light_evaluator) {
            graph_io_stream::streamEvaluateEdgePartitionFast(
                    config, graph_filename, vertexCut, replicas, replicationFactor, balance);
        } else {
            graph_io_stream::streamEvaluateEdgePartition(
                    config, graph_filename, vertexCut, replicas, replicationFactor, balance);
        }

        // update metrics
        flatbuffers::FlatBufferBuilder builder(1024);
        auto updatedMetrics = FreightInfo::CreatePartitionMetrics(builder,
                                                                  vertexCut,
                                                                  replicas,
                                                                  replicationFactor,
                                                                  balance);
        // Create a new EdgePartition with the updated metrics
        auto filenameOffset = builder.CreateString(filename);
        auto updatedEdgePartition = FreightInfo::CreateEdgePartition(
                builder,
                FreightInfo::CreateGraphMetadata(builder, filenameOffset, edgePartition->graph_metadata()->num_nodes(),
                                                 edgePartition->graph_metadata()->num_edges()),
                FreightInfo::CreatePartitionConfiguration(builder, edgePartition->partition_configuration()->k(),
                                                          edgePartition->partition_configuration()->seed(),
                                                          edgePartition->partition_configuration()->input_balance()),
                FreightInfo::CreateRunTime(builder, edgePartition->runtime()->io_time(),
                                           edgePartition->runtime()->mapping_time(),
                                           edgePartition->runtime()->total_time()),
                FreightInfo::CreateMemoryConsumption(builder, edgePartition->memory_consumption()->max_rss()),
                updatedMetrics
        );

        // Finish updating the FlatBuffer
        builder.Finish(updatedEdgePartition);
        file.seekp(0);
        file.write(reinterpret_cast<const char *>(builder.GetBufferPointer()), builder.GetSize());
        file.close();
    }

    std::cout << "Vertex Cut: " << vertexCut << std::endl;
    std::cout << "Replicas: " << replicas << std::endl;
    std::cout << "Replication Factor: " << replicationFactor << std::endl;
    std::cout << "Balance: " << balance << std::endl;

    return 0;
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string& fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        // Found a slash, extract the substring after the last slash
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    } else {
        // No slash found, just extract the substring before the last dot
        return fullPath.substr(0, lastDot);
    }
}


