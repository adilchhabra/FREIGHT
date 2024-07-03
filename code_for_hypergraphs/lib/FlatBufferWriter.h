/******************************************************************************
 * FlatBufferWriter.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

#ifndef KAHIP_FLATBUFFERWRITER_H
#define KAHIP_FLATBUFFERWRITER_H

#include <fstream>
#include <iostream>
#include <vector>
#include "Freight_Info_generated.h"
#include "partition/partition_config.h"

class FlatBufferWriter {
private:
    double buffer_io_time_;
    double global_mapping_time_;
    double total_time_;
    long maxRSS_;

    double cutNet_;
    double connectivity_;
    double balance_;

public:
    FlatBufferWriter()
            : buffer_io_time_(0.0), global_mapping_time_(0.0),
              total_time_(0.0), maxRSS_(0), cutNet_(0.0),
              connectivity_(0.0), balance_(0.0){}

    void updateResourceConsumption(double &buffer_io_time,
                                   double &global_mapping_time, double &total_time, long &maxRSS) {
        buffer_io_time_ = buffer_io_time;
        global_mapping_time_ = global_mapping_time;
        total_time_ = total_time;
        maxRSS_ = maxRSS;
    }

    void updateHypergraphPartitionMetrics(double &cutNet, double &connectivity,
                                          double &balance) {
        cutNet_ = cutNet;
        connectivity_ = connectivity;
        balance_ = balance;
    }

    // Function to extract the base filename without path and extension
    static std::string extractBaseFilename(const std::string& fullPath) {
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

    void write(const std::string& baseFilename, const PartitionConfig &config ) const {
        // output some information about the partition that we have computed
        flatbuffers::FlatBufferBuilder builder(1024);
        FreightInfo::GraphMetadataBuilder metadata_builder(builder);
        auto filenameOffset = builder.CreateString(baseFilename);
        metadata_builder.add_filename(filenameOffset);
        std::cout << "(Hyper)Graph: " << baseFilename << std::endl;
        if (config.edge_partition) {
            metadata_builder.add_num_nodes(config.total_edges);
            metadata_builder.add_num_edges(config.total_nodes);
            std::cout << "Nodes: " << config.total_edges << std::endl;
            std::cout << "Edges: " << config.total_nodes << std::endl;
        } else {
            metadata_builder.add_num_nodes(config.total_nodes);
            metadata_builder.add_num_edges(config.total_edges);
            std::cout << "Nodes: " << config.total_edges << std::endl;
            std::cout << "Edges: " << config.total_nodes << std::endl;
        }
        auto metadata = metadata_builder.Finish();
        builder.Finish(metadata);

        FreightInfo::PartitionConfigurationBuilder config_builder(builder);
        config_builder.add_k(config.k);
        config_builder.add_seed(config.seed);
        config_builder.add_input_balance(config.imbalance);
        if(config.edge_partition) {
            config_builder.add_edge_partition(1);
            std::cout << "Partition: Edge" << std::endl;
        } else {
            config_builder.add_edge_partition(0);
            std::cout << "Partition: Node" << std::endl;
        }
        config_builder.add_rle_length(config.rle_length);
        config_builder.add_kappa(config.kappa);
        std::cout << "Blocks (k): " << config.k << std::endl;
        std::cout << "Seed: " << config.seed << std::endl;
        std::cout << "Imbalance: " << config.imbalance << std::endl;
        if (config.rle_length > -1) {
            std::cout << "RLE Length: " << config.rle_length << std::endl;
        } else if (config.rle_length == -2) {
            std::cout << "External Memory Priority Queue" << std::endl;
        } else {
            std::cout << "RLE Length: None" << std::endl;
        }
        std::cout << "Kappa: " << config.kappa << std::endl;
        auto configdata = config_builder.Finish();
        builder.Finish(configdata);

        FreightInfo::RunTimeBuilder runtime_builder(builder);
        runtime_builder.add_io_time(buffer_io_time_);
        runtime_builder.add_mapping_time(global_mapping_time_);
        runtime_builder.add_total_time(total_time_);
        auto runtimedata = runtime_builder.Finish();
        builder.Finish(runtimedata);
        std::cout << "IO Time: " << buffer_io_time_ << std::endl;
        std::cout << "Mapping Time: " << global_mapping_time_ << std::endl;
        std::cout << "Total Time: " << total_time_ << std::endl;

        auto partition_metrics = FreightInfo::CreatePartitionMetrics(builder,
                                                                     cutNet_,
                                                                     connectivity_,
                                                                     0,
                                                                     -1,
                                                                     -1,
                                                                     0,
                                                                     balance_);
        if(!config.edge_partition) {
            std::cout << "Cut: " << cutNet_ << std::endl;
            std::cout << "Connectivity: " << connectivity_ << std::endl;
            std::cout << "Balance: " << balance_ << std::endl;
        }

        // Create MemoryConsumption
        if (maxRSS_ != -1) {
            std::cout << "Maximum Resident Set Size (KB): " << maxRSS_ << std::endl;
        }
        auto memory_consumption = FreightInfo::CreateMemoryConsumption(builder, maxRSS_);

        // Create Output
        FreightInfo::PartitionLogBuilder partition_builder(builder);
        partition_builder.add_graph_metadata(metadata);
        partition_builder.add_partition_configuration(configdata);
        partition_builder.add_runtime(runtimedata);
        partition_builder.add_memory_consumption(memory_consumption);
        partition_builder.add_metrics(partition_metrics);
        auto partition = partition_builder.Finish();
        builder.Finish(partition);

        //Step 4: Write to File
        const uint8_t *bufferPointer = builder.GetBufferPointer();
        int bufferSize = builder.GetSize();

        std::string outputFileNameStream;
        outputFileNameStream = config.output_path + baseFilename + "_" + std::to_string(config.k) + ".bin";
        const char *outputFileName = outputFileNameStream.c_str();
        if(config.write_results) {
            FILE *file = fopen(outputFileName, "wb");
            fwrite(bufferPointer, 1, bufferSize, file);
            fclose(file);
        }
    }
};

#endif //KAHIP_FLATBUFFERWRITER_H
