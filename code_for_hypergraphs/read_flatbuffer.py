import sys
import flatbuffers

# Assuming flatbuffers generated code is in the FreightInfo directory
from FreightInfo.PartitionLog import PartitionLog
from FreightInfo.GraphMetadata import GraphMetadata
from FreightInfo.PartitionConfiguration import PartitionConfiguration
from FreightInfo.RunTime import RunTime
from FreightInfo.MemoryConsumption import MemoryConsumption
from FreightInfo.PartitionMetrics import PartitionMetrics

def read_partition_log(file_path):
    with open(file_path, 'rb') as f:
        data = f.read()

    # Create a FlatBuffers buffer
    buf = flatbuffers.Builder(len(data))
    buf.Bytes = bytearray(data)
    
    # Get the root object for PartitionLog
    partition_log = PartitionLog.GetRootAsPartitionLog(buf.Bytes, 0)
    
    # Read GraphMetadata
    graph_metadata = partition_log.GraphMetadata()
    filename = graph_metadata.Filename().decode('utf-8') if graph_metadata.Filename() else None
    num_nodes = graph_metadata.NumNodes()
    num_edges = graph_metadata.NumEdges()

    # Read PartitionConfiguration
    partition_config = partition_log.PartitionConfiguration()
    k = partition_config.K()
    seed = partition_config.Seed()
    input_balance = partition_config.InputBalance()
    edge_partition = partition_config.EdgePartition()
    rle_length = partition_config.RleLength()
    kappa = partition_config.Kappa()

    # Read RunTime
    runtime = partition_log.Runtime()
    io_time = runtime.IoTime()
    mapping_time = runtime.MappingTime()
    total_time = runtime.TotalTime()

    # Read MemoryConsumption
    memory_consumption = partition_log.MemoryConsumption()
    max_rss = memory_consumption.MaxRss()

    # Read PartitionMetrics
    metrics = partition_log.Metrics()
    cut_net = metrics.CutNet()
    connectivity = metrics.Connectivity()
    edge_cut = metrics.EdgeCut()
    vertex_cut = metrics.VertexCut()
    replicas = metrics.Replicas()
    replication_factor = metrics.ReplicationFactor()
    balance = metrics.Balance()

    # Print all the data
    print("Graph Metadata")
    print("==============")
    print(f"Filename: {filename}")
    print(f"Number of Nodes: {num_nodes}")
    print(f"Number of Edges: {num_edges}")
    print("\nPartition Configuration")
    print("========================")
    print(f"K: {k}")
    print(f"Seed: {seed}")
    print(f"Input Balance: {input_balance}")
    print(f"Edge Partition: {edge_partition}")
    print(f"RLE Length: {rle_length}")
    print(f"Kappa: {kappa}")
    print("\nRun Time")
    print("========")
    print(f"IO Time: {io_time}")
    print(f"Mapping Time: {mapping_time}")
    print(f"Total Time: {total_time}")
    print("\nMemory Consumption")
    print("===================")
    print(f"Max RSS: {max_rss}")
    print("\nPartition Metrics")
    print("=================")
    print(f"Cut Net: {cut_net}")
    print(f"Connectivity: {connectivity}")
    print(f"Edge Cut: {edge_cut}")
    print(f"Vertex Cut: {vertex_cut}")
    print(f"Replicas: {replicas}")
    print(f"Replication Factor: {replication_factor}")
    print(f"Balance: {balance}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python read_partition_log.py <file_path>")
    else:
        file_path = sys.argv[1]
        read_partition_log(file_path)
