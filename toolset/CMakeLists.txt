ADD_EXECUTABLE(EdgeList2CSR.out edgelist_to_csr.cpp)
TARGET_LINK_LIBRARIES(EdgeList2CSR.out graph log)

ADD_EXECUTABLE(SortEdgeList.out sort_edgelist.cpp)
TARGET_LINK_LIBRARIES(SortEdgeList.out graph log)

ADD_EXECUTABLE(SortEdgeListWithWeight.out sort_edgelist_with_weight.cpp)
TARGET_LINK_LIBRARIES(SortEdgeListWithWeight.out graph log)

ADD_EXECUTABLE(SortDirectedEdgeList.out sort_directed_edgelist.cpp)
TARGET_LINK_LIBRARIES(SortDirectedEdgeList.out graph log)


ADD_EXECUTABLE(PartitionGraph.out partition_graph.cpp)
TARGET_LINK_LIBRARIES(PartitionGraph.out graph log)

ADD_EXECUTABLE(AssignLabel.out assign_label.cpp)

ADD_EXECUTABLE(AssignWeight.out assign_weight.cpp)
TARGET_LINK_LIBRARIES(AssignWeight.out graph log)

ADD_EXECUTABLE(AssignEdgeLabel.out assign_edge_label.cpp)
TARGET_LINK_LIBRARIES(AssignEdgeLabel.out graph log)
