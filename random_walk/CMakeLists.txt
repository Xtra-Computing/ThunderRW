SET(SFMT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/dependency/SFMT-src-1.5.1")
INCLUDE_DIRECTORIES(${SFMT_INCLUDE_DIR})

SET(CREEPER_SOURCE_FILES
        walker_generator.cpp ${SFMT_INCLUDE_DIR}/SFMT.c)


ADD_EXECUTABLE(deepwalk.out ${CREEPER_SOURCE_FILES} apps/deepwalk.cpp)
TARGET_LINK_LIBRARIES(deepwalk.out graph log io)
TARGET_COMPILE_DEFINITIONS(deepwalk.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(ppr.out ${CREEPER_SOURCE_FILES} apps/ppr.cpp)
TARGET_LINK_LIBRARIES(ppr.out graph log io)
TARGET_COMPILE_DEFINITIONS(ppr.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(metapath.out ${CREEPER_SOURCE_FILES} apps/metapath.cpp)
TARGET_LINK_LIBRARIES(metapath.out graph log io)
TARGET_COMPILE_DEFINITIONS(metapath.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(node2vec.out ${CREEPER_SOURCE_FILES} apps/node2vec.cpp)
TARGET_LINK_LIBRARIES(node2vec.out graph log io)
TARGET_COMPILE_DEFINITIONS(node2vec.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(genericwalk.out ${CREEPER_SOURCE_FILES} apps/genericwalk.cpp)
TARGET_LINK_LIBRARIES(genericwalk.out graph log io)
TARGET_COMPILE_DEFINITIONS(genericwalk.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(its_preprocess.out its_preprocess.cpp)
TARGET_LINK_LIBRARIES(its_preprocess.out graph log io)

ADD_EXECUTABLE(alias_preprocess.out alias_preprocess.cpp)
TARGET_LINK_LIBRARIES(alias_preprocess.out graph log io)

ADD_EXECUTABLE(rj_preprocess.out rejection_preprocess.cpp)
TARGET_LINK_LIBRARIES(rj_preprocess.out graph log io)

ADD_EXECUTABLE(graphwalker.out home_grown_graph_walker.cpp ${SFMT_INCLUDE_DIR}/SFMT.c)
TARGET_LINK_LIBRARIES(graphwalker.out graph log io)
TARGET_COMPILE_DEFINITIONS(graphwalker.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

ADD_EXECUTABLE(evaluation_impact.out impact_on_existing_systems.cpp ${SFMT_INCLUDE_DIR}/SFMT.c)
TARGET_LINK_LIBRARIES(evaluation_impact.out graph log io)
TARGET_COMPILE_DEFINITIONS(evaluation_impact.out PRIVATE SFMT_MEXP=19937 HAVE_SSE2=1)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
