cmake_minimum_required(VERSION 2.6)
project(rxpathtracer_gccO2)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg") #debug out
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") # -Wno-deprecated-declarations -Wno-write-strings")  #"-pg"
#gprof rxpathtracer_gccO2 gmon.out >rxpathtracer_gccO2.11.gmon.out.txt

find_package(MPI REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${MPI_INCLUDE_PATH})
add_executable(rxpathtracer_gccO2
			../src/RXPathTracer_Io.cpp
			../src/RXPathTracer_Math.cpp
			../src/RXPathTracer_Random.cpp
			../src/RXPathTracer_GrainPrinter.cpp
			../src/RXPathTracer_NucleationPrinter.cpp
			../src/RXPathTracer_TimeLogger.cpp
			../src/RXPathTracer_PhysicalConstants.cpp
			../src/RXPathTracer_Main.cpp
			../src/RXPathTracer_Kernel.cpp
)
target_link_libraries(rxpathtracer_gccO2  ${CMAKE_SOURCE_DIR}/src/vorosrc/src/libvoro++.a)
target_link_libraries(rxpathtracer_gccO2 ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(rxpathtracer_gccO2 PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(rxpathtracer_gccO2 PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if(VTK_LIBRARIES)
  target_link_libraries(rxpathtracer_gccO2 ${VTK_LIBRARIES})
else()
  target_link_libraries(rxpathtracer_gccO2 vtkHybrid vtkWidgets)
endif()
