cmake_minimum_required(VERSION 3.17)
project(kerr2Kitaev2)

set(CMAKE_CXX_STANDARD 17)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BOOST_ROOT /usr/local/include/boost)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-O1")
find_package(Boost 1.45.0)
find_package (Threads)
include_directories("/usr/include/eigen-3.4.0")

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(kerr2Kitaev2 main.cpp consts.hpp solver.hpp solver.cpp)
    target_link_libraries(kerr2Kitaev2 ${Boost_LIBRARIES})
    target_link_libraries(kerr2Kitaev2 Threads::Threads)

endif()
