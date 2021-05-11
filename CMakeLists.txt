cmake_minimum_required(VERSION 3.16)
project(geoLib)

set(CMAKE_CXX_STANDARD 17)


find_package(Boost MODULE REQUIRED filesystem system thread)
find_package(GTest REQUIRED CONFIG)



add_executable(${PROJECT_NAME}
        geoLibMain.cpp
        geographicLibrary.h
        )
target_link_libraries(${PROJECT_NAME}             PRIVATE
        Boost::system
        )
target_compile_features(${PROJECT_NAME}           PRIVATE
        cxx_std_17
        )


add_executable(geoLib_UnitTests
        geoLibTests.cpp
        geographicLibrary.h
        )
target_include_directories(geoLib_UnitTests      PRIVATE
        ${Boost_INCLUDE_DIRS}
        GoogleTest
        )
target_link_libraries(geoLib_UnitTests           PRIVATE
        Boost::system
        GTest::gtest
        GTest::gtest_main
        )
target_compile_features(geoLib_UnitTests         PRIVATE
        cxx_std_17
        )