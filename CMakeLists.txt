cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)


include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)

add_subdirectory (mechsystem)

add_executable (test_ode demos/test_ode.cpp)
target_link_libraries (test_ode PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

add_executable (demo_legendre demos/demo_legendre.cpp)
target_link_libraries (demo_legendre PUBLIC nanoblas)

add_executable (test_explicit_rk demos/test_explicit_rk.cpp)
target_link_libraries (test_explicit_rk PUBLIC nanoblas)

add_library(ASC_ode INTERFACE)
target_include_directories(ASC_ode INTERFACE src)

add_executable(test_pendulum demos/test_pendulum.cpp)
target_link_libraries(test_pendulum PRIVATE ASC_ode nanoblas)