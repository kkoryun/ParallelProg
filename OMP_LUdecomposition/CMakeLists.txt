cmake_minimum_required (VERSION 2.8.11)
project (ChM)

set(SRC LU.cpp Gradient.cpp)

add_executable(${PROJECT_NAME} ${SRC})