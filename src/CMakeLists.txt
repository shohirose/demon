cmake_minimum_required(VERSION 3.8)
project(demon2 CXX)

add_executable(1box
       	main.cpp
	functions.cpp
	particle.cpp
	)
target_compile_options(1box
	PUBLIC
	-Wall
	-Wextra
	-Wpedantic
	)
target_include_directories(1box
	PUBLIC
	${PROJECT_SOURCE_DIR}/../include
	)
target_compile_features(1box
	PUBLIC
	cxx_std_14
	)

