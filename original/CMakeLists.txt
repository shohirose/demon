cmake_minimum_required(VERSION 3.8)
project(demon C)

add_executable(1box 1box.c)
target_compile_options(1box
	PUBLIC
	-Wall
	-Wextra
	-Wpedantic
	)

add_executable(phase_omp phase_omp.c)
target_compile_options(phase_omp
	PUBLIC
	-Wall
	-Wextra
	-Wpedantic
	)

