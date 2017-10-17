##################
# Initialisation #
##################
# CMake requirements
# 3.1 is required for CMAKE_CXX_STANDARD
cmake_minimum_required(VERSION 3.1)

# set project name
cmake_policy(SET CMP0048 NEW)
project(setest VERSION 0.1 LANGUAGES CXX)

# on macOS, MACOSX_RPATH is enabled by default on more recent versions
# of CMake. Disable this behaviour, and let user enable it if need be.
cmake_policy(SET CMP0042 OLD)

# Set Release type for builds where CMAKE_BUILD_TYPE is unset
# This is usually a good default as this implictly enables
#
#   CXXFLAGS = -O3 -DNDEBUG
#
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Adhere to GNU filesystem layout conventions
include(GNUInstallDirs)

# Summarise found features at end
include(FeatureSummary)

################
# Dependencies #
################
# Request C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	message("OpenMP could not be found. ${PROJECT_NAME} will still work, albeit in single-threaded mode")
endif()

# Boost
find_package(Boost 1.59.0 COMPONENTS filesystem program_options system REQUIRED)

# config macros
configure_file(config.cmake.in config.h @ONLY)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")


#########
# Files #
#########
add_executable(
	setest
		src/defs.hpp
		src/load_input.hpp
		src/population.hpp
		src/population_impl.hpp
		src/sampling_impl.hpp
		src/substitution_model.hpp
		src/substitution_model_impl.hpp
		src/transmission_pair.hpp
		src/transmission_pair_impl.hpp
		src/types.hpp
		src/types_impl.hpp
		src/utility_functions.hpp
		src/main.cpp
)

target_include_directories(setest PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(setest ${Boost_LIBRARIES})

install(
	TARGETS setest
	DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(
	FILES README.md
	DESTINATION "${CMAKE_INSTALL_DOCDIR}")

install(FILES
	data/amino_acid.distribution
	data/amino_acid.matrix
	data/dna.distribution
	data/dna.matrix
	examples/dna_faulty.fasta
	examples/dna_nonsignificant.fasta
	examples/dna_significant.fasta
	examples/general_trait_faulty.tsv
	examples/general_trait_nonsignificant.tsv
	examples/general_trait_significant.tsv
	examples/protein_faulty.fasta
	examples/protein_nonsignificant.fasta
	examples/protein_significant.fasta
	DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")

#############
# Testsuite #
#############
# default disable testsuite
option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

if(BUILD_TESTING)
	# 1. copy over required data files
	# The source directory should _never_ be touched
	foreach(copyfile amino_acid.distribution amino_acid.matrix dna.distribution dna.matrix)
		configure_file(
			"${CMAKE_CURRENT_SOURCE_DIR}/data/${copyfile}"
			"${CMAKE_CURRENT_BINARY_DIR}/data/${copyfile}"
			COPYONLY)
	endforeach()

	# 2. instantiate tests
	foreach(test_var sequence_compile_test substitution_model_test substitution_model_sequence_space_test)
		add_executable(${test_var} testsuite/${test_var}.cpp)
		target_include_directories(${test_var} PRIVATE ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/src)
		target_link_libraries(${test_var} ${Boost_LIBRARIES})
		add_test(NAME ${test_var} COMMAND ${test_var})
	endforeach()
endif()

###########
# Summary #
###########
feature_summary(WHAT ALL)
