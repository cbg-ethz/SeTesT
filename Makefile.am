AM_CPPFLAGS = -I$(top_srcdir)/src $(BOOST_CPPFLAGS)

bin_PROGRAMS = setest

setest_SOURCES = \
	src/defs.hpp \
	src/load_input.hpp \
	src/population.hpp \
	src/population_impl.hpp \
	src/sampling_impl.hpp \
	src/substitution_model.hpp \
	src/substitution_model_impl.hpp \
	src/transmission_pair.hpp \
	src/transmission_pair_impl.hpp \
	src/types.hpp \
	src/types_impl.hpp \
	src/utility_functions.hpp \
	src/main.cpp

setest_LDADD = \
	$(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)

dist_pkgdata_DATA = \
	data/amino_acid.distribution \
	data/amino_acid.matrix \
	data/dna.distribution \
	data/dna.matrix \
	examples/dna_faulty.fasta \
	examples/dna_nonsignificant.fasta \
	examples/dna_significant.fasta \
	examples/general_trait_faulty.tsv \
	examples/general_trait_nonsignificant.tsv \
	examples/general_trait_significant.tsv \
	examples/protein_faulty.fasta \
	examples/protein_nonsignificant.fasta \
	examples/protein_significant.fasta

# CHECKS
TESTS = $(check_PROGRAMS)
check_PROGRAMS = \
	testsuite/sequence_compile_test \
	testsuite/substitution_model_test \
	testsuite/substitution_model_sequence_space_test

testsuite_sequence_compile_test_SOURCES = testsuite/sequence_compile_test.cpp

testsuite_substitution_model_test_SOURCES = testsuite/substitution_model_test.cpp
testsuite_substitution_model_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)

testsuite_substitution_model_sequence_space_test_SOURCES = testsuite/substitution_model_sequence_space_test.cpp
testsuite_substitution_model_sequence_space_test_LDADD = $(BOOST_LDFLAGS) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
