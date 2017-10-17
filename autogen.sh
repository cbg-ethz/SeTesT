#!/bin/bash

clean_files() {
	echo "Cleaning bootstrapped files"

	rm -rf .Rapp.history
	rm -rf .RData
	rm -rf .Rhistory

	rm -rf .deps/
	rm -rf .Tpo
	rm -rf *.dSYM
	rm -rf setest
	rm -rf setest-*.tar.bz2
	rm -rf setest-*/

	# Autotools
	rm -rf Makefile
	rm -rf Makefile.in
	rm -rf aclocal.m4
	rm -rf autom4te.cache/
	rm -rf compile
	rm -rf config.h*
	rm -rf configure
	rm -rf depcomp
	rm -rf install-sh
	rm -rf missing
	rm -rf stamp-h1
	rm -rf test-driver
	rm -rf *.log

	# source
	rm -rf src/*.o
	rm -rf src/.deps
	rm -rf src/.dirstamp

	# testsuite
	rm -rf testsuite/*.log
	rm -rf testsuite/*.trs
	rm -rf testsuite/sequence_compile_test
	rm -rf testsuite/substitution_model_test
	rm -rf testsuite/substitution_model_sequence_space_test
	rm -rf testsuite/*.o
	rm -rf testsuite/.deps
	rm -rf testsuite/.dirstamp

	# OS X cruft
	find . -name '.DS_Store' -type f -delete
}

if [[ "$1" == "--clean" ]]
then
	clean_files
	exit
fi

echo "Bootstrapping Autotools"
autoreconf -vif

if [[ "$1" == "--test" ]]
then
	echo "${DISTCHECK_CONFIGURE_FLAGS}"
	./configure ${DISTCHECK_CONFIGURE_FLAGS}
	DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS}" make distcheck
	clean_files
fi
