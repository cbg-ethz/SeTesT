#include <config.h>

#include <iostream>
#include <string>
#include <thread>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#define MAGIC_NUMBER 0xDEADBEEF

bool verbose = false;

enum class enum_data_type
{
	trait,
	protein,
	dna
};
enum class enum_run_mode
{
	simulate,
	test
};

#include "substitution_model.hpp"
#include "population.hpp"
#include "types.hpp"
#include "transmission_pair.hpp"

int main(int argc, char* argv[])
{
	std::cerr.imbue(std::locale("en_US.UTF-8"));

	// General configuration
	double lambda;
	count_type num_trials;
	double time;
	std::string regex_transmitter_str;
	std::string regex_recipient_str;
	unsigned num_threads;
	std::string parameter_file_name;
	uint64_t random_seed;

	// Simulation parameters
	uint32_t K;
	double alpha;
	count_type num_simulations;
	double p0;
	double f0;
	count_type transmitter_coverage;
	count_type recipient_coverage;

	/* 1) set up program options */
	// program options
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "Print this help")
		("verbose,v", "Show progress indicator while running the Monte Carlo simulations");

	// configuration options
	boost::program_options::options_description config("Configuration");
	config.add_options()
		("seed", boost::program_options::value<decltype(random_seed)>(&random_seed)->default_value(42), "Value of seed for deterministic run. A value of 0 will pick a random seed from some non-deterministic entropy source")
		(",l", boost::program_options::value<decltype(lambda)>(&lambda)->default_value(0.515, "0.515"), "Lambda bottleneck parameter")
		(",N", boost::program_options::value<decltype(num_trials)>(&num_trials)->default_value(1E6, "1,000,000"), "Number of Monte Carlo samples to calculate approximate p-value")
		("time", boost::program_options::value<decltype(time)>(&time)->default_value(100, "100"), "Time (in days) between (suspected) transmission event and sample date")
		(",T", boost::program_options::value<decltype(regex_transmitter_str)>(&regex_transmitter_str)->default_value("T/"), "Regex for the transmitter sequences")
		(",R", boost::program_options::value<decltype(regex_recipient_str)>(&regex_recipient_str)->default_value("R/"), "Regex for the recipient sequences")
		(",t", boost::program_options::value<decltype(num_threads)>(&num_threads)->default_value(std::thread::hardware_concurrency()), "Number of threads to use for simulations. Defaults to number of logical cores found")
		("protein", boost::program_options::value<decltype(parameter_file_name)>(&parameter_file_name)->implicit_value(""), "Load input as protein sequences.")
		("dna", boost::program_options::value<decltype(parameter_file_name)>(&parameter_file_name)->implicit_value(""), "Load input as peptide/protein sequences. By default " PACKAGE_NAME " assumes input to be a trait.");

	boost::program_options::options_description simulation("Simulation");
	simulation.add_options()
		(",s", "Runs simulations with fixed transmitter population")
		(",M", boost::program_options::value<decltype(alpha)>(&alpha)->implicit_value(1), "Runs simulations with Dir(alpha, ..., alpha) distribution of transmitter population")
		(",K", boost::program_options::value<decltype(K)>(&K)->default_value(2), "Number of traits in transmitter")
		(",S", boost::program_options::value<decltype(num_simulations)>(&num_simulations)->default_value(10000, "10,000"), "Number of meta-simulations to perform")
		(",p", boost::program_options::value<decltype(p0)>(&p0)->default_value(0, "1/K"), "Fraction in transmitter population of first trait")
		(",f", boost::program_options::value<decltype(f0)>(&f0)->default_value(1), "Transmissibility of first trait, where all other traits have 1. The null hypothesis H_0 is given with -f 1")
		("reads-transmitter", boost::program_options::value<decltype(transmitter_coverage)>(&transmitter_coverage)->default_value(10000, "10,000"), "Simulated number of total number of reads in transmitter")
		("reads-recipient", boost::program_options::value<decltype(recipient_coverage)>(&recipient_coverage)->default_value(10000, "10,000"), "Simulated number of total number of reads in recipient");

	// hidden options, i.e., input files
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("input-files", boost::program_options::value<std::vector<std::string>>(), "input files");

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(simulation).add(hidden);
	boost::program_options::options_description visible("Allowed options for " PACKAGE_NAME " " PACKAGE_VERSION);
	visible.add(generic).add(config).add(simulation);
	boost::program_options::positional_options_description p;
	p.add("input-files", -1);
	boost::program_options::variables_map global_options;

	/* 2) parse program options */
	try
	{
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), global_options);

		// show help options
		if (global_options.count("help"))
		{
			std::cerr << visible << '\n';
			exit(EXIT_SUCCESS);
		}

		boost::program_options::notify(global_options);
	}
	catch (boost::program_options::required_option& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		exit(EXIT_FAILURE);
	}
	catch (boost::program_options::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		exit(EXIT_FAILURE);
	}

#ifdef _OPENMP
	omp_set_num_threads(num_threads);
#endif

	std::vector<std::string> input_files;
	verbose = global_options.count("verbose");
	if (global_options.count("input-files"))
	{
		input_files = global_options["input-files"].as<std::vector<std::string>>();
	}

	// determine run mode
	enum_run_mode run_mode = enum_run_mode::test;
	if (global_options.count("-s") + global_options.count("-M") > 1)
	{
		std::cerr << "You have cannot have -s and -M enabled.\n";
		exit(EXIT_FAILURE);
	}
	else
	{
		if ((global_options.count("-s")) || (global_options.count("-M")))
		{
			run_mode = enum_run_mode::simulate;
		}
	}

	// determine data-type
	enum_data_type data_type;
	if (global_options.count("protein") + global_options.count("dna") > 1)
	{
		std::cerr << "You have cannot specify both --protein and --dna.\n";
		exit(EXIT_FAILURE);
	}

	if (global_options.count("protein"))
	{
		data_type = enum_data_type::protein;
	}
	else if (global_options.count("dna"))
	{
		data_type = enum_data_type::dna;
	}
	else
	{
		data_type = enum_data_type::trait;
	}

	if (run_mode == enum_run_mode::simulate)
	{
		if (K < 2)
		{
			std::cerr << "K has to be at least 2!\n";
			exit(EXIT_FAILURE);
		}

		if (lambda <= 0)
		{
			std::cerr << "lambda has to be strictly larger than 0!\n";
			exit(EXIT_FAILURE);
		}

		if (p0 == 0)
		{
			p0 = 1.0 / K;
		}

		if (p0 < 0)
		{
			std::cerr << "p0 has to be strictly positive!\n";
			exit(EXIT_FAILURE);
		}

		if (p0 >= 1)
		{
			std::cerr << "p0 has to be strictly smaller than 1!\n";
			exit(EXIT_FAILURE);
		}
	}

	/* 3) show parameters */
	if (verbose)
	{
		std::cerr << "Running in ";
		switch (data_type)
		{
			case enum_data_type::protein:
				std::cerr << "protein";
				break;

			case enum_data_type::dna:
				std::cerr << "DNA";
				break;

			case enum_data_type::trait:
				std::cerr << "trait";
				break;
		}
		std::cerr << " mode\n";

		std::cerr << std::setprecision(4)
				  << "  General parameters\n"
				  << "\tLambda:                  " << lambda << " (mu = " << lambda_to_mean(lambda) << ")\n"
				  << "\tNumber of trials:        " << num_trials << '\n'
				  << "\tTime transmit -> sample: " << time << " days\n"
				  << "\tTransmitter regex:       " << regex_transmitter_str << '\n'
				  << "\tRecipient regex:         " << regex_recipient_str << '\n'
				  << "\tNumber of threads:       " << num_threads << "\n\n";

		if (run_mode == enum_run_mode::simulate)
		{
			std::cerr
				<< "  Simulation parameters\n"
				<< "\tK:                       " << K << '\n'
				<< "\tNumber of simulations:   " << num_simulations << '\n'
				<< "\tp0:                      " << p0 << '\n'
				<< "\tf0:                      " << f0 << '\n'
				<< "\tReads in Transmitter:    " << transmitter_coverage << '\n'
				<< "\tReads in Recipient:      " << recipient_coverage << "\n\n";
		}
	}

	/* 4) start actual program */
	std::unique_ptr<transmission_pair> trans_pair;
	switch (data_type)
	{
		case enum_data_type::trait:
			// trait
			trans_pair = std::unique_ptr<transmission_pair>(new transmission_pair_impl<trait>);
			break;
		case enum_data_type::dna:
			// dna
			trans_pair = std::unique_ptr<transmission_pair>(new transmission_pair_impl<dna_sequence>(parameter_file_name, 2.75E-6));
			break;
		case enum_data_type::protein:
			// protein
			trans_pair = std::unique_ptr<transmission_pair>(new transmission_pair_impl<protein_sequence>(parameter_file_name, 4.6E-5));
			break;
	}

	const std::regex transmitter_regex(regex_transmitter_str);
	const std::regex recipient_regex(regex_recipient_str);
	switch (run_mode)
	{
		case enum_run_mode::simulate:
			std::cerr << "Running ";
			if (global_options.count("-s"))
			{
				std::cerr << "fixed-population simulations\n";
				trans_pair->simulate_fixed(K, num_simulations, num_trials, transmitter_coverage, recipient_coverage, lambda, p0, f0, random_seed);
			}
			else
			{
				std::cerr << "variable-population (alpha = " << alpha << ") simulations\n";
				trans_pair->simulate_variable(K, num_simulations, num_trials, transmitter_coverage, recipient_coverage, lambda, p0, f0, alpha, random_seed);
			}
			break;

		case enum_run_mode::test:
			std::cerr << "Running statistical test\n";
			if (input_files.size())
			{
				// load files from vector
				for (const auto& i : input_files)
				{
					if (!boost::filesystem::exists(i))
					{
						std::cerr << "File '" << i << "' does not exist! Exiting...\n";
						exit(EXIT_FAILURE);
					}

					std::ifstream input(i);
					trans_pair->init(input, transmitter_regex, recipient_regex, true);
					input.close();

					const double pvalue = trans_pair->pvalue(num_trials, time, lambda, random_seed);
					std::cerr << i << " - p: ";
					std::cout << pvalue << '\n';
				}
			}
			else
			{
				trans_pair->init(std::cin, transmitter_regex, recipient_regex, true);
				const double pvalue = trans_pair->pvalue(num_trials, time, lambda, random_seed);
				std::cerr << "p: ";
				std::cout << pvalue << '\n';
			}
			break;
	}

	exit(EXIT_SUCCESS);
}