#ifndef TRANSMISSION_PAIR_IMPL_HPP
#define TRANSMISSION_PAIR_IMPL_HPP

#include <regex>
#include <set>

#include <boost/math/tools/minima.hpp>

#include "load_input.hpp"

template <typename T>
transmission_pair_impl<T>::transmission_pair_impl() = default;

template <typename T>
transmission_pair_impl<T>::transmission_pair_impl(
	const transmission_pair_impl&)
	= default;

template <typename T>
transmission_pair_impl<T>::transmission_pair_impl(transmission_pair_impl&&) = default;

template <typename T>
transmission_pair_impl<T>& transmission_pair_impl<T>::
operator=(const transmission_pair_impl& other)
	= default;

template <typename T>
transmission_pair_impl<T>& transmission_pair_impl<T>::
operator=(transmission_pair_impl&& other)
	= default;

template <typename T>
void transmission_pair_impl<T>::init(std::istream& input_,
	const std::regex transmitter_regex_,
	const std::regex recipient_regex_,
	bool retain_unshared_)
{
	// 1. load data from stdin/file
	std::vector<T> transmitter_vec, recipient_vec;
	load_input(input_, transmitter_vec, recipient_vec, transmitter_regex_,
		recipient_regex_);

	// 2. sort vectors
	std::sort(transmitter_vec.begin(), transmitter_vec.end(), trait::compare);
	std::sort(recipient_vec.begin(), recipient_vec.end(), trait::compare);

	if (verbose)
	{
		std::cerr << "Input BEFORE initialization\n";
		std::cerr << "Transmitter:\n";
		for (const auto& i : transmitter_vec)
		{
			std::cerr << i;
		}
		std::cerr << "\nRecipient:\n";
		for (const auto& i : recipient_vec)
		{
			std::cerr << i;
		}
		std::cerr << '\n';
	}

	// 3. initialize traits/sequences
	init_(transmitter_vec, recipient_vec, retain_unshared_);

	if (verbose)
	{
		std::cerr << "Input AFTER initialization\n";
		std::cerr << "Transmitter:\n";
		for (const auto& i : m_transmitter_population)
		{
			std::cerr << i;
		}
		std::cerr << "\nRecipient:\n";
		for (const auto& i : m_recipient_population)
		{
			std::cerr << i;
		}
		std::cerr << '\n';
	}
}

/*
template <typename T>
void transmission_pair_impl<T>::set_avg_substitution_rate(double rate)
{
        substitution_model<T>::set_avg_substitution_rate(rate);
}
*/

template <typename T>
transmission_pair_impl<T>::transmission_pair_impl(
	const std::string& file_stem_parameters_, double rate)
	: substitution_model<T>(file_stem_parameters_, rate) {}

template <typename T>
double transmission_pair_impl<T>::pvalue(count_type num_trials, double time,
	double lambda, uint64_t seed) const
{
	uint32_t K = m_transmitter_population.size();
	double dist = this->weighted_distance(m_transmitter_population,
		m_recipient_population, time);
	std::default_random_engine rng(seed);
	return perform_test(K, m_transmitter_population.begin(),
		m_transmitter_population.m_total_counts,
		m_recipient_population.m_total_counts, dist, lambda,
		num_trials, rng, T::visitor);
}

template <typename T>
void transmission_pair_impl<T>::simulate_fixed(
	uint32_t K, count_type num_simulations, count_type num_trials,
	count_type transmitter_coverage, count_type recipient_coverage,
	double lambda, double p0, double f0, uint64_t seed) const
{
	const double p_i = (1 - p0) / (K - 1);

	double vec_p_T[K];
	vec_p_T[0] = p0;
	for (count_type i = 1; i < K; ++i)
	{
		vec_p_T[i] = p_i;
	}

	std::cerr << "Transmitter\n"
			  << std::fixed << std::setprecision(3);
	std::copy(vec_p_T, vec_p_T + K,
		std::ostream_iterator<double>(std::cerr, "  "));
	std::cerr << '\n';

	double vec_p_R[K];
	const double denom = f0 * p0 + (1 - p0);
	vec_p_R[0] = p0 * f0 / denom;
	for (count_type i = 1; i < K; ++i)
	{
		vec_p_R[i] = p_i / denom;
	}

	std::cerr << "Recipient\n"
			  << std::fixed << std::setprecision(3);
	std::copy(vec_p_R, vec_p_R + K,
		std::ostream_iterator<double>(std::cerr, "  "));
	std::cerr << '\n';

	const std::size_t num_elements_row_float = 1;
	std::unique_ptr<double[]> vec_float(
		new double[num_simulations * num_elements_row_float]);

	const std::size_t num_elements_row_sample = 3 * K + 1;
	std::unique_ptr<count_type[]> vec_sample(
		new count_type[num_simulations * num_elements_row_sample]);

#pragma omp parallel
	{
		int omp_seed = seed *
#ifdef _OPENMP
			(omp_get_thread_num() + 1)
#else
			1
#endif
			;
		std::default_random_engine rng(omp_seed);

		double dist;

#pragma omp for schedule(static, 10)
		for (uint32_t i = 0; i < num_simulations; ++i)
		{
#ifndef NDEBUG
#pragma omp critical
			{
				std::cerr << i << '\n';
			}
#endif

			double* const float_offset = vec_float.get() + i * num_elements_row_float;

			count_type* const sample_offset = vec_sample.get() + i * num_elements_row_sample;

			count_type* const bottleneck_size = sample_offset + (0);
			count_type* const vec_X_T = sample_offset + (1);
			count_type* const vec_X_I = sample_offset + (1 + K);
			count_type* const vec_X_R = sample_offset + (1 + 2 * K);

			// generate samples
			*bottleneck_size = generate_transmission_sample(
				K, vec_p_T, vec_p_R, vec_X_T, vec_X_I, vec_X_R, transmitter_coverage,
				recipient_coverage, lambda, rng);

			dist = distance(K, vec_X_T, vec_X_R, transmitter_coverage, recipient_coverage);
			float_offset[0] = perform_test(K, vec_X_T, transmitter_coverage, recipient_coverage, dist, lambda, num_trials, rng);
		}
	}

	std::stringstream output_file_name;
	output_file_name << std::fixed << std::setprecision(2) << "sim_l=" << lambda
					 << "_p0=" << p0 << "_f0=" << f0 << ".txt";

	std::ofstream output(output_file_name.str());
	output << "p\tB";

	const char field_sep = '\t';

	for (const auto& i : { "XT_", "XI_", "XR_" })
	{
		for (uint32_t j = 0; j < K; ++j)
		{
			output << field_sep << i << j;
		}
	}
	output << '\n';

	const double* float_offset = vec_float.get();
	const count_type* sample_offset = vec_sample.get();
	output << std::setprecision(4);
	for (uint32_t i = 0; i < num_simulations; float_offset += num_elements_row_float, sample_offset += num_elements_row_sample, ++i)
	{
		output << float_offset[0];
		for (auto it = sample_offset; it != sample_offset + num_elements_row_sample; ++it)
		{
			output << field_sep << (*it);
		}
		output << '\n';
	}
}

template <typename T>
void transmission_pair_impl<T>::simulate_variable(
	uint32_t K, count_type num_simulations, count_type num_trials,
	count_type transmitter_coverage, count_type recipient_coverage,
	double lambda, double p0, double f0, double alpha, uint64_t seed) const
{
	double vec_alpha[K];
	std::fill(vec_alpha, vec_alpha + K, alpha);

	std::cerr << "alpha vector for Dirichlet\n"
			  << std::fixed << std::setprecision(3);
	std::copy(vec_alpha, vec_alpha + K, std::ostream_iterator<double>(std::cerr, "  "));
	std::cerr << '\n';

	const std::size_t num_elements_row_float = 2 * K + 1;
	std::unique_ptr<double[]> vec_float(new double[num_simulations * num_elements_row_float]);

	const std::size_t num_elements_row_sample = 3 * K + 1;
	std::unique_ptr<count_type[]> vec_sample(new count_type[num_simulations * num_elements_row_sample]);

#pragma omp parallel
	{
		int omp_seed = seed *
#ifdef _OPENMP
			(omp_get_thread_num() + 1)
#else
			1
#endif
			;
		std::default_random_engine rng(omp_seed);

		double dist;

#pragma omp for schedule(static, 10)
		for (uint32_t i = 0; i < num_simulations; ++i)
		{
#ifndef NDEBUG
#pragma omp critical
			{
				std::cerr << i << '\n';
			}
#endif

			double* const float_offset = vec_float.get() + i * num_elements_row_float;

			// generate transmitter population
			double* const vec_p_T = float_offset + (1);
			gen_rand_dirichlet(vec_alpha, vec_alpha + K, vec_p_T, rng);

			// generate recipient population
			double* const vec_p_R = float_offset + (1 + K);

			const double p0 = vec_p_T[0];
			const double denom = f0 * p0 + (1 - p0);

			vec_p_R[0] = p0 * f0 / denom;
			for (std::size_t i = 1; i < K; ++i)
			{
				vec_p_R[i] = vec_p_T[i] / denom;
			}

			count_type* const sample_offset = vec_sample.get() + i * num_elements_row_sample;

			count_type* const bottleneck_size = sample_offset + (0);
			count_type* const vec_X_T = sample_offset + (1);
			count_type* const vec_X_I = sample_offset + (1 + K);
			count_type* const vec_X_R = sample_offset + (1 + 2 * K);

			// generate samples
			*bottleneck_size = generate_transmission_sample(
				K, vec_p_T, vec_p_R, vec_X_T, vec_X_I, vec_X_R, transmitter_coverage,
				recipient_coverage, lambda, rng);

			dist = distance(K, vec_X_T, vec_X_R, transmitter_coverage,
				recipient_coverage);
			float_offset[0] = perform_test(K, vec_X_T, transmitter_coverage, recipient_coverage,
				dist, lambda, num_trials, rng);
		}
	}

	std::stringstream output_file_name;
	output_file_name << std::fixed << std::setprecision(2) << "sim_l=" << lambda
					 << "_p0=" << p0 << "_f0=" << f0 << "_alpha=" << alpha
					 << ".txt";

	std::ofstream output(output_file_name.str());
	output << "p";

	const char field_sep = '\t';

	for (const auto& i : { "pT_", "pR_", "B", "XT_", "XI_", "XR_" })
	{
		for (uint32_t j = 0; j < K; ++j)
		{
			output << field_sep << i;

			if (i[0] == 'B')
			{
				break;
			}

			output << j;
		}
	}
	output << '\n';

	const double* float_offset = vec_float.get();
	const count_type* sample_offset = vec_sample.get();
	output << std::setprecision(4);
	for (uint32_t i = 0; i < num_simulations;
		 float_offset += num_elements_row_float,
				  sample_offset += num_elements_row_sample, ++i)
	{
		output << float_offset[0];
		for (auto it = float_offset + 1;
			 it != float_offset + num_elements_row_float; ++it)
		{
			output << field_sep << (*it);
		}
		for (auto it = sample_offset; it != sample_offset + num_elements_row_sample;
			 ++it)
		{
			output << field_sep << (*it);
		}
		output << '\n';
	}
}

/* private members */
template <typename T>
void transmission_pair_impl<T>::init_(std::vector<T>& transmitter_vec,
	std::vector<T>& recipient_vec,
	bool retain_unshared_)
{
	// 1. extract valid loci
	auto locus_extractor = [](const std::vector<T>& cont, std::set<char>& bases,
		std::size_t pos) -> bool {
		bases.clear();
		bool has_gap = false;

		for (const auto& i : cont)
		{
			if ((i.m_name[pos] == '-') || (i.m_name[pos] == '*'))
			{
				has_gap = true;
			}
			else
			{
				if (!T::is_valid_elem(i.m_name[pos]))
				{
					std::cerr << "Sequence '" << i.m_name << "' contains invalid char '"
							  << i.m_name[pos] << "' at position " << pos << "\n";
					std::cerr << "p: ";
					std::cout << "NA\n";
					exit(EXIT_SUCCESS);
				}
				bases.insert(i.m_name[pos]);
			}
		}
		return has_gap;
	};

	std::set<char> transmitter_bases;
	std::set<char> recipient_bases;
	std::set<char> common_bases;
	bool has_gap, has_shared;
	std::vector<std::size_t> include_loci;
	std::vector<std::size_t> gap_loci;
	std::vector<std::size_t> unshared_base_loci;

	const auto length = transmitter_vec.front().m_name.length();
	for (std::size_t j = 0; j < length; ++j)
	{
		has_gap = (locus_extractor(transmitter_vec, transmitter_bases, j) || locus_extractor(recipient_vec, recipient_bases, j));
		common_bases.clear();
		std::set_intersection(transmitter_bases.begin(), transmitter_bases.end(),
			recipient_bases.begin(), recipient_bases.end(),
			std::inserter(common_bases, common_bases.begin()));

		has_shared = !common_bases.empty();

		if ((!has_gap) && ((has_shared || retain_unshared_)))
		{
			include_loci.push_back(j);
		}
		else
		{
			if (has_gap)
			{
				gap_loci.push_back(j);
			}
			else
			{
				unshared_base_loci.push_back(j);
			}
		}
	}

	if (verbose)
	{
		typedef std::pair<std::string, std::vector<std::size_t>&> my_pair;
		for (const auto& i :
			{ my_pair("Retained", include_loci), my_pair("     Gap", gap_loci),
				my_pair("Unshared", unshared_base_loci) })
		{
			std::cerr << i.first << " loci: ";
			std::copy(i.second.begin(), i.second.end(),
				std::ostream_iterator<std::size_t>(std::cerr, " "));
			std::cerr << '\n';
		}
	}

	if (include_loci.size() == 0)
	{
		std::cerr
			<< "There are no loci that do not"
			<< (retain_unshared_ ? "" : " either") << " contain gaps"
			<< (retain_unshared_
					   ? ""
					   : " or have any shared bases between transmitter and recipient")
			<< ".\n";
		std::cerr << "p: ";
		std::cout << "NA\n";
		exit(EXIT_SUCCESS);
	}

	// 2. build final population
	m_transmitter_population.swap(transmitter_vec);
	m_recipient_population.swap(recipient_vec);

	auto concatenate_seq = [&include_loci](population<T>& pop) {
		std::map<std::string, T> pop_map;

		std::string temp;
		count_type sum = 0;
		for (auto& j : pop)
		{
			sum += j.m_count;

			temp.clear();
			for (const auto& k : include_loci)
			{
				temp.push_back(j.m_name[k]);
			}

			auto it = pop_map.insert(std::pair<std::string, T>(temp, j));
			if (it.second)
			{
				it.first->second.m_name = temp;
			}
			else
			{
				it.first->second.m_count += j.m_count;
			}
		}

		pop.clear();
		pop.reserve(pop_map.size());
		for (const auto& i : pop_map)
		{
			pop.push_back(i.second);
		}

		pop.m_total_counts = sum;
	};
	concatenate_seq(m_transmitter_population);
	concatenate_seq(m_recipient_population /*, m_transmitter_population*/);

	if (verbose)
	{
		std::cerr << "Transmitter Population: (Total: "
				  << m_transmitter_population.m_total_counts << ")\n";
		int j = 0;
		for (const auto& i : m_transmitter_population)
		{
			++j;
			std::cerr << i << '\n';
		}

		std::cerr << "\nRecipient Population: (Total: "
				  << m_recipient_population.m_total_counts << ")\n";
		j = 0;
		for (const auto& i : m_recipient_population)
		{
			++j;
			std::cerr << i << '\n';
		}
	}
}

template <>
void transmission_pair_impl<trait>::init_(std::vector<trait>& transmitter_vec,
	std::vector<trait>& recipient_vec,
	bool retain_unshared_)
{
	// 1. check that recipient traits form a subset of the transmitter traits
	if (!std::includes(transmitter_vec.cbegin(), transmitter_vec.cend(),
			recipient_vec.cbegin(), recipient_vec.cend()))
	{
		std::cerr
			<< "ERROR: Recipient traits are not a subset of transmitter traits!\n";
		std::cerr << "p: ";
		std::cout << "NA\n";
		exit(EXIT_SUCCESS);
	}

	// 2. build new recipient vector
	std::vector<trait> new_rec_vec;
	new_rec_vec.reserve(transmitter_vec.size());
	for (const auto& i : transmitter_vec)
	{
		const auto it = find_if(recipient_vec.cbegin(), recipient_vec.cend(),
			[&i](const trait& val) {
				return val.m_name == i.m_name;
			});
		if (it == recipient_vec.cend())
		{
			// not found
			new_rec_vec.emplace_back(i.m_name, 0);
		}
		else
		{
			// found
			new_rec_vec.push_back(*it);
		}
	}

	// 3. build final population
	m_transmitter_population.swap(transmitter_vec);
	m_recipient_population.swap(new_rec_vec);

	auto normalisation = [](population<trait>& pop) {
		count_type sum = 0;
		for (const auto& j : pop)
		{
			sum += j.m_count;
		}

		pop.m_total_counts = sum;
	};
	normalisation(m_transmitter_population);
	normalisation(m_recipient_population);
}

#endif /* TRANSMISSION_PAIR_IMPL_HPP */
