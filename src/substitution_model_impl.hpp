#ifndef SUBSTITUTION_MODEL_IMPL_HPP
#define SUBSTITUTION_MODEL_IMPL_HPP

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include <algorithm>
#include <numeric>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/odeint.hpp>

template <std::size_t N>
inline double calculate_average_substitution_rate(const double(&sub_mat)[N][N],
	const double(&prob_vec)[N])
{
	double rates[N][N];
	double sum = 0, norm;
	for (std::size_t j = 0; j < N; ++j)
	{
		norm = 0;

		for (std::size_t i = 0; i < N; ++i)
		{
			rates[i][j] = prob_vec[i] * sub_mat[i][j];

			if (i != j)
			{
				norm += rates[i][j];
			}
		}

		rates[j][j] = -norm;
		sum += prob_vec[j] * norm;
	}

#ifndef NDEBUG
	std::cerr << "Sum of columns:\n";
	for (std::size_t j = 0; j < N; ++j)
	{
		norm = 0;

		for (std::size_t i = 0; i < N; ++i)
		{
			norm += rates[i][j];
		}

		std::cerr << j << '\t' << norm << '\n';
	}
#endif

	return sum;
}

// ctor
template <typename T, bool log_space>
void substitution_model<T, log_space>::init(const std::string& file_stem)
{
	const std::string rate_matrix_file_name = file_stem + ".matrix";
	std::cout << rate_matrix_file_name << '\n';
	if (!boost::filesystem::exists(rate_matrix_file_name))
	{
		std::cerr << "Cannot open rate matrix file!\n" << std::endl;
		exit(EXIT_FAILURE);
	}

	const std::string stationary_distribution_file_name = file_stem + ".distribution";
	if (!boost::filesystem::exists(stationary_distribution_file_name))
	{
		std::cerr << "Cannot open stationary distribution file!\n" << std::endl;
		exit(EXIT_FAILURE);
	}

	auto file_loader = [](const std::string& input_file, double matrix[][N],
		const std::size_t matrix_rows)
	{
		std::string temp;
		std::ifstream input(input_file);
		std::size_t real_row = 0, line_num = 0;

		std::vector<std::string> split_vec;

		if (input.is_open())
		{
			while ((input.good()) && (real_row < matrix_rows))
			{
				++line_num;
				std::getline(input, temp);
				temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

				boost::trim_if(temp, boost::is_any_of("\t "));
				boost::split(split_vec, temp, boost::is_any_of("\t "),
					boost::token_compress_on);
				const std::size_t num_columns = split_vec.size();

				if ((num_columns == 1) && (!split_vec.front().length()))
				{
					// empty line
					continue;
				}

				if (split_vec.front()[0] == '#')
				{
					// comment line
					continue;
				}

				if (num_columns != N)
				{
					std::cerr << "In file '" << input_file << "', line number "
							  << line_num << " contains " << num_columns
							  << " columns, expected " << N << "!\n";
					exit(EXIT_FAILURE);
				}

				for (std::size_t i = 0; i < N; ++i)
				{
					try
					{
						matrix[real_row][i] = boost::lexical_cast<double>(split_vec[i]);
					}
					catch (boost::bad_lexical_cast&)
					{
						std::cerr << "In file '" << input_file << "', value '"
								  << split_vec[i] << "' in row " << line_num << ", column "
								  << i + 1 << " is invalid!\n";
						exit(EXIT_FAILURE);
					}
				}

				++real_row;
			}

			if (real_row != matrix_rows)
			{
				std::size_t missing_rows = matrix_rows - real_row;
				std::cerr << "In file '" << input_file << "', missing " << missing_rows
						  << " row" << (missing_rows != 1 ? "s" : "") << "!\n";
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			std::cerr << "Input file '" << input_file << "' could not be opened!\n";
			exit(EXIT_FAILURE);
		}
		input.close();

		// normalise rate matrix
		if (matrix_rows == 1)
		{
			// vector
			const double denom = 1.0 / std::accumulate(matrix[0], matrix[0] + N, 0.0);
			for (std::size_t i = 0; i < N; ++i)
			{
				matrix[0][i] *= denom;
			}
		}
		else
		{
			// matrix
			double mean;
			for (std::size_t i = 0; i < N; ++i)
			{
				// set diagonals to 0
				matrix[i][i] = 0;

				for (std::size_t j = i + 1; j < N; ++j)
				{
					mean = (matrix[i][j] + matrix[j][i]) / 2;

					// make symmetric
					matrix[i][j] = mean;
					matrix[j][i] = mean;
				}
			}
		}
	};

	// 1. load rate matrix
	file_loader(rate_matrix_file_name, m_data_matrix, N);

	// 2. load equilibrium distribution
	file_loader(stationary_distribution_file_name, &m_data_vector, 1);
	std::cout << "Sum of probability vector: "
			  << std::accumulate(m_data_vector, m_data_vector + N, 0.0) << '\n';

	auto display_matrix = [](const std::string& input_file,
		const double matrix[][N],
		const std::size_t matrix_rows)
	{
		std::cout << (matrix_rows == 1 ? "Vector" : "Matrix") << " from '"
				  << input_file << "':\n" << std::fixed << std::setprecision(4);
		for (std::size_t i = 0; i < matrix_rows; ++i)
		{
			for (std::size_t j = 0; j < N; ++j)
			{
				std::cout << std::right << std::setw(8) << matrix[i][j];
			}
			std::cout << '\n';
		}
	};

	// 3. Normalise complete rate matrix
	double sub_rate = calculate_average_substitution_rate(m_data_matrix, m_data_vector);
	std::cout << "Average substitution rate before normalisation: " << sub_rate
			  << '\n';
	for (double* i = &m_data_matrix[0][0]; i < &m_data_matrix[0][0] + N * N;
		 ++i)
	{
		(*i) /= sub_rate;
	}

	sub_rate = calculate_average_substitution_rate(m_data_matrix, m_data_vector);
	std::cout << "Average substitution rate after normalisation: " << sub_rate
			  << '\n';

	display_matrix(rate_matrix_file_name, m_data_matrix, N);
	display_matrix(stationary_distribution_file_name, &m_data_vector, 1);
}

template <typename T, bool log_space>
substitution_model<T, log_space>::substitution_model(std::string file_stem,
	double rate)
	: m_rate(rate), m_time(-1), m_cached(false)
{
	if (file_stem.empty())
	{
		constexpr std::array<const char*, 2> lookup_directories{
			{ "data", DATA_PATH }
		};

		constexpr const char* prefix = (std::is_same<T, protein_sequence>::value ? "amino_acid" : "dna");

		bool found_file = false;
		for (const auto& i : lookup_directories)
		{
			file_stem = i;
			file_stem.push_back('/');
			file_stem.append(prefix);
			found_file = boost::filesystem::exists(file_stem + ".matrix");
#ifndef NDEBUG
			std::cout << "Looking for parameter files in " << i << "/... "
					  << (found_file ? "yes" : "no") << '\n';
#endif
			if (found_file)
			{
				break;
			}
		}

		if (found_file == false)
		{
			std::cerr << "Could not locate parameter " << prefix << ".matrix!"
					  << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	std::cout << "Average substitution rate = " << m_rate << '\n';

	init(file_stem);
}

// set average base substitution rate
template <typename T, bool log_space>
void substitution_model<T, log_space>::set_avg_substitution_rate(double rate)
{
	m_rate = rate;
	m_cached = false;
}

// transition probabilities
template <typename T, bool log_space>
inline double substitution_model<T, log_space>::logP_from_to_seq(
	const std::string& from_sequence, const std::string& to_sequence,
	double time) const
{
	if (from_sequence.length() != to_sequence.length())
	{
		throw std::invalid_argument(
			std::string("Both strings have to be of equal length! Exiting...\n"));
	}

	if (log_space == true)
	{
		// log base
		const auto from_cend = from_sequence.cend();
		double result = 0;
		for (auto from_it = from_sequence.cbegin(), to_it = to_sequence.cbegin();
			 from_it != from_cend; ++from_it, ++to_it)
		{
			result += logP_from_to_base(*from_it, *to_it, time);
		}

		return result;
	}
	else
	{
		// non-log base
		const double p = P_from_to_seq(from_sequence, to_sequence, time);
		return (p ? std::log(p) : log_lower_bound);
	}
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::logP_from_to_seq(
	const std::string& from_sequence, const std::string& to_sequence) const
{
	return logP_from_to_seq(from_sequence, to_sequence, m_time);
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::P_from_to_seq(
	const std::string& from_sequence, const std::string& to_sequence,
	double time) const
{
	if (from_sequence.length() != to_sequence.length())
	{
		throw std::invalid_argument(
			std::string("Both strings have to be of equal length! Exiting...\n"));
	}

	if (log_space == true)
	{
		// log base
		return std::exp(logP_from_to_seq(from_sequence, to_sequence, time));
	}
	else
	{
		// non-log base
		const auto from_cend = from_sequence.cend();
		double result = 0;
		for (auto from_it = from_sequence.cbegin(), to_it = to_sequence.cbegin();
			 from_it != from_cend; ++from_it, ++to_it)
		{
			result *= P_from_to_base(*from_it, *to_it, time);
		}

		return result;
	}
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::P_from_to_seq(
	const std::string& from_sequence, const std::string& to_sequence) const
{
	return P_from_to_seq(from_sequence, to_sequence, m_time);
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::logP_from_to_base(
	char from_base, char to_base, double time) const
{
	if (log_space == true)
	{
		// log base
		return get_entry_impl(from_base, to_base, time);
	}
	else
	{
		// non-log base
		const double ij = P_from_to_base(from_base, to_base, time);
		return (ij ? std::log(ij) : log_lower_bound);
	}
}

template <typename T, bool log_space>
inline double
substitution_model<T, log_space>::logP_from_to_base(char from_base,
	char to_base) const
{
	return logP_from_to_base(from_base, to_base, m_time);
}

template <typename T, bool log_space>
inline double
substitution_model<T, log_space>::P_from_to_base(char from_base, char to_base,
	double time) const
{
	if (log_space == true)
	{
		// log base
		return std::exp(logP_from_to_base(from_base, to_base, time));
	}
	else
	{
		// non-log base
		return get_entry_impl(from_base, to_base, time);
	}
}

template <typename T, bool log_space>
inline double
substitution_model<T, log_space>::P_from_to_base(char from_base,
	char to_base) const
{
	return P_from_to_base(from_base, to_base, m_time);
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::logP_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence,
	double time) const
{
	if (from_sequence.length() != to_sequence.length())
	{
		throw std::invalid_argument(
			std::string("Both strings have to be of equal length! Exiting...\n"));
	}

	const auto from_cend = from_sequence.cend();
	double result = 0;
	for (auto from_it = from_sequence.cbegin(), to_it = to_sequence.cbegin();
		 from_it != from_cend; ++from_it, ++to_it)
	{
		result += logP_from_to_base(*from_it, *to_it, time) + logP_marginal_base(*from_it) - logP_marginal_base(*to_it);
	}
	return result;
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::logP_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence) const
{
	return logP_of_parent_given_child(from_sequence, to_sequence, m_time);
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::P_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence,
	double time) const
{
	return std::exp(logP_of_parent_given_child(from_sequence, to_sequence, time));
}

template <typename T, bool log_space>
inline double substitution_model<T, log_space>::P_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence) const
{
	return P_of_parent_given_child(from_sequence, to_sequence, m_time);
}

// weighted test statistic
template <typename T, bool log_space>
double substitution_model<T, log_space>::weighted_distance(
	const population<T>& transmitter_population,
	const population<T>& recipient_population, double time) const
{
	const uint32_t K = transmitter_population.size();
	std::vector<double> weighted_counts(K, 0.0);
	double normalisation;

	const auto it_end_recip = recipient_population.cend();
	const auto it_end_trans = transmitter_population.cend();

	P_of_parent_given_child("A", "A", time);
	std::cout << std::right << std::scientific << std::setprecision(2);
	int col_width = transmitter_population.front().m_name.length() + 8;

	std::cout << std::setw(col_width) << "";
	for (const auto& i : transmitter_population)
	{
		std::cout << std::setw(col_width) << i.m_name;
	}
	std::cout << '\n';

	for (auto it = recipient_population.cbegin(); it != it_end_recip; ++it)
	{
		normalisation = 0;

		std::cout << std::setw(col_width) << it->m_name;
		for (const auto& i : transmitter_population)
		{
			normalisation += P_of_parent_given_child(i.m_name, it->m_name, time);
			std::cout << std::setw(col_width)
					  << P_of_parent_given_child(i.m_name, it->m_name, time);
		}
		std::cout << "\tSum: " << normalisation << '\n';

		if (normalisation == 0)
		{
			std::cerr << "ERROR: Could not match up transmitter and recipient "
						 "haplotypes!\n";
			exit(EXIT_FAILURE);
		}

		auto it_weighted = weighted_counts.begin();
		for (auto it_trans = transmitter_population.cbegin();
			 it_trans != it_end_trans; ++it_weighted, ++it_trans)
		{
			(*it_weighted) += P_of_parent_given_child(it_trans->m_name, it->m_name, time) / normalisation * it->m_count;
		}
	}

	count_type int_sum = 0;
	std::vector<count_type> rounded_weighted_counts(
		transmitter_population.size());
	std::transform(weighted_counts.cbegin(), weighted_counts.cend(),
		rounded_weighted_counts.begin(),
		[&int_sum](const double d) -> count_type
		{
			count_type result = std::round(d);
			int_sum += result;
			return result;
		});

	double dist = distance(
		K, transmitter_population.cbegin(), rounded_weighted_counts.cbegin(),
		transmitter_population.m_total_counts, int_sum, T::visitor);

#ifndef NDEBUG
	std::cout << "Vector of Transmitter Counts:\n";
	for (const auto& i : transmitter_population)
	{
		std::cout << '\t' << i.m_count;
	}

	std::cout << "\nVector of (weighted) Recipient Counts:\n";
	for (const auto& i : weighted_counts)
	{
		std::cout << '\t' << i;
	}

	std::cout << "\nVector of (weighted) Recipient Counts:\n";
	for (const auto& i : rounded_weighted_counts)
	{
		std::cout << '\t' << i;
	}
	std::cout << "\nDistance: " << dist << '\n';
#endif

	std::cout.unsetf(std::ios_base::scientific);

	// exit(0);
	return dist;
}

// partial specializations
template <bool log_space>
inline double substitution_model<trait, log_space>::P_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence,
	double time) const
{
	return (from_sequence == to_sequence);
}

template <bool log_space>
inline double substitution_model<trait, log_space>::P_of_parent_given_child(
	const std::string& from_sequence, const std::string& to_sequence) const
{
	return P_of_parent_given_child(from_sequence, to_sequence, 0.0);
}

template <bool log_space>
inline double substitution_model<trait, log_space>::weighted_distance(
	const population<trait>& transmitter_population,
	const population<trait>& recipient_population, double time) const
{
	return distance(
		transmitter_population.size(), transmitter_population.cbegin(),
		recipient_population.cbegin(), transmitter_population.m_total_counts,
		recipient_population.m_total_counts, trait::visitor, trait::visitor);
}

// debug display
template <typename T, bool log_space>
void substitution_model<T, log_space>::display(double time,
	bool display_raw) const
{
	get_entry_impl(T::index_to_elem(0), T::index_to_elem(0), time);

	std::cout << "Transition " << (log_space == true ? "(logarithmic) " : "")
			  << "tables for t = " << time << '\n' << std::fixed
			  << std::setprecision(3);
	double sum, average_sub_rate = 0;
	char from_elem, to_elem;

	for (std::size_t j = 0; j < N; ++j)
	{
		from_elem = T::index_to_elem(j);

		std::cout << "P(" << from_elem << " -> \n";
		sum = 0;

		for (std::size_t i = 0; i < N; ++i)
		{
			to_elem = T::index_to_elem(i);
			const double entry_ij = P_from_to_base(from_elem, to_elem, time);

			std::cout << T::index_to_elem(i) << ':' << entry_ij << ' ';
			sum += entry_ij;
		}

		average_sub_rate -= P_marginal_base(from_elem) * P_from_to_base(from_elem, from_elem, time);

		std::cout << '\n' << std::setw(4) << std::right << "sum"
				  << ": " << sum << "\n\n";
	}

	std::cout << "Average substitutions per time per site: " << average_sub_rate
			  << '\n';
}

/* private */
template <typename T, bool log_space>
double substitution_model<T, log_space>::get_entry_impl(char from_base,
	char to_base,
	double time) const
{
	if (time < 0)
	{
		throw std::invalid_argument(
			std::string("Time has to be non-negative! Exiting...\n"));
	}

	if ((time != m_time) || (m_cached == false))
#pragma omp critical
	{
#ifndef NDEBUG
		std::cout << "Resetting cache (time[" << time << "] != m_time[" << m_time
				  << "])\n";
#endif

		// 1. initialize
		std::array<double, N * N> p{};
		for (std::size_t i = 0; i < N; ++i)
		{
			p[i * (N + 1)] = 1;
		}

		// 2. create rate matrix
		matrix_type rates;
		double sum;
		for (std::size_t j = 0; j < N; ++j)
		{
			sum = 0;

			for (std::size_t i = 0; i < N; ++i)
			{
				rates[i][j] = m_rate * m_data_vector[i] * m_data_matrix[i][j];
				sum += rates[i][j];
			}

			rates[j][j] = -sum;
		}

		// 3. right-hand side of ODE
		auto rhs = [&rates](const std::array<double, N * N>& p,
			std::array<double, N * N>& dpdt,
			const double t) -> void
		{
			// std::cerr << t << '\n';

			for (std::size_t block = 0; block < N; ++block)
			{
				for (std::size_t i = 0; i < N; ++i)
				{
					dpdt[block * N + i] = 0;

					for (std::size_t j = 0; j < N; ++j)
					{
						dpdt[block * N + i] += rates[i][j] * p[block * N + j];
					}
				}
			}
		};

		// 4. perform integration
		/*
    boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled( 1.0E-6 , 1.0E-6,
    boost::numeric::odeint::runge_kutta_dopri5<std::array<double, N*N>>()),
            rhs,
            p,
            0.0,
            time,
            0.01);
    */

		boost::numeric::odeint::integrate(rhs, p, 0.0, time, 0.1);

		// 5. set transition matrix
		for (std::size_t i = 0; i < N; ++i)
		{
			for (std::size_t j = 0; j < N; ++j)
			{
				const double& ij = p[j * N + i];
				m_probability_matrix_time[i][j] = (log_space == true ? (ij ? std::log(ij) : log_lower_bound) : ij);
			}

			m_marginal_probability[i] = (log_space == true ? std::log(m_data_vector[i]) : m_data_vector[i]);
		}

		m_time = time;
		m_cached = true;
	}

	std::size_t from_base_i = T::elem_to_index(from_base);
	std::size_t to_base_i = T::elem_to_index(to_base);

	return m_probability_matrix_time[to_base_i][from_base_i];
}

template <typename T, bool log_space>
inline double
substitution_model<T, log_space>::logP_marginal_base(char elem) const
{
	if (log_space == true)
	{
		// log base
		const std::size_t elem_i = T::elem_to_index(elem);
		return m_marginal_probability[elem_i];
	}
	else
	{
		// non-log base
		const double p = P_marginal_base(elem);
		return (p ? std::log(p) : log_lower_bound);
	}
}

template <typename T, bool log_space>
inline double
substitution_model<T, log_space>::P_marginal_base(char elem) const
{
	if (log_space == true)
	{
		// log base
		return std::exp(logP_marginal_base(elem));
	}
	else
	{
		// non-log base
		const std::size_t elem_i = T::elem_to_index(elem);
		return m_marginal_probability[elem_i];
	}
}

#endif /* SUBSTITUTION_MODEL_HPP */