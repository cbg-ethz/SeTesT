#ifndef SAMPLING_IMPL_HPP
#define SAMPLING_IMPL_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <utility>

#include <boost/math/tools/roots.hpp>

/* helper functions */
count_type generate_zero_trun_poisson(double lambda,
	std::default_random_engine& rng)
{
	// https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html
	const double u = std::uniform_real_distribution<>()(rng);
	const double t_first_event = -std::log(1 - u * (1 - std::exp(-lambda)));

	const double remain = lambda - t_first_event;
	return std::poisson_distribution<>(remain)(rng) + 1;
}

double mean_to_lambda(double mu)
{
	boost::uintmax_t max_iter = 1000;

	auto result = boost::math::tools::toms748_solve([mu](double lambda) {
		return (mu - lambda / (1 - std::exp(-lambda)));
	},
		1E-5, 1E5, boost::math::tools::eps_tolerance<double>(30), max_iter);

	return (result.first + result.second) / 2;
}

double lambda_to_mean(double lambda)
{
	return lambda / (1 - std::exp(-lambda));
}

template <typename T>
constexpr const T& identity(const T& value)
{
	return value;
}

/* L2-norm between two probabilities */
template <typename InputIt1, typename InputIt2, typename V1, typename V2>
inline double distance(uint32_t K, InputIt1 T_iter, InputIt2 R_iter,
	double T_total, double R_total, V1 T_visitor,
	V2 R_visitor)
{
	double result = 0, diff;

	const InputIt1 end_it = T_iter + K;
	for (; T_iter != end_it; ++T_iter, ++R_iter)
	{
		diff = T_visitor(*T_iter) / T_total - R_visitor(*R_iter) / R_total;
		result += diff * diff;
	}

	return result;
}

template <typename InputIt1, typename InputIt2, typename V>
inline double distance(uint32_t K, InputIt1 T_iter, InputIt2 R_iter,
	double T_total, double R_total, V visitor)
{
	using input_type2 = typename std::iterator_traits<InputIt2>::value_type;
	return distance(K, T_iter, R_iter, T_total, R_total, visitor,
		identity<input_type2>);
}

template <typename InputIt1, typename InputIt2>
inline double distance(uint32_t K, InputIt1 T_iter, InputIt2 R_iter,
	double T_total, double R_total)
{
	using input_type1 = typename std::iterator_traits<InputIt1>::value_type;
	using input_type2 = typename std::iterator_traits<InputIt2>::value_type;
	return distance(K, T_iter, R_iter, T_total, R_total, identity<input_type1>,
		identity<input_type2>);
}

/* Dirichlet */
void gen_rand_dirichlet(const double* input_begin_iter,
	const double* input_end_iter, double* output_begin_iter,
	std::default_random_engine& rng)
{
	double sum = 0;

	// 1. first sample from gamma distribution
	double* output_end_iter = std::transform(input_begin_iter, input_end_iter, output_begin_iter,
		[&sum, &rng](const double& val) -> double {
			assert(val > 0);

			double x = std::gamma_distribution<>(val, 1.0)(rng);
			sum += x;
			return x;
		});

	// 2. then normalize the vector
	std::for_each(output_begin_iter, output_end_iter,
		[&sum](double& val) {
			val /= sum;
		});
}

/* Multinomial */
template <typename InputIt, typename OutputIt, typename V>
void gen_rand_multinomial(InputIt input_begin_iter, InputIt input_end_iter,
	OutputIt output_begin_iter, count_type N,
	std::default_random_engine& rng, double sum,
	V visitor)
{
	typedef typename std::iterator_traits<InputIt>::value_type input_type;

	if (!sum)
	{
		sum = std::accumulate(input_begin_iter, input_end_iter, 0.0,
			[visitor](double d, const input_type& val) {
				return d + visitor(val);
			});
	}

	// std::cerr << "Initial Sum: " << sum << '\n';
	// std::cerr << "Initial   N: " << N << '\n';

	OutputIt output_end_iter = std::transform(input_begin_iter, input_end_iter - 1, output_begin_iter,
		[&](const input_type& val) -> count_type {
			// std::cerr << "sum: " << sum << '\n';

			assert(sum >= 0);

			if (N)
			{
				const double value = visitor(val);
				double p = value / sum;

				if (p < 0)
				{
					p = 0;
				}
				if (p > 1)
				{
					p = 1;
				}

				count_type X_i = std::binomial_distribution<>(N, p)(rng);

				// std::cerr << "val: " << value << '\n';
				// std::cerr << "  N: " << N << '\n';
				// std::cerr << "  p: " << p << '\n';
				// std::cerr << "  X: " << X_i << '\n' << '\n';

				sum -= value;
				N -= X_i;

				return X_i;
			}
			else
			{
				return 0;
			}
		});

	// std::cerr << '\n' << '\n';

	*output_end_iter = N;
}

template <typename InputIt, typename OutputIt>
void gen_rand_multinomial(InputIt input_begin_iter, InputIt input_end_iter,
	OutputIt output_begin_iter, count_type N,
	std::default_random_engine& rng)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	gen_rand_multinomial(input_begin_iter, input_end_iter, output_begin_iter, N,
		rng, 0.0, identity<input_type>);
}

template <typename InputIt, typename OutputIt>
void gen_rand_multinomial(InputIt input_begin_iter, InputIt input_end_iter,
	OutputIt output_begin_iter, count_type N,
	std::default_random_engine& rng, double sum)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	gen_rand_multinomial(input_begin_iter, input_end_iter, output_begin_iter, N,
		rng, sum, identity<input_type>);
}

/* Generate transmitter sample */
template <typename InputIt, typename OutputIt, typename V>
void gen_rand_transmitter_sample(uint32_t K, InputIt p_T_iter, OutputIt T_iter,
	count_type transmitter_coverage,
	std::default_random_engine& rng, V visitor)
{
	gen_rand_multinomial(p_T_iter, p_T_iter + K, T_iter, transmitter_coverage,
		rng, 1.0, visitor);
}

template <typename InputIt, typename OutputIt>
void gen_rand_transmitter_sample(uint32_t K, InputIt p_T_iter, OutputIt T_iter,
	count_type transmitter_coverage,
	std::default_random_engine& rng)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	gen_rand_transmitter_sample(K, p_T_iter, T_iter, transmitter_coverage, rng,
		identity<input_type>);
}

/* Generate recipient sample */
template <typename InputIt, typename OutputIt, typename V>
count_type
gen_rand_recipient_sample(uint32_t K, InputIt p_R_iter,
	OutputIt Bottleneck_iter, OutputIt R_iter,
	count_type recipient_coverage, double lambda,
	std::default_random_engine& rng, V visitor)
{
	// 1st step: generate bottlenecked sample
	count_type bottleneck_size = generate_zero_trun_poisson(lambda, rng);
	gen_rand_multinomial(p_R_iter, p_R_iter + K, Bottleneck_iter, bottleneck_size,
		rng, 1.0);

	// 2nd step: generate actual sample after bottleneck
	gen_rand_multinomial(Bottleneck_iter, Bottleneck_iter + K, R_iter,
		recipient_coverage, rng, bottleneck_size);

	return bottleneck_size;
}

template <typename InputIt, typename OutputIt>
count_type gen_rand_recipient_sample(uint32_t K, InputIt p_R_iter,
	OutputIt Bottleneck_iter, OutputIt R_iter,
	count_type recipient_coverage,
	double lambda,
	std::default_random_engine& rng)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	return gen_rand_recipient_sample(K, p_R_iter, Bottleneck_iter, R_iter,
		recipient_coverage, lambda, rng,
		identity<input_type>);
}

/* Generate both transmitter and recipient sample */
template <typename InputIt, typename OutputIt, typename V>
count_type
generate_transmission_sample(uint32_t K, InputIt p_T_iter, InputIt p_R_iter,
	OutputIt T_iter, OutputIt Bottleneck_iter,
	OutputIt R_iter, count_type transmitter_coverage,
	count_type recipient_coverage, double lambda,
	std::default_random_engine& rng, V visitor)
{
	// transmitter sample
	gen_rand_transmitter_sample(K, p_T_iter, T_iter, transmitter_coverage, rng,
		visitor);

	// recipient sample
	return gen_rand_recipient_sample(K, p_R_iter, Bottleneck_iter, R_iter,
		recipient_coverage, lambda, rng, visitor);
}

template <typename InputIt, typename OutputIt>
count_type
generate_transmission_sample(uint32_t K, InputIt p_T_iter, InputIt p_R_iter,
	OutputIt T_iter, OutputIt Bottleneck_iter,
	OutputIt R_iter, count_type transmitter_coverage,
	count_type recipient_coverage, double lambda,
	std::default_random_engine& rng)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	return generate_transmission_sample(K, p_T_iter, p_R_iter, T_iter,
		Bottleneck_iter, R_iter,
		transmitter_coverage, recipient_coverage,
		lambda, rng, identity<input_type>);
}

/* actual test fucntions */
template <typename InputIt, typename V>
double perform_test(uint32_t K, InputIt T_iter, count_type T_total,
	count_type R_total, double test_statistic, double lambda,
	count_type num_trials, std::default_random_engine& rng,
	V visitor)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	count_type num_more_extremes = 0;

	double p_T_MLE[K];
	std::transform(T_iter, T_iter + K, p_T_MLE,
		[T_total, visitor](const input_type& val) -> double {
			return visitor(val) / static_cast<double>(T_total);
		});

	count_type temp_bottleneck[K];
	count_type temp_recip[K];
	double rand_dist;

	for (count_type i = 0; i < num_trials; ++i)
	{
		gen_rand_recipient_sample(K, p_T_MLE, temp_bottleneck, temp_recip, R_total,
			lambda, rng);
		rand_dist = distance(K, T_iter, temp_recip, T_total, R_total, visitor);

		num_more_extremes += (rand_dist >= test_statistic);
	}

	return static_cast<double>(num_more_extremes) / num_trials;
}

template <typename InputIt>
double perform_test(uint32_t K, InputIt T_iter, count_type T_total,
	count_type R_total, double test_statistic, double lambda,
	count_type num_trials, std::default_random_engine& rng)
{
	using input_type = typename std::iterator_traits<InputIt>::value_type;
	return perform_test(K, T_iter, T_total, R_total, test_statistic, lambda,
		num_trials, rng, identity<input_type>);
}

#endif /* SAMPLING_IMPL_HPP */