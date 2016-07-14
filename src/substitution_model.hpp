#ifndef SUBSTITUTION_MODEL_HPP
#define SUBSTITUTION_MODEL_HPP

#include "defs.hpp"

#include <array>
#include <string>
#include <functional>
#include <utility>

#include "types.hpp"
#include "population.hpp"

template <typename T, bool log_space = true>
class substitution_model
{
public:
	static constexpr std::size_t N = T::N;
	static constexpr double log_lower_bound = -30;

	typedef double matrix_type[N][N];
	typedef double vector_type[N];

	// ctor + initializer
	void init(const std::string& file_stem);
	substitution_model(std::string file_stem, double rate);

	// set average base substitution rate
	void set_avg_substitution_rate(double rate);

	// transition probabilities
	inline double logP_from_to_seq(const std::string& from_sequence,
		const std::string& to_sequence,
		double time) const;
	inline double logP_from_to_seq(const std::string& from_sequence,
		const std::string& to_sequence) const;
	inline double P_from_to_seq(const std::string& from_sequence,
		const std::string& to_sequence,
		double time) const;
	inline double P_from_to_seq(const std::string& from_sequence,
		const std::string& to_sequence) const;

	inline double logP_from_to_base(char from_elem, char to_elem,
		double time) const;
	inline double logP_from_to_base(char from_elem, char to_elem) const;
	inline double P_from_to_base(char from_elem, char to_elem, double time) const;
	inline double P_from_to_base(char from_elem, char to_elem) const;

	inline double logP_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence,
		double time) const;
	inline double
	logP_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence) const;
	inline double P_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence,
		double time) const;
	inline double P_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence) const;

	// weighted test statistic
	double weighted_distance(const population<T>& transmitter_population,
		const population<T>& recipient_population,
		double time) const;

	// debug display
	void display(double time, bool display_raw = false) const;

private:
	double get_entry_impl(char from_base, char to_base, double time) const;

	inline double logP_marginal_base(char elem) const;
	inline double P_marginal_base(char elem) const;

	/* data values */
	double m_rate;
	matrix_type m_data_matrix;
	vector_type m_data_vector;

	/* cached values to improve lookup times */
	mutable matrix_type m_probability_matrix_time;
	mutable vector_type m_marginal_probability;

	mutable double m_time;
	mutable bool m_cached;
};

template <bool log_space>
class substitution_model<trait, log_space>
{
public:
	// empty class, as "traits" are too generic to have a substitution model
	inline double P_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence,
		double time) const;
	inline double P_of_parent_given_child(const std::string& from_sequence,
		const std::string& to_sequence) const;

	double weighted_distance(const population<trait>& transmitter_population,
		const population<trait>& recipient_population,
		double time) const;
};

#include "substitution_model_impl.hpp"

#endif /* SUBSTITUTION_MODEL_HPP */