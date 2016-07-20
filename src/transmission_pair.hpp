#ifndef TRANSMISSION_PAIR_HPP
#define TRANSMISSION_PAIR_HPP

#include "defs.hpp"

#include <memory>
#include <regex>
#include <string>
#include <type_traits>
#include <vector>

#include "population.hpp"
#include "substitution_model.hpp"
#include "types.hpp"

class transmission_pair
{
public:
	virtual void init(std::istream& input_, const std::regex transmitter_regex_,
		const std::regex recipient_regex_,
		bool retain_unshared_)
		= 0;

	virtual double pvalue(count_type num_trials, double time, double lambda,
		uint64_t seed) const = 0;

	virtual void simulate_fixed(uint32_t K, count_type num_simulations,
		count_type num_trials,
		count_type transmitter_coverage,
		count_type recipient_coverage, double lambda,
		double p0, double f0, uint64_t seed) const = 0;

	virtual void simulate_variable(uint32_t K, count_type num_simulations,
		count_type num_trials,
		count_type transmitter_coverage,
		count_type recipient_coverage, double lambda,
		double p0, double f0, double alpha,
		uint64_t seed) const = 0;

	virtual ~transmission_pair() = default;
};

template <typename T>
class transmission_pair_impl : public transmission_pair,
							   private substitution_model<T>
{
public:
	transmission_pair_impl();
	transmission_pair_impl(const transmission_pair_impl&);
	transmission_pair_impl(transmission_pair_impl&&);
	transmission_pair_impl& operator=(const transmission_pair_impl& other);
	transmission_pair_impl& operator=(transmission_pair_impl&& other);

	transmission_pair_impl(const std::string& file_stem_parameters_, double rate);

	virtual void init(std::istream& input_, const std::regex transmitter_regex_,
		const std::regex recipient_regex_,
		bool retain_unshared_) override;

	virtual double pvalue(count_type num_trials, double time, double lambda,
		uint64_t seed) const override;

	virtual void simulate_fixed(uint32_t K, count_type num_simulations,
		count_type num_trials,
		count_type transmitter_coverage,
		count_type recipient_coverage, double lambda,
		double p0, double f0,
		uint64_t seed) const override;

	virtual void simulate_variable(uint32_t K, count_type num_simulations,
		count_type num_trials,
		count_type transmitter_coverage,
		count_type recipient_coverage, double lambda,
		double p0, double f0, double alpha,
		uint64_t seed) const override;

private:
	void init_(std::vector<T>& transmitter_population_,
		std::vector<T>& recipient_population_, bool retain_unshared_);

	population<T> m_transmitter_population;
	population<T> m_recipient_population;
};

#include "transmission_pair_impl.hpp"

#endif /* TRANSMISSION_PAIR_HPP */