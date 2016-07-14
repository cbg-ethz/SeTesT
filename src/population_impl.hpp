#ifndef POPULATION_IMPL_HPP
#define POPULATION_IMPL_HPP

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <random>

#include <boost/math/tools/roots.hpp>

#include "population.hpp"
#include "sampling_impl.hpp"

// Constructors
template <typename T>
population<T>::population() = default;

template <typename T>
population<T>::population(const population& population_) = default;

template <typename T>
population<T>::population(population&& population_) = default;

template <typename T>
population<T>& population<T>::
operator=(const population<T>& population_)
	= default;

template <typename T>
population<T>& population<T>::operator=(population<T>&& population_) = default;

template <typename T>
population<T>::population(const std::vector<T>& vec_)
	: std::vector<T>(vec_),
	  m_total_counts(std::accumulate(
		  this->cbegin(), this->cend(), 0ul,
		  [](count_type A, const T& B)
		  {
			  return A + B.m_count;
		  })){};

template <typename T>
population<T>::population(std::vector<T>&& vec_)
	: std::vector<T>(std::move(vec_)),
	  m_total_counts(std::accumulate(
		  this->cbegin(), this->cend(), 0ul,
		  [](count_type A, const T& B)
		  {
			  return A + B.m_count;
		  })){};

template <typename T>
population<T>::population(uint32_t K, count_type total_counts_)
	: std::vector<T>(K), m_total_counts(total_counts_){};

// member functions
template <typename T>
double population<T>::distance(const population& recipient) const
{
	assert(this->size() == recipient.size());

	double result = 0, diff;
	const double trans_total = this->m_total_counts;
	const double recip_total = recipient.m_total_counts;
	const auto end_it = this->cend();
	for (auto i = this->cbegin(), j = recipient.cbegin(); i != end_it; ++i, ++j)
	{
		diff = (i->m_count / trans_total) - (j->m_count / recip_total);
		result += diff * diff;
	}
	return result;
}

#endif /* POPULATION_IMPL_HPP */