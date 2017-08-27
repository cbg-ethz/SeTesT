#ifndef POPULATION_HPP
#define POPULATION_HPP

#include "defs.hpp"

#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "types.hpp"

template <typename T>
class population : private std::vector<T>
{
public:
	// ctor
	population();
	population(const population& population_);
	population(population&& population_);
	population& operator=(const population& population_);
	population& operator=(population&& population_);

	population(const std::vector<T>& vec_);
	population(std::vector<T>&& vec_);

	population(uint32_t K, count_type total_counts_);

	// exposed member functions of std::vector<T>
	using typename std::vector<T>::value_type;
	using std::vector<T>::vector;
	using std::vector<T>::resize;
	using std::vector<T>::begin;
	using std::vector<T>::end;
	using std::vector<T>::cbegin;
	using std::vector<T>::cend;
	using std::vector<T>::swap;
	using std::vector<T>::clear;
	using std::vector<T>::reserve;
	using std::vector<T>::size;
	using std::vector<T>::push_back;
	using std::vector<T>::front;

	// print
	friend std::ostream& operator<<(std::ostream& output,
		const population& population_) noexcept
	{
		std::for_each(population_.cbegin(), population_.cend(),
			[&output](const T& val) {
				output << val << '\n';
			});
		return output;
	}

	double distance(const population& recipient) const;

	count_type m_total_counts;
};

#include "population_impl.hpp"

#endif /* POPULATION_HPP */