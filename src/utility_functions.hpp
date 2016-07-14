#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

/*
 * Copyright (c) 2016 David Seifert
 *
 * This file is part of ngshmmalign
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */

#include <string>
#include <array>
#include <cmath>

template <typename T>
inline std::string number_to_sequence(std::size_t number, int L)
{
	std::string result(L, 'A');
	int j;

	for (int i = L - 1; i >= 0; --i)
	{
		j = std::pow(T::N, i);
		result[L - 1 - i] = T::index_to_elem(number / j);
		number %= j;
	}

	return result;
}

template <typename T>
inline std::size_t sequence_to_number(const std::string& sequence)
{
	const int L = sequence.length();
	std::size_t result = 0;
	for (int i = 0; i < L; ++i)
	{
		result += std::pow(T::N, i) * T::elem_to_index(L - 1 - i);
	}

	return result;
}

#endif /* UTILITY_FUNCTIONS_HPP */