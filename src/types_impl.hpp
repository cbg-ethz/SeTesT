#ifndef TYPES_IMPL_HPP
#define TYPES_IMPL_HPP

#include <vector>
#include <iostream>
#include <limits>

#include <boost/algorithm/string.hpp>

///////////
// trait //
///////////

// ctor
trait::trait() = default;
trait::trait(const trait&) = default;
trait::trait(trait&&) = default;

trait::trait(const std::string& name_, count_type count_)
	: m_name(name_), m_count(count_) {}
trait::trait(std::string&& name_, count_type count_)
	: m_name(std::move(name_)), m_count(count_) {}

// assignment operators
trait& trait::operator=(const trait& other) = default;
trait& trait::operator=(trait&& other) = default;

template <typename T,
	typename std::enable_if<std::is_integral<typename std::remove_reference<T>::type>::value,
			  int>::type>
trait& trait::operator=(T&& count_)
{
	m_count = count_;
	return *this;
}

// visitors
inline count_type trait::visitor(const trait& trait_) { return trait_.m_count; }

inline bool trait::compare(const trait& A, const trait& B)
{
	return A.m_name < B.m_name;
}

// misc
inline bool trait::operator<(const trait& b) const
{
	return this->m_name < b.m_name;
}

std::ostream& operator<<(std::ostream& output, const trait& trait_) noexcept
{
	return output << "Name:\t" << trait_.m_name << '\n' << "Count:\t"
				  << trait_.m_count << '\n';
}

//////////////
// sequence //
//////////////

// ctor
sequence::sequence() = default;
sequence::sequence(const sequence&) = default;
sequence::sequence(sequence&&) = default;
sequence::sequence(std::string&& name_, count_type count_, std::string&& id_)
	: trait(std::move(name_), count_), m_id(std::move(id_)) {}

// assignment operators
sequence& sequence::operator=(const sequence& other) = default;
sequence& sequence::operator=(sequence&& other) = default;

// misc
// template <typename T, typename std::enable_if<std::is_same<std::string,
// T>::value, int>::type>
// sequence::sequence(T&& id, T&& seq, count_type count_) :
// trait(std::forward<T>(seq), count_), m_id(std::forward<T>(id)) {}

std::ostream& operator<<(std::ostream& output,
	const sequence& sequence_) noexcept
{
	return output << "ID:\t" << sequence_.m_id << '\n'
				  << static_cast<const trait&>(sequence_);
}

//////////////////////
// protein_sequence //
//////////////////////
constexpr std::array<char, 20> protein_sequence::valid_elem;

// ctor
protein_sequence::protein_sequence() = default;
protein_sequence::protein_sequence(const protein_sequence&) = default;
protein_sequence::protein_sequence(protein_sequence&&) = default;
protein_sequence::protein_sequence(std::string&& name_, count_type count_,
	std::string&& id_)
	: sequence(std::move(name_), count_, std::move(id_)) {}

// assignment operators
protein_sequence& protein_sequence::
operator=(const protein_sequence& other)
	= default;
protein_sequence& protein_sequence::
operator=(protein_sequence&& other)
	= default;

// misc
template <std::size_t N>
inline char index_to_elem_impl(std::size_t i,
	const std::array<char, N>& valid_chars)
{
#ifndef NDEBUG
	if (i > N)
	{
		throw std::invalid_argument(
			std::string("index 'i' is too large, has to be strictly smaller than " + std::to_string(N) + "! Exiting...\n"));
	}
#endif
	return valid_chars[i];
}

inline char protein_sequence::index_to_elem(std::size_t i)
{
	return index_to_elem_impl(i, valid_elem);
}

inline std::size_t protein_sequence::elem_to_index(char elem)
{
#ifndef NDEBUG
	if (elem < 65)
	{
		throw std::invalid_argument(
			std::string("index 'i' is too small, has to be strictly larger than "
						"64! Exiting...\n"));
	}
	if (elem > 89)
	{
		throw std::invalid_argument(
			std::string("index 'i' is too large, has to be strictly smaller than "
						"90! Exiting...\n"));
	}
#endif

	constexpr std::size_t max_range = std::numeric_limits<std::size_t>::max();
	constexpr std::array<std::size_t, 25> elem_map = { {
		0 /* A */, max_range, 1 /* C */, 2 /* D */, 3 /* E */, 4 /* F */,
		5 /* G */, 6 /* H */, 7 /* I */, max_range, 8 /* K */, 9 /* L */,
		10 /* M */, 11 /* N */, max_range, 12 /* P */, 13 /* Q */, 14 /* R */,
		15 /* S */, 16 /* T */, max_range, 17 /* V */, 18 /* W */, max_range,
		19 /* Y */
	} };

	const std::size_t result = elem_map[elem - 'A'];

#ifndef NDEBUG
	if (result == max_range)
	{
		throw std::invalid_argument(std::string("Element ") + elem + " is not a valid amino acid! Exiting...\n");
	}
#endif

	return result;
}

template <std::size_t N>
inline bool is_valid_elem_impl(char elem,
	const std::array<char, N>& valid_chars)
{
	return std::binary_search(valid_chars.cbegin(), valid_chars.cend(), elem);
}

inline bool protein_sequence::is_valid_elem(char elem)
{
	return is_valid_elem_impl(elem, valid_elem);
}

//////////////////
// dna_sequence //
//////////////////
constexpr std::array<char, 4> dna_sequence::valid_elem;

// ctor
dna_sequence::dna_sequence() = default;
dna_sequence::dna_sequence(const dna_sequence&) = default;
dna_sequence::dna_sequence(dna_sequence&&) = default;
dna_sequence::dna_sequence(std::string&& name_, count_type count_,
	std::string&& id_)
	: sequence(std::move(name_), count_, std::move(id_)) {}

// assignment operators
dna_sequence& dna_sequence::operator=(const dna_sequence& other) = default;
dna_sequence& dna_sequence::operator=(dna_sequence&& other) = default;

// misc
inline char dna_sequence::index_to_elem(std::size_t i)
{
	return index_to_elem_impl(i, valid_elem);
}

inline std::size_t dna_sequence::elem_to_index(char elem)
{
#ifndef NDEBUG
	if (elem < 65)
	{
		throw std::invalid_argument(
			std::string("index 'i' is too small, has to be strictly larger than "
						"64! Exiting...\n"));
	}
	if (elem > 84)
	{
		throw std::invalid_argument(
			std::string("index 'i' is too large, has to be strictly smaller than "
						"85! Exiting...\n"));
	}
#endif

	constexpr std::size_t max_range = std::numeric_limits<std::size_t>::max();
	constexpr std::array<std::size_t, 20> elem_map = { {
		0 /* A */, max_range, 1 /* C */, max_range, max_range, max_range,
		2 /* G */, max_range, max_range, max_range, max_range, max_range,
		max_range, max_range, max_range, max_range, max_range, max_range,
		max_range, 3 /* T */
	} };

	const std::size_t result = elem_map[elem - 'A'];

#ifndef NDEBUG
	if (result == max_range)
	{
		throw std::invalid_argument(std::string("Element ") + elem + " is not a valid amino acid! Exiting...\n");
	}
#endif

	return result;
}

inline bool dna_sequence::is_valid_elem(char elem)
{
	return is_valid_elem_impl(elem, valid_elem);
}

#endif /* TYPES_IMPL_HPP */