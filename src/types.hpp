#ifndef TYPES_HPP
#define TYPES_HPP

#include "defs.hpp"

#include <array>
#include <memory>
#include <string>
#include <type_traits>

struct trait
{
	std::string m_name;
	count_type m_count;

	// ctor
	trait();
	trait(const trait&);
	trait(trait&&);

	trait(const std::string& name_, count_type count_);
	trait(std::string&& name_, count_type count_);

	// assignment operators
	trait& operator=(const trait& other);
	trait& operator=(trait&& other);

	template <
		typename T,
		typename std::enable_if<std::is_integral<typename std::remove_reference<T>::type>::value,
			int>::type
		= 0>
	trait& operator=(T&& count_);

	// visitors
	static inline count_type visitor(const trait& trait_);
	static inline bool compare(const trait& A, const trait& B);

	// misc
	inline bool operator<(const trait& b) const;
	friend std::ostream& operator<<(std::ostream& output,
		const trait& trait_) noexcept;
};

struct sequence : public trait
{
	std::string m_id;

	// ctor
	using trait::trait;
	sequence();
	sequence(const sequence&);
	sequence(sequence&&);
	sequence(std::string&& name_, count_type count_, std::string&& id_);

	// assignment operators
	using trait::operator=;
	sequence& operator=(const sequence& other);
	sequence& operator=(sequence&& other);

	// visitors
	using trait::visitor;

	// misc
	friend std::ostream& operator<<(std::ostream& output,
		const sequence& sequence_) noexcept;
};

struct protein_sequence : public sequence
{
	static constexpr std::size_t N = 20;
	static constexpr std::array<char, N> valid_elem = {
		{ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
			'R', 'S', 'T', 'V', 'W', 'Y' }
	};

	// ctor
	using sequence::sequence;
	protein_sequence();
	protein_sequence(const protein_sequence&);
	protein_sequence(protein_sequence&&);
	protein_sequence(std::string&& name_, count_type count_, std::string&& id_);

	// assignment operators
	using sequence::operator=;
	protein_sequence& operator=(const protein_sequence& other);
	protein_sequence& operator=(protein_sequence&& other);

	// visitors
	using sequence::visitor;

	// misc
	static inline char index_to_elem(std::size_t i);
	static inline std::size_t elem_to_index(char elem);

	static inline bool is_valid_elem(char elem);
};

struct dna_sequence : public sequence
{
	static constexpr std::size_t N = 4;
	static constexpr std::array<char, N> valid_elem = { { 'A', 'C', 'G', 'T' } };

	// ctor
	using sequence::sequence;
	dna_sequence();
	dna_sequence(const dna_sequence&);
	dna_sequence(dna_sequence&&);
	dna_sequence(std::string&& name_, count_type count_, std::string&& id_);

	// assignment operators
	using sequence::operator=;
	dna_sequence& operator=(const dna_sequence& other);
	dna_sequence& operator=(dna_sequence&& other);

	// visitors
	using sequence::visitor;

	// misc
	static inline char index_to_elem(std::size_t i);
	static inline std::size_t elem_to_index(char elem);

	static inline bool is_valid_elem(char elem);
};

#include "types_impl.hpp"

#endif /* TYPES_HPP */