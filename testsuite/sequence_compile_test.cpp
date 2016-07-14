#undef NDEBUG

#include <types.hpp>

int main(int argc, char* argv[])
{
	// 1) test amino acids
	protein_sequence prot;
	std::cout << "Protein Test\n";

	for (const char c1 : { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' })
	{
		std::size_t i1 = prot.elem_to_index(c1);
		const char c2 = prot.index_to_elem(i1);
		std::size_t i2 = prot.elem_to_index(c2);

		std::cout
			<< "c1 : " << c1 << '\n'
			<< "i1 : " << i1 << '\n'
			<< "c2 : " << c2 << '\n'
			<< "i2 : " << i2 << '\n'
			<< '\n';

		if ((c1 != c2) || (i1 != i2))
		{
			return EXIT_FAILURE;
		}
	}

	// 2) test DNA bases
	dna_sequence dna;
	std::cout << "DNA Test\n";

	for (const char c1 : { 'A', 'C', 'G', 'T' })
	{
		std::size_t i1 = dna.elem_to_index(c1);
		const char c2 = dna.index_to_elem(i1);
		std::size_t i2 = dna.elem_to_index(c2);

		std::cout
			<< "c1 : " << c1 << '\n'
			<< "i1 : " << i1 << '\n'
			<< "c2 : " << c2 << '\n'
			<< "i2 : " << i2 << '\n'
			<< '\n';

		if ((c1 != c2) || (i1 != i2))
		{
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}