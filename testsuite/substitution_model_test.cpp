#include <config.h>

#undef NDEBUG

#include <iostream>
#include <iomanip>
#include <cmath>

#include "types.hpp"
#include "substitution_model.hpp"

int main(int argc, char* argv[])
{
	// 1) test amino acids
	substitution_model<protein_sequence, true> prot_log_space("./data/amino_acid", 4.6E-5);

	prot_log_space.display(0);
	prot_log_space.display(0.1);
	prot_log_space.display(1.0);
	prot_log_space.display(1E4);

	// 2) check divergences between logarithmic and non-logarithmic models
	substitution_model<protein_sequence, false> prot_normal_space("./data/amino_acid", 4.6E-5);
	double ij, first_val, second_val, sum = 0;
	std::cout << std::fixed << std::setprecision(3);
	for (const auto i : protein_sequence::valid_elem)
	{
		for (const auto j : protein_sequence::valid_elem)
		{
			first_val = prot_log_space.P_from_to_base(i, j, 0.5);
			second_val = prot_normal_space.P_from_to_base(i, j, 0.5);

			ij = 2 * std::fabs(first_val - second_val) / (std::fabs(first_val) + std::fabs(second_val));
			sum += ij;
			std::cout << std::setw(7) << std::right << ij;
		}

		std::cout << '\n';
	}

	std::cout << "Sum of differences: " << std::setprecision(15) << sum << '\n';

	return (sum < std::numeric_limits<float>::epsilon() ? EXIT_SUCCESS : EXIT_FAILURE);
}