#include <config.h>

#undef NDEBUG

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

bool verbose = true;

#include "substitution_model.hpp"
#include "types.hpp"
#include "utility_functions.hpp"

int main(int argc, char* argv[])
{
	const int L = 2;

	double sum_p = 0, p;
	double sum_q = 0, q;

	const double t = (argc > 1 ? std::atoi(argv[1]) : 100);

	std::string protein;
	substitution_model<protein_sequence, false> prot_model("./data/amino_acid", 4.6E-5);

	const std::string orig_protein(L, 'Y');

	std::cout << std::fixed << std::setprecision(8);

	for (int i = 0; i < 400; ++i)
	{
		protein = number_to_sequence<protein_sequence>(i, L);
		p = prot_model.P_of_parent_given_child(protein, orig_protein, t);
		sum_p += p;

		q = prot_model.P_of_parent_given_child(orig_protein, protein, t);
		sum_q += q;

		std::cout << protein << ": " << p << " (" << q << ")\n";
	}

	std::cout << "Sum: " << sum_p << " (incorrect sum when switching parent and child: " << sum_q << ")\n";

	return (((sum_p > 0.9) && (sum_p < 1.1)) ? EXIT_SUCCESS : EXIT_FAILURE);
}