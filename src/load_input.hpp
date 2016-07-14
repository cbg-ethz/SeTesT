#ifndef FASTA_HPP
#define FASTA_HPP

#include <cctype>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename T>
void load_input(std::istream& input, std::vector<T>& transmitter_vec,
	std::vector<T>& recipient_vec,
	const std::regex& transmitter_regex,
	const std::regex& recipient_regex)
{
	transmitter_vec.clear();
	recipient_vec.clear();

	std::string id, seq, temp;
	uint32_t last_line = 1, line_num = 0;

	auto add_to_vec = [&](const std::string& id, std::string&& seq)
	{
		std::vector<std::string> split_vec;
		boost::split(split_vec, id, boost::is_any_of("_"),
			boost::token_compress_on);
		switch (split_vec.size())
		{
			case 1:
				std::cerr << "ERROR: Too few fields on line " << last_line << '\n';
				exit(EXIT_FAILURE);

			case 2:
				break;

			default:
				std::cerr << "ERROR: Too many fields on line " << last_line << '\n';
				exit(EXIT_FAILURE);
		}

		count_type count;
		try
		{
			count = boost::lexical_cast<count_type>(split_vec.back());
		}
		catch (boost::bad_lexical_cast&)
		{
			std::cerr << "ERROR: Count argument '" << split_vec.back()
					  << "' is not an integral/floating point value! Aborting.\n";
			exit(EXIT_FAILURE);
		}

		const int result = 2 * std::regex_search(split_vec.front(), transmitter_regex) + 1 * std::regex_search(split_vec.front(), recipient_regex);

		switch (result)
		{
			case 0:
				std::cerr << "WARNING: Entry on line " << last_line
						  << " is unrecognized as either transmitter or recipient\n";
				break;

			case 1:
				recipient_vec.emplace_back(
					std::move(seq), count,
					std::regex_replace(split_vec.front(), recipient_regex, ""));
				break;

			case 2:
				transmitter_vec.emplace_back(
					std::move(seq), count,
					std::regex_replace(split_vec.front(), transmitter_regex, ""));
				break;

			case 3:
				std::cerr
					<< "ERROR: Entry on line " << last_line
					<< " matches both transmitter and recipient regular expressions\n";
				exit(EXIT_FAILURE);
				break;
		}
	};

	while (input.good())
	{
		std::getline(input, temp);
		temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

		++line_num;
		if (temp.empty())
			continue;

		if (temp[0] == '>')
		{
			// identifier
			if (!(seq.empty()))
			{
				add_to_vec(id, std::move(seq));
			}
			last_line = line_num;
			id = temp.substr(1);
			seq.clear();
		}
		else
		{
			// sequence
			std::transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
			seq.append(temp);
		}
	}

	// add the last line
	if (!(seq.empty()))
	{
		add_to_vec(id, std::move(seq));
	}
}

/*
template <typename T>
std::vector<T> fasta_read(const std::string& input_file)
{
        std::ifstream input(input_file);

        if (input.is_open())
        {
                return fasta_read<T>(input);
        }
        else
        {
                std::cerr << "ERROR: Input file '" << input_file << "' could not
be opened!\n";
                exit(EXIT_FAILURE);
        }
}
*/

void load_input(std::istream& input, std::vector<trait>& transmitter_vec,
	std::vector<trait>& recipient_vec,
	const std::regex& transmitter_regex,
	const std::regex& recipient_regex)
{
	transmitter_vec.clear();
	recipient_vec.clear();

	std::vector<std::string> split_vec;
	std::string temp;
	uint32_t line_num = 0;

	while (input.good())
	{
		std::getline(input, temp);
		temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

		++line_num;
		if (temp.empty())
		{
			continue;
		}

		boost::split(split_vec, temp, boost::is_any_of("\t"),
			boost::token_compress_on);
		if (split_vec.size() == 3)
		{
			count_type count;
			try
			{
				count = boost::lexical_cast<count_type>(split_vec[2]);
			}
			catch (boost::bad_lexical_cast&)
			{
				std::cerr << "ERROR: Count argument '" << split_vec[2]
						  << "' is not an integral/floating point value! Aborting.\n";
				exit(EXIT_FAILURE);
			}

			const int result = 2 * std::regex_match(split_vec[0], transmitter_regex) + 1 * std::regex_match(split_vec[0], recipient_regex);

			switch (result)
			{
				case 0:
					std::cerr << "WARNING: Entry on line " << line_num
							  << " is unrecognized as either transmitter or recipient\n";
					break;

				case 1:
					recipient_vec.emplace_back(std::move(split_vec[1]), count);
					break;

				case 2:
					transmitter_vec.emplace_back(std::move(split_vec[1]), count);
					break;

				case 3:
					std::cerr
						<< "ERROR: Entry on line " << line_num
						<< " matches both transmitter and recipient regular expressions\n";
					exit(EXIT_FAILURE);
					break;
			}
		}
		else
		{
			std::cerr << "Line '" << line_num
					  << "' does not contain 3 tab-delimited fields.\n";
			exit(EXIT_FAILURE);
		}
	}
}

#endif /* FASTA_HPP */
