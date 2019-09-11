/*
 * main.cpp
 *
 *  Created on: Sep 11, 2019
 *      Author: Sarah Lutteropp
 */

#include <stddef.h>
#include <cstdlib>
#include <random>
#include <vector>

#include "SimpleMVPTree.hpp"

double distance(const unsigned int& a, const unsigned int& b) {
	return (double) abs(a - b);
}

int main() {
	std::vector<unsigned int> numbers;

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> dist(1, 999999999); // distribution in range [1, 999999999]

	for (size_t i = 0; i < 100000; ++i) {
		numbers.push_back(dist(rng));
	}

	SimpleMVPTree<unsigned int> smvp(numbers, distance);

	unsigned int query = 10;
	double min_distance = smvp.search_mindist(query);
	unsigned int mindist_element = smvp.search_mindist_elem(query);

	std::cout << "Everything worked fine\n";

	return 0;
}
