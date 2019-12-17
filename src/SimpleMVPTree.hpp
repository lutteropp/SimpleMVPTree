/*
 * SimpleMVPTree.hpp
 *
 * Created on: Sep 11, 2019
 * Author: Sarah Lutteropp
 */

#pragma once

#include <bits/move.h>
#include <stddef.h>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template<typename T>
class SimpleMVPTree {
public:
	SimpleMVPTree(const std::vector<T>& data, std::function<double(const T&, const T&)> dist) :
			data { data }, tau { std::numeric_limits<double>::max() }, distance { dist }, root { nullptr } {
		num_vantage_points = ceil(log(data.size() / 2));
		vp_indices.resize(num_vantage_points);
		distToVP.resize(data.size());
		for (size_t i = 0; i < distToVP.size(); ++i) {
			distToVP[i].resize(num_vantage_points);
		}
		items.resize(data.size());
		for (size_t i = 0; i < items.size(); ++i) {
			items[i] = i;
		}

		root = buildFromPoints(0, items.size());
	}

	double search_mindist(const T& target) {
		tau = std::numeric_limits<double>::max();

		std::vector<double> vp_dist(num_vantage_points);
		for (size_t i = 0; i < num_vantage_points; ++i) {
			vp_dist[i] = distance(data[vp_indices[i]], target);
			if (vp_dist[i] < tau) {
				tau = vp_dist[i];
				closest_idx = vp_indices[i];
			}
			if (tau == 0) {
				return 0;
			}
		}
		search(root.get(), target, vp_dist, 0);

		// Uncomment to check if the result is correct (this is very slow!)
		/*size_t min_idx = 0;
		 double min = std::numeric_limits<double>::infinity();
		 for (size_t i = 0; i < _items.size(); ++i) {
		 double dist = distance(target, _data[i]);
		 if (dist < min) {
		 min = dist;
		 min_idx = i;
		 }
		 }
		 if (_tau != min) {
		 std::cout << "ERROR!!! THE RESULT IS WRONG!!!\n";
		 std::cout << "_tau: " << _tau << "\n";
		 std::cout << "min: " << min << "\n";
		 std::cout << "witness: " << _data[min_idx] << "\n";
		 throw std::runtime_error("Stopping now");
		 }*/

		return tau;
	}

	T search_mindist_elem(const T& target) {
		search_mindist(target);
		return data[closest_idx];
	}

private:
	struct DirectoryNode {
		std::unique_ptr<DirectoryNode> left;
		std::unique_ptr<DirectoryNode> right;
		unsigned int s1, e1;
		unsigned int s2, e2;
		double threshold; // median distance
	};

	unsigned int num_vantage_points;
	size_t closest_idx;
	double tau;
	std::function<double(const T&, const T&)> distance;
	std::unique_ptr<DirectoryNode> root;
	const std::vector<T>& data;
	std::vector<unsigned int> items;
	std::vector<unsigned int> vp_indices;
	std::vector<std::vector<double> > distToVP;

	unsigned int findVantagePoint(unsigned int lower, unsigned int upper) {
		// chose a random node
		unsigned int cand = (unsigned int) ((double) rand() / RAND_MAX * (upper - lower - 1)) + lower;
		// find the node that is farthest away from this node as VP point
		double maxDist = 0;
		unsigned int vpPoint = lower;

		static const unsigned int NUM_SAMPLES = 200;
		static const unsigned int MAX_COLLISIONS = 5;
		std::unordered_set<unsigned int> samples;
		unsigned int actCollisions = 0;
		while (samples.size() < (upper - lower - 1) && samples.size() < NUM_SAMPLES) {
			// try to randomly pick a sample
			unsigned int sam = (unsigned int) ((double) rand() / RAND_MAX * (upper - lower - 1)) + lower;
			if (sam != cand) {
				if (samples.find(sam) != samples.end()) {
					actCollisions++;
					if (actCollisions >= MAX_COLLISIONS) {
						break;
					}
				}
				samples.insert(sam);
			}
		}

		for (unsigned int i : samples) {
			double dist = distance(data[items[cand]], data[items[i]]);
			if (dist > maxDist) {
				maxDist = dist;
				vpPoint = i;
			}
		}
		return vpPoint;
	}

	std::unique_ptr<DirectoryNode> makeDirectory(unsigned int lower, unsigned int upper, unsigned int actVPIndex) {
		std::unique_ptr<DirectoryNode> node = std::make_unique<DirectoryNode>();
		if (lower >= upper) {
			return node;
		}

		unsigned int median = (upper + lower) / 2;
		std::nth_element(items.begin() + lower, items.begin() + median, items.begin() + upper, [&](unsigned int i, unsigned int j) {
			return distToVP[i][actVPIndex] < distToVP[j][actVPIndex];
		});
		node->s1 = lower;
		node->e1 = median;
		node->s2 = median;
		node->e2 = upper;
		node->threshold = distToVP[items[median]][actVPIndex];

		if (actVPIndex < num_vantage_points - 1) {
			node->left = makeDirectory(node->s1, node->e1, actVPIndex + 1);
			node->right = makeDirectory(node->s2, node->e2, actVPIndex + 1);
		}

		return node;
	}

	std::unique_ptr<DirectoryNode> buildFromPoints(unsigned int lower, unsigned int upper) {
		std::unique_ptr<DirectoryNode> root = std::make_unique<DirectoryNode>();
		if (upper == lower) {
			return root;
		}
		unsigned int actVPIndex = 0;

		unsigned int maxDist = 0;
		unsigned int maxDistIdx = 0;

		// first Vantage point is a bit different from the rest:
		unsigned int vp = findVantagePoint(0, upper);

		vp_indices[0] = items[vp];
		std::swap(items[lower], items[vp]);
		for (size_t j = lower + 1; j < upper; ++j) {
			double dist = distance(data[items[j]], data[items[lower]]);
			distToVP[items[j]][0] = dist;
			if (dist > maxDist) {
				maxDist = dist;
				maxDistIdx = items[j];
			}
		}
		// select remaining Vantage points
		for (size_t i = 1; i < num_vantage_points; ++i) {
			unsigned int vp = maxDistIdx;
			vp_indices[i] = items[vp];
			std::swap(items[lower + i], items[vp]);
			maxDist = 0;
			maxDistIdx = 0;
			for (size_t j = lower + i + 1; j < upper; ++j) {
				double dist = distance(data[items[j]], data[items[lower + i]]);
				distToVP[items[j]][i] = dist;
				// look at the minimum distance from vantage points so far
				for (size_t k = 1; k < i; ++k) {
					dist = std::min(dist, distToVP[items[j]][k]);
				}
				if (dist > maxDist) {
					maxDist = dist;
					maxDistIdx = items[j];
				}
			}
		}
		// reorder the items and store it in the directory.
		unsigned int median = (upper + lower + num_vantage_points) / 2;
		// partition around the median distance from first VP
		std::nth_element(items.begin() + lower + num_vantage_points, items.begin() + median, items.begin() + upper,
				[&](unsigned int i, unsigned int j) {
					return distToVP[i][actVPIndex] < distToVP[j][actVPIndex];
				});
		root->s1 = lower + num_vantage_points;
		root->e1 = median;
		root->s2 = median;
		root->e2 = upper;
		root->threshold = distToVP[items[median]][actVPIndex];

		if (actVPIndex < num_vantage_points - 1) {
			root->left = makeDirectory(root->s1, root->e1, actVPIndex + 1);
			root->right = makeDirectory(root->s2, root->e2, actVPIndex + 1);
		}
		return root;
	}

	void iterate(T target, const std::vector<double>& vp_dist, unsigned int lower, unsigned int upper) {
		for (size_t i = lower; i < upper; ++i) {
			bool okay = true;
			for (size_t j = 0; j < num_vantage_points; ++j) {
				if ((vp_dist[j] > distToVP[items[i]][j] + tau) || (distToVP[items[i]][j] > vp_dist[j] + tau)) {
					okay = false;
					break;
				}
			}
			if (okay) {
				double dist = distance(data[items[i]], target);
				if (dist < tau) {
					tau = dist;
					closest_idx = items[i];
				}
				if (tau == 0) {
					return;
				}
			}
		}
	}

	void search(DirectoryNode* node, T target, const std::vector<double>& vp_dist, unsigned int actVPIdx) {
		if (node == nullptr) {
			return;
		}
		if (actVPIdx == num_vantage_points - 1) {
			if (vp_dist[actVPIdx] <= node->threshold + tau) { // search in left side
				iterate(target, vp_dist, node->s1, node->e1);
			}
			if (vp_dist[actVPIdx] + tau >= node->threshold) { // search in right side
				iterate(target, vp_dist, node->s2, node->e2);
			}
		} else { // more search in directory structure
			if (vp_dist[actVPIdx] <= node->threshold + tau) { // search in left side
				search(node->left.get(), target, vp_dist, actVPIdx + 1);
			}
			if (vp_dist[actVPIdx] + tau >= node->threshold) { // search in right side
				search(node->right.get(), target, vp_dist, actVPIdx + 1);
			}
		}
	}
};

