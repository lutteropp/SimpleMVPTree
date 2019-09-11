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
			_data(data), _tau(std::numeric_limits<double>::max()), distance(dist), _root(NULL) {
		_num_vantage_points = ceil(log(data.size() / 2));
		_vp_indices.resize(_num_vantage_points);
		_distToVP.resize(data.size());
		for (size_t i = 0; i < _distToVP.size(); ++i) {
			_distToVP[i].resize(_num_vantage_points);
		}
		_items.resize(data.size());
		for (size_t i = 0; i < _items.size(); ++i) {
			_items[i] = i;
		}

		_root = buildFromPoints(0, _items.size());
	}

	~SimpleMVPTree() {
		delete _root;
	}

	double search_mindist(const T& target) {
		_tau = std::numeric_limits<double>::max();

		std::vector<double> vp_dist(_num_vantage_points);
		for (size_t i = 0; i < _num_vantage_points; ++i) {
			vp_dist[i] = distance(_data[_vp_indices[i]], target);
			if (vp_dist[i] < _tau) {
				_tau = vp_dist[i];
				_closest_idx = _vp_indices[i];
			}
			if (_tau == 0) {
				return 0;
			}
		}
		search(_root, target, vp_dist, 0);

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

		return _tau;
	}

	T search_mindist_elem(const T& target) {
		search_mindist(target);
		return _data[_closest_idx];
	}

private:
	unsigned int _num_vantage_points;
	const std::vector<T>& _data;
	std::vector<unsigned int> _items;
	double _tau;
	size_t _closest_idx;
	std::vector<unsigned int> _vp_indices;
	std::vector<std::vector<double> > _distToVP;
	std::function<double(const T&, const T&)> distance;

	typedef struct DirectoryNode {
		DirectoryNode* left;
		DirectoryNode* right;
		unsigned int s1, e1;
		unsigned int s2, e2;
		double threshold; // median distance
	} DirectoryNode;

	DirectoryNode* _root;

	struct DistanceComparator {
		const std::vector<std::vector<double> >& _distToVP;
		unsigned int _vpIdx;

		DistanceComparator(const std::vector<std::vector<double> >& distToVP, unsigned int vpIdx) :
				_distToVP(distToVP), _vpIdx(vpIdx) {
		}
		bool operator()(const unsigned int& a, const unsigned int& b) {
			return _distToVP[a][_vpIdx] < _distToVP[b][_vpIdx];
		}
	};

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
			double dist = distance(_data[_items[cand]], _data[_items[i]]);
			if (dist > maxDist) {
				maxDist = dist;
				vpPoint = i;
			}
		}
		return vpPoint;
	}

	DirectoryNode* makeDirectory(unsigned int lower, unsigned int upper, unsigned int actVPIndex) {
		DirectoryNode* node = new DirectoryNode();
		if (lower >= upper) {
			return node;
		}

		unsigned int median = (upper + lower) / 2;
		std::nth_element(_items.begin() + lower, _items.begin() + median, _items.begin() + upper,
				DistanceComparator(_distToVP, actVPIndex));
		node->s1 = lower;
		node->e1 = median;
		node->s2 = median;
		node->e2 = upper;
		node->threshold = _distToVP[_items[median]][actVPIndex];

		if (actVPIndex < _num_vantage_points - 1) {
			node->left = makeDirectory(node->s1, node->e1, actVPIndex + 1);
			node->right = makeDirectory(node->s2, node->e2, actVPIndex + 1);
		}

		return node;
	}

	DirectoryNode* buildFromPoints(unsigned int lower, unsigned int upper) {
		DirectoryNode* root = new DirectoryNode();
		if (upper == lower) {
			return root;
		}
		unsigned int actVPIndex = 0;

		unsigned int maxDist = 0;
		unsigned int maxDistIdx = 0;

		// first Vantage point is a bit different from the rest:
		unsigned int vp = findVantagePoint(0, upper);

		_vp_indices[0] = _items[vp];
		std::swap(_items[lower], _items[vp]);
		for (size_t j = lower + 1; j < upper; ++j) {
			double dist = distance(_data[_items[j]], _data[_items[lower]]);
			_distToVP[_items[j]][0] = dist;
			if (dist > maxDist) {
				maxDist = dist;
				maxDistIdx = _items[j];
			}
		}
		// select remaining Vantage points
		for (size_t i = 1; i < _num_vantage_points; ++i) {
			unsigned int vp = maxDistIdx;
			_vp_indices[i] = _items[vp];
			std::swap(_items[lower + i], _items[vp]);
			maxDist = 0;
			maxDistIdx = 0;
			for (size_t j = lower + i + 1; j < upper; ++j) {
				double dist = distance(_data[_items[j]], _data[_items[lower + i]]);
				_distToVP[_items[j]][i] = dist;
				// look at the minimum distance from vantage points so far
				for (size_t k = 1; k < i; ++k) {
					dist = std::min(dist, _distToVP[_items[j]][k]);
				}
				if (dist > maxDist) {
					maxDist = dist;
					maxDistIdx = _items[j];
				}
			}
		}
		// reorder the items and store it in the directory.
		unsigned int median = (upper + lower + _num_vantage_points) / 2;
		// partition around the median distance from first VP
		std::nth_element(_items.begin() + lower + _num_vantage_points, _items.begin() + median, _items.begin() + upper,
				DistanceComparator(_distToVP, actVPIndex));
		root->s1 = lower + _num_vantage_points;
		root->e1 = median;
		root->s2 = median;
		root->e2 = upper;
		root->threshold = _distToVP[_items[median]][actVPIndex];

		if (actVPIndex < _num_vantage_points - 1) {
			root->left = makeDirectory(root->s1, root->e1, actVPIndex + 1);
			root->right = makeDirectory(root->s2, root->e2, actVPIndex + 1);
		}
		return root;
	}

	void iterate(T target, const std::vector<double>& vp_dist, unsigned int lower, unsigned int upper) {
		for (size_t i = lower; i < upper; ++i) {
			bool okay = true;
			for (size_t j = 0; j < _num_vantage_points; ++j) {
				if ((vp_dist[j] > _distToVP[_items[i]][j] + _tau) || (_distToVP[_items[i]][j] > vp_dist[j] + _tau)) {
					okay = false;
					break;
				}
			}
			if (okay) {
				double dist = distance(_data[_items[i]], target);
				if (dist < _tau) {
					_tau = dist;
					_closest_idx = _items[i];
				}
				if (_tau == 0) {
					return;
				}
			}
		}
	}

	void search(DirectoryNode* node, T target, const std::vector<double>& vp_dist, unsigned int actVPIdx) {
		if (node == NULL) {
			return;
		}
		if (actVPIdx == _num_vantage_points - 1) {
			if (vp_dist[actVPIdx] <= node->threshold + _tau) { // search in left side
				iterate(target, vp_dist, node->s1, node->e1);
			}
			if (vp_dist[actVPIdx] + _tau >= node->threshold) { // search in right side
				iterate(target, vp_dist, node->s2, node->e2);
			}
		} else { // more search in directory structure
			if (vp_dist[actVPIdx] <= node->threshold + _tau) { // search in left side
				search(node->left, target, vp_dist, actVPIdx + 1);
			}
			if (vp_dist[actVPIdx] + _tau >= node->threshold) { // search in right side
				search(node->right, target, vp_dist, actVPIdx + 1);
			}
		}
	}
};

