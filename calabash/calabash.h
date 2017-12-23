#pragma once

#ifndef __CALABASH_CALABASH_H__
#define __CALABASH_CALABASH_H__

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <unordered_set>
#include <cstdlib>

struct Edge {
	int head, tail;
	double weight;

	Edge() {}
	Edge(const int h, const int t, const double w) :
		head(h), tail(t), weight(w) {}
	~Edge() {}
};

class Edges {
public:
	Edges();
	Edges(const Edges &rhs): _edges(rhs._edges), _numNodes(rhs._numNodes) {}
	Edges(std::istream&);
	~Edges() {}

	size_t num_edges() { return _edges->size(); }
	int num_nodes() const { return _numNodes; }
	Edges &push_back(const Edge &e);

	// const iterators
	auto begin() const { return _edges->begin(); }
	auto end() const { return _edges->end(); }

private:
	std::shared_ptr<std::vector<Edge>> _edges;
	int _numNodes;
};

class Diamond {
public:
	Diamond() {};
	Diamond(const Edges &edges, const int numNodes);
	Diamond(const std::vector<int> &states, const Edges &edges, const int numNodes);
	~Diamond() {}

	const double get_power() { return _power; }
	const int get_state(const int i) {
		int pos_i = i, neg_i = -i;
		if (_state.find(pos_i) != _state.end()) {
			return 1;
		} else if (_state.find(neg_i) != _state.end()) {
			return -1;
		} else {
			return 0;
		}
	}
	std::string &get_state_output();

private:
	Eigen::MatrixXd _edgeMat;
	Eigen::MatrixXd _scoreMat;
	std::unordered_set<int> _state;
	Edges _edges;
	int _numNodes;
	double _power;

	double calc_power(Eigen::MatrixXd &graphMat,
		Eigen::MatrixXd &scoreMat,
		const std::unordered_set<int> &states,
		const int numNodes);
	double update_power();
	void set_init_state();
	void set_state(const std::vector<int> &states);
	void set_state(const int idx);
};


#endif // !__CALABASH_CALABASH_H__
