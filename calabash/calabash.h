#pragma once

#ifndef __CALABASH_CALABASH_H__
#define __CALABASH_CALABASH_H__

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_set>
#include <cstdlib>

// macros to determin if an node state is in state list
// p, e: std::pair<int, double>(idx_of_node, weight_of_edge),
// which usually return from map h2t (head to tail) or t2h from Edges
#define PAIR_IN_GRAPH(p, s) (s.find(p.first) != s.end())
#define IN_GRAPH(e, s) ((e.head == 0 || s.find(e.head) != s.end()) && \
	s.find(e.tail) != s.end())

struct Edge {
	int head, tail;
	double weight;

	Edge() {}
	//Edge(const Edge &rhs) :
	//	head(rhs.head), tail(rhs.tail), weight(rhs.weight) {}
	Edge(const int h, const int t, const double w) :
		head(h), tail(t), weight(w) {}
	~Edge() {}
};

class Edges {
public:
	Edges();
	//Edges(const Edges &rhs): _edges(rhs._edges), _h2t(rhs._h2t), _t2h(rhs._t2h),
	//	_numNodes(rhs._numNodes) {}
	Edges(std::istream&);
	~Edges() {}

	// query methods
	size_t num_edges() { return _edges->size(); }
	int num_nodes() const { return _numNodes; }
	auto query_tail(const int h);
	auto query_head(const int t);

	// write methods
	Edges &push_back(const Edge &e);

	// const iterators
	auto begin() const { return _edges->begin(); }
	auto end() const { return _edges->end(); }

private:
	std::shared_ptr<std::vector<Edge>> _edges;
	// use secondary map to store [h->t](w)
	std::shared_ptr<std::map<int, std::map<int, double>>> _h2t;
	std::shared_ptr<std::map<int, std::map<int, double>>> _t2h;
	int _numNodes;
};

class Diamond {
public:
	Diamond() {}
	//Diamond(const Diamond &rhs) : _edgeMat(rhs._edgeMat), _scoreMat(rhs._scoreMat),
	//	_state(rhs._state), _edges(rhs._edges), _numNodes(rhs._numNodes), _power(rhs._power) {}
	Diamond(std::istream&);
	Diamond(const Edges &edges, const int numNodes);
	Diamond(const std::vector<int> &states, const Edges &edges, const int numNodes);
	~Diamond() {}

	double get_power() { return _power; }
	int get_num_nodes() { return _numNodes; }
	std::string get_state_output(std::ostream &os);

	void set_state(const std::vector<int> &states);
	void set_state(const int idx);

private:
	Eigen::MatrixXd _edgeMat;
	Eigen::MatrixXd _scoreMat;
	std::shared_ptr<std::unordered_set<int>> _state;
	Edges _edges;
	int _numNodes;
	double _power;

	double calc_power(Eigen::MatrixXd &graphMat,
		Eigen::MatrixXd &scoreMat,
		const std::unordered_set<int> &states,
		const int numNodes);
	double calc_power(const Eigen::MatrixXd &graphMat, Eigen::MatrixXd &scoreMat, const int numNodes);
	void set_init_state();

	// update methods
	void update_edge_mat(const int state);
};

#endif // !__CALABASH_CALABASH_H__
