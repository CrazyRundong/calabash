#include <sstream>
#include "calabash.h"

Edges::Edges() : _edges(std::make_shared<std::vector<Edge>>()), _numNodes(0) { }
Edges::Edges(std::istream &ins) : _edges(std::make_shared<std::vector<Edge>>()), _numNodes(0) {
	Edge e;
	int h, t;
	double w;
	std::string line;
	ins >> _numNodes;
	while (std::getline(ins, line)) {
		if (!line.length())
			continue;
		std::stringstream tokens(line);
		tokens >> h >> t >> w;
		e = { h, t, w };
		_edges->push_back(e);
	}
}

Edges &Edges::push_back(const Edge &e) {
	this->_edges->push_back(e);
	return *this;
}

// Diamond methods

#define IN_GRAPH(e, s) ((e.head == 0 || s.find(e.head) != s.end()) && \
	s.find(e.tail) != s.end())

Diamond::Diamond(const Edges &edges, const int numNodes) :
	_edgeMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_scoreMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_edges(edges),
	_numNodes(numNodes) {
	set_init_state();
}

Diamond::Diamond(const std::vector<int> &states, const Edges &edges, const int numNodes) :
	_edgeMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_scoreMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_state(std::unordered_set<int>(states.begin(), states.end())),
	_edges(edges),
	_numNodes(numNodes) {
	for (auto const &edge : edges) {
		if (IN_GRAPH(edge, _state)) {
			_edgeMat(abs(edge.head), abs(edge.tail)) = edge.weight;
		}
	}
	update_power();
}

void Diamond::set_init_state() {
	std::unordered_set<int> states;
	Eigen::MatrixXd tempGraph, tempScore;
	double power_pos, power_neg;

	for (int i = 1; i < _numNodes + 1; ++i) {
		tempGraph.resize(i + 1, i + 1);
		tempScore.resize(i + 1, i + 1);

		std::unordered_set<int> state_pos(states);
		std::unordered_set<int> state_neg(states);
		state_pos.insert(i);
		state_neg.insert(-i);

		power_pos = calc_power(tempGraph, tempScore, state_pos, i);
		power_neg = calc_power(tempGraph, tempScore, state_neg, i);
		states = power_pos > power_neg ? state_pos : state_neg;
	}
	_state = states;
	_edgeMat = tempGraph;
	_scoreMat = tempScore;
	_power = std::max(power_pos, power_neg);
}

double Diamond::calc_power(Eigen::MatrixXd &graphMat,
	Eigen::MatrixXd &scoreMat,
	const std::unordered_set<int> &states,
	const int numNodes) {
	double power;
	graphMat.setZero();
	scoreMat.setZero();

	// TODO: maybe very slow here
	for (auto const &edge : _edges) {
		if (IN_GRAPH(edge, states)) {
			graphMat(abs(edge.head), abs(edge.tail)) = edge.weight;
		}
	}

	scoreMat = -1. * graphMat;
	for (int i = 0; i < numNodes + 1; ++i) {
		scoreMat(i, i) += graphMat.col(i).sum();
	}
	power = scoreMat.block(1, 1, numNodes, numNodes).determinant();

	return power;
}

void Diamond::set_state(const std::vector<int> &states) {
	_state.clear();
	for (auto const &state : states) {
		_state.insert(state);
	}
}

inline void Diamond::set_state(const int idx) {
	_state.erase(-idx);
	_state.insert(idx);
}

// functions to calculate power
// only called after state_update, edge_update
double Diamond::update_power() {
	_scoreMat = -1. * _edgeMat;
	for (int i = 0; i < _numNodes + 1; ++i) {
		_scoreMat(i, i) = _edgeMat.col(i).sum();
	}
	_power = _scoreMat.norm();

	return _power;
}

