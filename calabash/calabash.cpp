#include "calabash.h"

/*
 * Edge methods
 */

Edges::Edges() :
	_edges(std::make_shared<std::vector<Edge>>()),
	_h2t(std::make_shared<std::map<int, std::map<int, double>>>()),
	_t2h(std::make_shared<std::map<int, std::map<int, double>>>()),
	_numNodes(0) {}

Edges::Edges(std::istream &ins) :
	_edges(std::make_shared<std::vector<Edge>>()),
	_h2t(std::make_shared<std::map<int, std::map<int, double>>>()),
	_t2h(std::make_shared<std::map<int, std::map<int, double>>>()),
	_numNodes(0) {
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
		// not very modern, but works
		(*_h2t)[h].insert({ t, w });
		(*_t2h)[t].insert({ h, w });
	}
}

inline auto Edges::query_tail(const int h) {
	return _h2t->at(h);
}

inline auto Edges::query_head(const int t) {
	return _t2h->at(t);
}

Edges &Edges::push_back(const Edge &e) {
	this->_edges->push_back(e);
	(*_h2t)[e.head].insert({ e.tail, e.weight });
	(*_t2h)[e.tail].insert({ e.head, e.weight });
	return *this;
}

/*
 * Diamond methods
 */

Diamond::Diamond(std::istream &ism) :
	_edges(std::make_shared<Edges>(ism)) {
	_edgeMat = Eigen::MatrixXd::Zero(_edges->num_nodes() + 1, _edges->num_nodes() + 1);
	_scoreMat = Eigen::MatrixXd::Zero(_edges->num_nodes() + 1, _edges->num_nodes() + 1);
	set_init_state();
}

Diamond::Diamond(const Edges &edges, const int numNodes) :
	_edges(std::make_shared<Edges>(edges)) {
	_edgeMat = Eigen::MatrixXd::Zero(_edges->num_nodes() + 1, _edges->num_nodes() + 1);
	_scoreMat = Eigen::MatrixXd::Zero(_edges->num_nodes() + 1, _edges->num_nodes() + 1);
	set_init_state();
}

Diamond::Diamond(const std::vector<int> &states, const Edges &edges, const int numNodes) :
	_edgeMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_scoreMat(Eigen::MatrixXd::Zero(numNodes + 1, numNodes + 1)),
	_state(new std::unordered_set<int>(states.begin(), states.end())),
	_edges(std::make_shared<Edges>(edges)) {
	for (auto const &edge : edges) {
		if (IN_GRAPH(edge, (*_state))) {
			_edgeMat(abs(edge.head), abs(edge.tail)) = edge.weight;
		}
	}
	_power = calc_power(_edgeMat, _scoreMat, _edges->num_nodes());
}

std::string Diamond::get_state_output() {
	std::ostringstream fmtstr;
	for (int i = 1; i < get_num_nodes() + 1; ++i) {
		if (_state->find(i) != _state->end()) {
			fmtstr << std::showpos << i << " ";
		} else {
			fmtstr << (-i) << " ";
		}
	}

	return fmtstr.str();
}

void Diamond::set_state(const std::vector<int> &states) {
	_state->clear();
	for (auto const &state : states) {
		_state->insert(state);
		update_edge_mat(state);
	}
	_power = calc_power(_edgeMat, _scoreMat, _edges->num_nodes());
}

void Diamond::set_state(const int state) {
	_state->erase(-state);
	_state->insert(state);
	update_edge_mat(state);
	_power = calc_power(_edgeMat, _scoreMat, _edges->num_nodes());
}

// private methods:

void Diamond::set_init_state() {
	std::unordered_set<int> states;
	Eigen::MatrixXd tempGraph, tempScore;
	double power_pos, power_neg;

	for (int i = 1; i < _edges->num_nodes() + 1; ++i) {
		tempGraph.resize(i + 1, i + 1);
		tempScore.resize(i + 1, i + 1);
		tempGraph.setZero();
		tempScore.setZero();

		std::unordered_set<int> state_pos(states);
		std::unordered_set<int> state_neg(states);
		state_pos.insert(i);
		state_neg.insert(-i);

		power_pos = calc_power(tempGraph, tempScore, state_pos, i);
		power_neg = calc_power(tempGraph, tempScore, state_neg, i);
		states = power_pos > power_neg ? state_pos : state_neg;
	}
	_state.reset(new std::unordered_set<int>(states));
	_edgeMat = tempGraph;
	_scoreMat = tempScore;
	_power = std::max(power_pos, power_neg);
}

double Diamond::calc_power(Eigen::MatrixXd &graphMat,
	Eigen::MatrixXd &scoreMat,
	const std::unordered_set<int> &states,
	const int numNodes) {
	double power;

	// a faster graphMat building process
	// tairPair: (tailIdx, weight)
	for (auto const &tailPair : _edges->query_tail(0)) {
		if (PAIR_IN_GRAPH(tailPair, states)) {
			graphMat(0, std::abs(tailPair.first)) = tailPair.second;
		}
	}
	for (auto const head : states) {
		for (auto const &tailPair : _edges->query_tail(head)) {
			if (PAIR_IN_GRAPH(tailPair, states)) {
				graphMat(std::abs(head), std::abs(tailPair.first)) = tailPair.second;
			}
		}
	}
	power = calc_power(graphMat, scoreMat, numNodes);

	return power;
}

// functions to calculate power
// only called after _edgeMat modified
double Diamond::calc_power(const Eigen::MatrixXd &graphMat, Eigen::MatrixXd &scoreMat, const int numNodes) {
	double power;
	scoreMat = -1. * graphMat;
	for (int i = 0; i < numNodes + 1; ++i) {
		scoreMat(i, i) += graphMat.col(i).sum();
	}
	power = scoreMat.block(1, 1, numNodes, numNodes).determinant();

	return power;
}

void Diamond::update_edge_mat(const int state) {
	int idx = std::abs(state);
	// update edges starts with state
	for (auto const &tailPair : _edges->query_tail(state)) {
		if (PAIR_IN_GRAPH(tailPair, (*_state))) {
			_edgeMat(idx, std::abs(tailPair.first)) = tailPair.second;
		}
	}
	// update edges end with state
	for (auto const &headPair : _edges->query_head(state)) {
		if (PAIR_IN_GRAPH(headPair, (*_state))) {
			_edgeMat(std::abs(headPair.first), idx) = headPair.second;
		}
	}
}

/*
 * Solver classes
 */

// RandomWalkingSolver
RandomWalkingSolver::RandomWalkingSolver(const std::string &file) :
	_solver_list(std::make_shared<std::vector<Diamond>>()) {
	std::ifstream graphFile(file);
	std::shared_ptr<Diamond> init_solution = std::make_shared<Diamond>(graphFile);
	_solver_list->push_back(*init_solution);
	_best_power = _solver_list->back().get_power();
	_best_solution.reset(new Diamond(_solver_list->back()));
}

double RandomWalkingSolver::solve(const int steps) {
	double solved_power = 0.;
	double step_best_power = -1.;
	int num_nodes = _best_solution->get_num_nodes();
	std::unique_ptr<Diamond> step_best(new Diamond(*_best_solution));

	// TODO: parallelize this!
	for (int t = 0; t < steps; ++t) {
		_solver_list->clear();
		for (int i = 1; i < num_nodes + 1; ++i) {
			_solver_list->push_back(Diamond(*step_best));
			_solver_list->back().set_state(i);
			if (_solver_list->back().get_power() > step_best_power) {
				step_best_power = _solver_list->back().get_power();
				step_best.reset(new Diamond(_solver_list->back()));
			}
		}
		solved_power = step_best_power;
	}

	_best_power = solved_power;
	_best_solution.reset(step_best.release());

	return solved_power;
}

std::ostream& RandomWalkingSolver::print_result(std::ostream &iosm) {
	iosm << _best_solution->get_state_output();
	return iosm;
}
