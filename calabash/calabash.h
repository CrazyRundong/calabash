#pragma once

#ifndef __CALABASH_CALABASH_H__
#define __CALABASH_CALABASH_H__

// optional: mkl library
// #define EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_set>
#include <cstdlib>

#ifdef _DEBUG_DIAMOND_
#include <iomanip>
#include <ctime>
#endif

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
    Edge(const int h, const int t, const double w) :
        head(h), tail(t), weight(w) {}
    ~Edge() {}
};

/*
 * Edges
 * instances share the same vector, hash storage
 */
class Edges {
public:
    Edges();
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

/*
 * Diamond
 * each Diamond instance gets unique state, while share same Edges
 */
class Diamond {
public:
    Diamond() {}
    Diamond(const Diamond &rhs) : _edgeMat(rhs._edgeMat), _scoreMat(rhs._scoreMat),
        _state(new std::unordered_set<int>(*rhs._state)), _edges(rhs._edges), _power(rhs._power) {}
    Diamond(std::istream&);
    Diamond(const Edges &edges, const int numNodes);
    Diamond(const std::vector<int> &states, const Edges &edges, const int numNodes);
    ~Diamond() {}

    double get_power() { return _power; }
    int get_num_nodes() { return _edges->num_nodes(); }
    std::string get_state_output();

    void set_state(const std::vector<int> &states);
    void set_state(const int idx);

private:
    // Eigen::Matrix uses deep-copy constructor
    Eigen::MatrixXd _edgeMat;
    Eigen::MatrixXd _scoreMat;
    std::unique_ptr<std::unordered_set<int>> _state;
    std::shared_ptr<Edges> _edges;
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

/*
 * Solver virtual class
 */
class Solver {
public:
    Solver() : _best_power(0.), _best_solution(nullptr) {}
    virtual ~Solver() {}
    virtual std::ostream& print_result(std::ostream&) = 0;
    virtual double solve(const int) = 0;
protected:
    double _best_power;
    std::unique_ptr<Diamond> _best_solution;
};

class RandomWalkingSolver : public Solver {
public:
    RandomWalkingSolver() = default;
    RandomWalkingSolver(const std::string &file);

    std::ostream& print_result(std::ostream&) override;
    double solve(const int steps) override;
private:
    std::shared_ptr<std::vector<Diamond>> _solver_list;
};

#endif // !__CALABASH_CALABASH_H__
