/**
 * @file cost_model.hpp
 * @brief Cost models for edit distance computations.
 *
 * CostModel is an abstract interface with virtual clone(); the two
 * concrete models (UnitCostModel and MatrixCostModel) derive from it and
 * override consume() and is_monge(). Solver owns its model through a
 * cloned std::unique_ptr (the only unique_ptr user in the codebase);
 * POAGraph, BinaryLifter and NaiveSolver still hold CostModel& references,
 * so a referenced model object must outlive every object built from it.
 * DEFAULT_COST_MODEL provides a permanent fallback for call sites that
 * just want unit costs.
 */
#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hlp_grep {

/**
 * @brief Abstract interface for edit operation costs.
 *
 * Only the concrete subclasses are constructible by callers; every
 * consumer stores a CostModel& reference.
 */
class CostModel {
public:
	virtual ~CostModel() = default;

	/**
	 * @brief Returns an owned deep copy of this model under the abstract
	 *        interface (prototype pattern; one line per subclass).
	 */
	virtual std::unique_ptr<CostModel> clone() const = 0;

	/**
	 * @brief Cost of inserting @p b into a sequence.
	 *
	 * Only MatrixCostModel distinguishes characters; UnitCostModel
	 * returns its scalar cost for every character.
	 */
	virtual int ins(char b) const = 0;

	/**
	 * @brief Cost of deleting @p a from a sequence.
	 *
	 * Only MatrixCostModel distinguishes characters; UnitCostModel
	 * returns its scalar cost for every character.
	 */
	virtual int del(char a) const = 0;

	/**
	 * @brief Returns the cost of consuming @p b with @p g.
	 *
	 * @param g Graph character being consumed.
	 * @param b Character consumed against it.
	 * @return The substitution cost of matching @p g with @p b.
	 */
	virtual int consume(char g, char b) const = 0;

	/**
	 * @brief Reports whether the model produces Monge chain matrices.
	 *
	 * A cost check without building any DP tables (see the concrete
	 * implementations).
	 */
	virtual bool is_monge() const = 0;

protected:
	CostModel() = default;
};

/**
 * @brief Unit-style cost model without a substitution matrix.
 *
 * Consuming equal characters costs @p match; any distinct pair costs
 * @p mismatch. Defaults to the unit-cost model (match 0, mismatch 1,
 * insertion/deletion 1).
 */
class UnitCostModel : public CostModel {
public:
	/**
	 * @brief Constructs the model with the given costs.
	 *
	 * @param ins      Insertion cost (default 1).
	 * @param del      Deletion cost (default 1).
	 * @param match    Cost of consuming equal characters (default 0).
	 * @param mismatch Cost of consuming distinct characters (default 1).
	 */
	UnitCostModel(int ins = 1, int del = 1, int match = 0, int mismatch = 1)
	    : ins_cost(ins), del_cost(del), match(match), mismatch(mismatch) {
		monge = match >= 0 && mismatch >= 0 && ins_cost >= 0 &&
		        del_cost >= 0;
	}

	int ins(char) const override {
		return ins_cost;
	}

	int del(char) const override {
		return del_cost;
	}

	int consume(char g, char b) const override {
		return g == b ? match : mismatch;
	}

	/**
	 * @brief Simple O(1) Monge check, precomputed in the constructor.
	 *
	 * The single-edge chain blocks decompose into a scaled window minimum
	 * plus a linear term, which is Monge for any nonnegative cost
	 * assignment; negative costs are reported as non-Monge.
	 */
	bool is_monge() const override {
		return monge;
	}

	std::unique_ptr<CostModel> clone() const override {
		return std::make_unique<UnitCostModel>(*this);
	}

	int match;    ///< Cost of consuming equal characters.
	int mismatch; ///< Cost of consuming distinct characters.

private:
	int ins_cost; ///< Cost of inserting a character.
	int del_cost; ///< Cost of deleting a character.
	bool monge = false; ///< Precomputed is_monge() value.
};

/// Default unit-cost model: reference fallback for callers without their
/// own model object (see CostModel).
inline const UnitCostModel DEFAULT_COST_MODEL{};

/**
 * @brief Substitution-matrix cost model.
 *
 * Substitution costs come from an explicit |alphabet| x |alphabet| matrix;
 * characters outside the alphabet fall back to index 0. Insertion and
 * deletion costs are set per character; an empty cost vector defaults to
 * cost 1 for every character.
 */
class MatrixCostModel : public CostModel {
public:
	/**
	 * @brief Constructs the model with the given substitution matrix and
	 *        per-character insertion/deletion costs.
	 *
	 * @param alphabet Characters of the alphabet, defining the row/column
	 *                 order of @p matrix.
	 * @param matrix Substitution cost matrix: `matrix[i][j]` is the cost of
	 *               matching `alphabet[i]` with `alphabet[j]`. Characters not
	 *               present in @p alphabet are looked up at index 0.
	 *               Defaults to the unit-cost model over the DNA alphabet
	 *               (0 on the diagonal, 1 elsewhere).
	 * @param ins_costs Per-character insertion costs, in alphabet order;
	 *                  empty to use cost 1 for every character.
	 * @param del_costs Per-character deletion costs, in alphabet order;
	 *                  empty to use cost 1 for every character.
	 */
	MatrixCostModel(std::string alphabet = "AGCT",
	                std::vector<std::vector<int>> matrix = {},
	                std::vector<int> ins_costs = {},
	                std::vector<int> del_costs = {}) {
		this->alphabet = std::move(alphabet);
		if (matrix.empty()) {
			this->matrix.assign(this->alphabet.size(),
			                    std::vector<int>(this->alphabet.size(), 1));
			for (std::size_t i = 0; i < this->alphabet.size(); ++i)
				this->matrix[i][i] = 0;
		} else {
			this->matrix = std::move(matrix);
		}
		if (ins_costs.empty())
			ins_costs.assign(this->alphabet.size(), 1);
		if (del_costs.empty())
			del_costs.assign(this->alphabet.size(), 1);
		if (ins_costs.size() != this->alphabet.size() ||
		    del_costs.size() != this->alphabet.size())
			throw std::invalid_argument(
			    "MatrixCostModel: per-character cost vectors must have "
			    "alphabet size");
		this->ins_costs = std::move(ins_costs);
		this->del_costs = std::move(del_costs);
		for (std::size_t i = 0; i < this->alphabet.size(); ++i)
			index[static_cast<unsigned char>(this->alphabet[i])] = i;
		monge = !has_negative_cycle();
	}

	int ins(char b) const override {
		return ins_costs[index_of(b)];
	}

	int del(char a) const override {
		return del_costs[index_of(a)];
	}

	int consume(char g, char b) const override {
		return matrix[index_of(g)][index_of(b)];
	}

	/**
	 * @brief True if no negative-cost cycle exists in the scoring digraph.
	 *
	 * Runs the closure over the alphabet-plus-gap edit digraph and checks
	 * whether any walk can return to its start node at negative total
	 * cost (including single-character self-loops, so a negative diagonal
	 * entry counts as a cycle). No edit-distance DP tables are built.
	 */
	bool has_negative_cycle() const {
		const std::vector<std::vector<int>> dist = closure();
		for (std::size_t i = 0; i <= alphabet.size(); ++i)
			if (dist[i][i] < 0)
				return true;
		return false;
	}

	/**
	 * @brief Replaces every cost with its shortest-path closure.
	 *
	 * Builds the scoring digraph over the alphabet plus a gap node (rows x
	 * gap = del, gap x columns = ins, everything else = consume) and
	 * replaces each of the three cost tables with the minimum cost of a
	 * walk between its endpoints: matrix[a][b] costs a walk a -> b,
	 * ins[b] a walk gap -> b and del[a] a walk a -> gap. After the call no
	 * two-step edit chain is cheaper than the direct edit, which is the
	 * premise behind the Monge guarantees of the solver's chain blocks.
	 *
	 * @throws std::invalid_argument if has_negative_cycle() reports a
	 *         negative cost cycle.
	 */
	void make_monge() {
		if (has_negative_cycle())
			throw std::invalid_argument(
			    "MatrixCostModel::make_monge: negative-cost cycle");
		const std::vector<std::vector<int>> dist = closure();
		for (std::size_t i = 0; i < alphabet.size(); ++i)
			for (std::size_t j = 0; j < alphabet.size(); ++j)
				matrix[i][j] = dist[i][j];
		for (std::size_t i = 0; i < alphabet.size(); ++i) {
			del_costs[i] = dist[i][alphabet.size()];
			ins_costs[i] = dist[alphabet.size()][i];
		}
		monge = !has_negative_cycle();
	}

	/**
	 * @brief True unless the scoring digraph has a negative-cost cycle.
	 *
	 * Such a cycle is the only way a model can produce non-Monge chain
	 * blocks: with cycle-free costs the block recurrence decomposes into
	 * a linear term plus a running window minimum, which is Monge, and
	 * min-plus compositions preserve that. The value is precomputed in
	 * the constructor and refreshed by make_monge().
	 */
	bool is_monge() const override {
		return monge;
	}

	std::unique_ptr<CostModel> clone() const override {
		return std::make_unique<MatrixCostModel>(*this);
	}

private:
	/**
	 * @brief Floyd-Warshall closure over the scoring digraph.
	 *
	 * Nodes are the alphabet characters plus one gap node at index
	 * |alphabet|; direct edges are sub(a, b), del(a) (a -> gap) and
	 * ins(b) (gap -> b). Char diagonals start at the direct consume cost
	 * (so the result keeps the direct value unless a walk beats it); the
	 * gap node starts unreachable from itself, exposing gap cycles as
	 * alternating insert/delete walks.
	 */
	std::vector<std::vector<int>> closure() const {
		const std::size_t sigma = alphabet.size();
		constexpr int INF = std::numeric_limits<int>::max() / 2;
		std::vector<std::vector<int>> dist(
		    sigma + 1, std::vector<int>(sigma + 1, INF));
		for (std::size_t i = 0; i < sigma; ++i)
			for (std::size_t j = 0; j < sigma; ++j)
				dist[i][j] = matrix[i][j];
		for (std::size_t i = 0; i < sigma; ++i)
			dist[i][sigma] = del_costs[i]; // walk i -> gap: delete i
		for (std::size_t j = 0; j < sigma; ++j)
			dist[sigma][j] = ins_costs[j]; // walk gap -> j: insert j
		// No direct gap -> gap edge; walks through a character are allowed.
		for (std::size_t mid = 0; mid <= sigma; ++mid) {
			if (dist[mid][mid] == INF)
				continue;
			for (std::size_t i = 0; i <= sigma; ++i) {
				if (dist[i][mid] == INF)
					continue;
				for (std::size_t j = 0; j <= sigma; ++j)
					if (dist[mid][j] != INF)
						dist[i][j] = std::min(
						    dist[i][j], dist[i][mid] + dist[mid][j]);
			}
		}
		return dist;
	}

	std::size_t index_of(char c) const {
		return index[static_cast<unsigned char>(c)];
	}

	std::string alphabet; ///< Alphabet, defining matrix order.
	std::vector<std::vector<int>> matrix; ///< Substitution cost matrix.
	std::vector<int> ins_costs; ///< Per-character insertion costs.
	std::vector<int> del_costs; ///< Per-character deletion costs.
	std::size_t index[256] = {}; ///< char -> alphabet index.
	bool monge = false; ///< Precomputed is_monge() value.
};

} // namespace hlp_grep
