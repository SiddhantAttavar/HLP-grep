/**
 * @file cost_model.hpp
 * @brief Cost model for edit distance computations.
 *
 * CostModel holds the costs of the basic edit operations (insertion,
 * deletion and match/mismatch substitution). Consumers such as POAGraph,
 * BinaryLifter, NaiveSolver and Solver store a CostModel& reference, so a
 * referenced model must outlive every object built from it.
 * DEFAULT_COST_MODEL provides a permanent fallback for call sites that
 * just want unit costs.
 */
#pragma once

#include <stdexcept>

namespace hlp_grep {

/**
 * @brief Edit operation costs.
 *
 * Consuming equal characters costs @p match; any distinct pair costs
 * @p mismatch. Defaults to the unit-cost model (match 0, mismatch 1,
 * insertion/deletion 1).
 */
class CostModel {
public:
	/**
	 * @brief Constructs the model with the given costs.
	 *
	 * Rejects any model that does not satisfy the structural precondition
	 * behind the argmin-staircase composition used by DistMatrix and
	 * BinaryLifter (see DistMatrix::min_plus_product):
	 *   match == 0,
	 *   ins >= 0,
	 *   del >= 0,
	 *   ins + del >= mismatch
	 * (the last implies mismatch >= 0, so match == 0 makes every
	 * substitution cost at least the match cost).
	 *
	 * @param ins      Insertion cost (default 1).
	 * @param del      Deletion cost (default 1).
	 * @param match    Cost of consuming equal characters (default 0).
	 * @param mismatch Cost of consuming distinct characters (default 1).
	 * @throws std::invalid_argument if any precondition is violated.
	 */
	CostModel(int ins = 1, int del = 1, int match = 0, int mismatch = 1)
	    : match(match), mismatch(mismatch), ins_cost(ins), del_cost(del) {
		if (match != 0)
			throw std::invalid_argument(
			    "CostModel: match cost must be 0");
		if (ins_cost < 0)
			throw std::invalid_argument(
			    "CostModel: insertion cost must be nonnegative");
		if (del_cost < 0)
			throw std::invalid_argument(
			    "CostModel: deletion cost must be nonnegative");
		if (ins_cost + del_cost < mismatch)
			throw std::invalid_argument(
			    "CostModel: ins + del must be >= mismatch");
	}

	/** @brief Cost of inserting a character. */
	int ins() const {
		return ins_cost;
	}

	/** @brief Cost of deleting a character. */
	int del() const {
		return del_cost;
	}

	/**
	 * @brief Returns the cost of consuming @p b with @p g.
	 *
	 * @param g Graph character being consumed.
	 * @param b Character consumed against it.
	 * @return The substitution cost of matching @p g with @p b.
	 */
	int consume(char g, char b) const {
		return g == b ? match : mismatch;
	}

	int match;    ///< Cost of consuming equal characters.
	int mismatch; ///< Cost of consuming distinct characters.

private:
	int ins_cost; ///< Cost of inserting a character.
	int del_cost; ///< Cost of deleting a character.
};

/// Default unit-cost model: reference fallback for callers without their
/// own model object.
inline const CostModel DEFAULT_COST_MODEL{};

} // namespace hlp_grep
