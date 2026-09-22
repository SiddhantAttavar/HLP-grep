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
	 * @param ins      Insertion cost (default 1).
	 * @param del      Deletion cost (default 1).
	 * @param match    Cost of consuming equal characters (default 0).
	 * @param mismatch Cost of consuming distinct characters (default 1).
	 */
	CostModel(int ins = 1, int del = 1, int match = 0, int mismatch = 1)
	    : match(match), mismatch(mismatch), ins_cost(ins), del_cost(del) {
		monge = match >= 0 && mismatch >= 0 && ins_cost >= 0 &&
		        del_cost >= 0;
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

	/**
	 * @brief Reports whether the model produces Monge chain matrices.
	 *
	 * The single-edge chain blocks decompose into a scaled window minimum
	 * plus a linear term, which is Monge for any nonnegative cost
	 * assignment; negative costs are reported as non-Monge. The value is
	 * precomputed in the constructor.
	 */
	bool is_monge() const {
		return monge;
	}

	int match;    ///< Cost of consuming equal characters.
	int mismatch; ///< Cost of consuming distinct characters.

private:
	int ins_cost; ///< Cost of inserting a character.
	int del_cost; ///< Cost of deleting a character.
	bool monge = false; ///< Precomputed is_monge() value.
};

/// Default unit-cost model: reference fallback for callers without their
/// own model object.
inline const CostModel DEFAULT_COST_MODEL{};

} // namespace hlp_grep
