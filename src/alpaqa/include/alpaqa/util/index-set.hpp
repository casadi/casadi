#pragma once

#include <alpaqa/config/config.hpp>

namespace alpaqa::detail {

template <Config Conf>
struct IndexSet {
    USING_ALPAQA_CONFIG(Conf);

    IndexSet(length_t N, length_t n) : N{N}, n{n}, storage{N + N * n} {}

    length_t N;
    length_t n;
    indexvec storage;

    auto sizes() { return storage.segment(0, N); }
    auto sizes() const { return storage.segment(0, N); }
    auto indices() { return storage.segment(N, N * n); }
    auto indices() const { return storage.segment(N, N * n); }

    crindexvec indices(index_t i) const {
        length_t nJ = sizes()(i);
        return indices().segment(n * i, nJ);
    }

    crindexvec compl_indices(index_t i) const {
        length_t nJ = sizes()(i);
        length_t nK = n - nJ;
        return indices().segment(n * i + nJ, nK);
    }

    template <class F>
    void update(const F &condition) {
        // Evaluate the condition for all indices in the given 'time_step',
        // append the indices evaluating to true to 'out', and return the number
        // of indices that were appended.
        auto build_Jt = [&](index_t time_step, index_t *out) {
            index_t j = 0; // index into the array of inactive indices
            for (index_t c = 0; c < n; ++c) { // components within time step
                if (condition(time_step, c))  // if the component is active,
                    out[j++] = c; // append the index of this component to J
            }
            return j; // return the number of elements in J
        };
        // Compute the complement of the index set defined by the range
        // 'in' and 'n_in'. Append the indices not in the input to 'out'.
        auto complement = [&](const index_t *in, length_t n_in, index_t *out) {
            length_t c = 0; // components within time step
            length_t k = 0; // index into the array of active indices
            // iterate over the array with indices 'in'
            for (index_t i = 0; i < n_in; ++i) { //
                index_t j = in[i];
                for (; c < j; ++c) // for all indices not in J
                    out[k++] = c;  // append the index of this component
                ++c;               // skip indices in J, i.e. c == j
            }
            // add final indices not in J
            for (; c < n; ++c)
                out[k++] = c;
        };

        auto sizes    = this->sizes();
        auto *indices = this->indices().data();
        for (index_t t = 0; t < N; ++t) { // time steps
            // Generate the set of inactive indices, J
            index_t num_J = build_Jt(t, indices);
            sizes(t)      = num_J; // save number of inactive indices for later
            // Generate the complement, the set of active indices, K
            complement(indices, num_J, indices + num_J);
            // Prepare for next time step
            indices += n;
        }
    }
};

} // namespace alpaqa::detail