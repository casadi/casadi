#pragma once

#include <alpaqa/config/config.hpp>
#include <span>

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

    // Compute the complement of the index set `in`. Write the indices not in
    // the input to 'out'.
    static void compute_complement(std::span<const index_t> in,
                                   std::span<index_t> out) {
        compute_complement(in, out.data(),
                           static_cast<length_t>(in.size() + out.size()));
    }
    static void compute_complement(std::span<const index_t> in,
                                   std::span<index_t> out, length_t n) {
        assert(in.size() + out.size() == static_cast<size_t>(n));
        compute_complement(in, out.data(), n);
    }
    static void compute_complement(crindexvec in, rindexvec out, length_t n) {
        assert(in.size() + out.size() == n);
        compute_complement(std::span{in.data(), static_cast<size_t>(in.size())},
                           out.data(), n);
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

        auto sizes    = this->sizes();
        auto *indices = this->indices().data();
        for (index_t t = 0; t < N; ++t) { // time steps
            // Generate the set of inactive indices, J
            index_t num_J = build_Jt(t, indices);
            sizes(t)      = num_J; // save number of inactive indices for later
            std::span J{indices, static_cast<size_t>(num_J)};
            // Generate the complement, the set of active indices, K
            compute_complement(J, indices + num_J, n);
            // Prepare for next time step
            indices += n;
        }
    }

  private:
    static void compute_complement(std::span<const index_t> in, index_t *out,
                                   length_t n) {
        length_t c = 0; // components within time step
        length_t k = 0; // index into the array of active indices
        // iterate over the array with indices 'in'
        for (index_t j : in) { //
            for (; c < j; ++c) // for all indices not in J
                out[k++] = c;  // append the index of this component
            ++c;               // skip indices in J, i.e. c == j
        }
        // add final indices not in J
        for (; c < n; ++c)
            out[k++] = c;
    }
};

} // namespace alpaqa::detail