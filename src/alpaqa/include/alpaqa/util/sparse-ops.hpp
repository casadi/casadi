#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#include <alpaqa/util/enumerate.hpp>
#include <alpaqa/util/iter-adapter.hpp>
#include <alpaqa/util/set-intersection.hpp>

#include <chrono>
#include <iostream>
#include <ranges>
#include <span>

#include <Eigen/Sparse>

namespace alpaqa::util {

namespace detail {

/// Returns a range over the row indices in the given column of @p sp_mat that
/// are also in @p mask.
/// The range consists of the full Eigen InnerIterators (row, column, value).
template <class SpMat, class MaskVec>
auto select_rows_in_col(const SpMat &sp_mat, MaskVec &mask, auto column) {
    // Make a range that iterates over all matrix elements in the given column:
    using row_iter_t = typename SpMat::InnerIterator;
    util::iter_range_adapter<row_iter_t> col_range{{sp_mat, column}};
    // Projector that extracts the row index from an element of that range:
    static constexpr auto proj_row = [](const row_iter_t &it) {
        return static_cast<typename MaskVec::value_type>(it.row());
    };
    // Compute the intersection between the matrix elements and the mask:
    auto intersection = util::iter_set_intersection(std::move(col_range),
                                                    std::ranges::ref_view{mask},
                                                    std::less{}, proj_row);
    // Extract just the iterator to the matrix element (dropping the mask):
    auto extract_eigen_iter = []<class T>(T &&tup) -> decltype(auto) {
        return std::get<0>(std::forward<T>(tup));
    };
    return std::views::transform(std::move(intersection), extract_eigen_iter);
}

/// Like @ref select_rows_in_col, but returns a range of tuples containing the
/// Eigen InnerIterator and a linear index into the mask.
template <class SpMat, class MaskVec>
auto select_rows_in_col_iota(const SpMat &sp_mat, MaskVec &mask, auto column) {
    // Make a range that iterates over all matrix elements in the given column:
    using row_iter_t = typename SpMat::InnerIterator;
    util::iter_range_adapter<row_iter_t> col_range{{sp_mat, column}};
    // Projector that extracts the row index from an element of that range:
    static constexpr auto proj_row = [](const row_iter_t &it) {
        return static_cast<typename MaskVec::value_type>(it.row());
    };
    // Make a range of tuples of the index into the mask and the mask value:
    auto iota_mask = util::enumerate(std::ranges::ref_view{mask});
    // Projector that extracts the mask value from such a tuple:
    static constexpr auto proj_mask = [](const auto &tup) -> decltype(auto) {
        return std::get<1>(tup);
    };
    // Compute the intersection between the matrix elements and the mask:
    auto intersection =
        util::iter_set_intersection(std::move(col_range), std::move(iota_mask),
                                    std::less{}, proj_row, proj_mask);
    // Extract the iterator to the matrix element and the index into the mask:
    auto extract_eigen_iter_and_index = []<class T>(T && tup)
        requires(std::is_rvalue_reference_v<T &&>)
    {
        auto &[eigen_iter, enum_tup] = tup;
        auto &mask_index             = std::get<0>(enum_tup);
        return std::tuple{std::move(eigen_iter), std::move(mask_index)};
    };
    return std::views::transform(std::move(intersection),
                                 extract_eigen_iter_and_index);
}

} // namespace detail

/// R += R_full(mask,mask)
template <class SpMat, class Mat, class MaskVec>
void sparse_add_masked(const SpMat &R_full, Mat &&R, const MaskVec &mask) {
    // Iterate over all columns in the mask
    for (auto [ci, c] : util::enumerate(mask))
        // Iterate over rows in intersection of mask and sparse column
        for (auto [r, ri] : detail::select_rows_in_col_iota(R_full, mask, c))
            R(ri, ci) += r.value();
}

/// S += S_full(mask,:)
template <class SpMat, class Mat, class MaskVec>
void sparse_add_masked_rows(const SpMat &S_full, Mat &&S, const MaskVec &mask) {
    using index_t = typename SpMat::Index;
    // Iterate over all columns
    for (index_t c = 0; c < S_full.cols(); ++c)
        // Iterate over rows in intersection of mask and sparse column
        for (auto [r, ri] : detail::select_rows_in_col_iota(S_full, mask, c))
            S(ri, c) += r.value();
}

/// out += R(mask_J,mask_K) * v(mask_K);
template <class SpMat, class CVec, class Vec, class MaskVec>
void sparse_matvec_add_masked_rows_cols(const SpMat &R, const CVec &v,
                                        Vec &&out, const MaskVec &mask_J,
                                        const MaskVec &mask_K) {
    using detail::select_rows_in_col_iota;
    // Iterate over all columns in the mask K
    for (auto c : mask_K)
        // Iterate over rows in intersection of mask J and sparse column
        for (auto &&[r, ri] : select_rows_in_col_iota(R, mask_J, c))
            out(ri) += r.value() * v(c);
}

/// out += S(mask,:)ᵀ * v(mask);
template <class SpMat, class CVec, class Vec, class MaskVec>
void sparse_matvec_add_transpose_masked_rows(const SpMat &S, const CVec &v,
                                             Vec &&out, const MaskVec &mask) {
    using index_t = typename SpMat::Index;
    // Iterate over all rows of Sᵀ
    for (index_t c = 0; c < S.cols(); ++c)
        // Iterate over columns in intersection of mask K and sparse row
        for (auto r : detail::select_rows_in_col(S, mask, c))
            out(c) += r.value() * v(r.row());
}

#if __cpp_lib_ranges_zip >= 202110L && __cpp_lib_ranges_enumerate >= 202302L

template <Config Conf>
void convert_triplets_to_ccs(const auto &rows, const auto &cols,
                             rindexvec<Conf> inner_idx,
                             rindexvec<Conf> outer_ptr,
                             index_t<Conf> idx_0 = 0) {
    USING_ALPAQA_CONFIG(Conf);
    // Inner indices: simply the row indices
    assert(std::size(rows) == std::size(inner_idx));
    auto cvt_indices = [&](auto i) { return static_cast<index_t>(i) - idx_0; };
    std::ranges::ref_view rows_vw = rows;
    std::ranges::transform(rows_vw, std::begin(inner_idx), cvt_indices);
    // Outer indices: need to count the number of nonzeros per column
    auto cols_iter = std::begin(cols);
    for (auto &&[i, outer] : std::views::enumerate(outer_ptr)) {
        cols_iter = std::lower_bound(cols_iter, std::end(cols), i + idx_0);
        outer     = static_cast<index_t>(cols_iter - std::begin(cols));
    }
}

/// Sort the (row, column, value) triplets, column first, then row.
template <class... Ts>
void sort_triplets(Ts &&...triplets) {
    // Sort the indices (column first, then row)
    auto cmp = [](const auto &a, const auto &b) {
        return std::tie(std::get<1>(a), std::get<0>(a)) <
               std::tie(std::get<1>(b), std::get<0>(b));
    };
    auto indices = std::views::zip(std::ranges::ref_view{triplets}...);
    auto t0      = std::chrono::steady_clock::now();
    std::ranges::sort(indices, cmp);
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Sorting took: "
              << std::chrono::duration<double>{t1 - t0}.count() * 1e6
              << " µs\n";
}

#endif

} // namespace alpaqa::util