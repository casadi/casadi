#include <Eigen/Jacobi>
#include <cstddef>
#include <alpaqa/util/ringbuffer.hpp>
#include <alpaqa/util/vec.hpp>
#include <type_traits>

namespace alpaqa {

/// Incremental QR factorization using modified Gram-Schmidt with
/// reorthogonalization.
///
/// Computes A = QR while allowing efficient removal of the first
/// column of A or adding new columns at the end of A.
class LimitedMemoryQR {
  public:
    LimitedMemoryQR() = default;

    /// @param  n
    ///         The size of the vectors, the number of rows of A.
    /// @param  m
    ///         The maximum number of columns of A.
    ///
    /// The maximum dimensions of Q are n×m and the maximum dimensions of R are
    /// m×m.
    LimitedMemoryQR(size_t n, size_t m) : Q(n, m), R(m, m) {}

    size_t n() const { return Q.rows(); }
    size_t m() const { return Q.cols(); }

    /// Add the given column to the right.
    template <class VecV>
    void add_column(const VecV &v) {
        assert(q_idx < m());

        auto q = Q.col(q_idx);
        auto r = R.col(r_idx_end);

        // Modified Gram-Schmidt to make q orthogonal to Q
        q = v;
        for (size_t i = 0; i < q_idx; ++i) {
            real_t s = Q.col(i).dot(q);
            r(i)     = s;
            q -= s * Q.col(i);
        }

        // Compute the norms of orthogonalized q and original v
        real_t norm_q = q.norm();
        real_t norm_v = v.norm();

        // If ‖q‖ is significantly smaller than ‖v‖, perform
        // reorthogonalization
        real_t η = 0.7;
        while (norm_q < η * norm_v) {
            ++reorth_count;
            for (size_t i = 0; i < q_idx; ++i) {
                real_t s = Q.col(i).dot(q);
                r(i) += s;
                q -= s * Q.col(i);
            }
            norm_v = norm_q;
            norm_q = q.norm();
        }

        // Normalize q such that new matrix (Q q) remains orthogonal (i.e. has
        // orthonormal columns)
        r(q_idx) = norm_q;
        q /= norm_q;

        // Increment indices, add a column to Q and R.
        ++q_idx;
        r_idx_end = r_succ(r_idx_end);
    }

    /// Remove the leftmost column.
    void remove_column() {
        assert(q_idx > 0);

        // After removing the first colomn of the upper triangular matrix R,
        // it becomes upper Hessenberg. Givens rotations are used to make it
        // triangular again.
        Eigen::JacobiRotation<real_t> G;
        size_t r = 0;                   // row index of R
        size_t c = r_succ(r_idx_start); // column index of R in storage
        while (r < q_idx - 1) {
            // Compute the Givens rotation that makes the subdiagonal element
            // of column c or R zero.
            G.makeGivens(R(r, c), R(r + 1, c), &R(r, c));
            // Apply it to the remaining columns of R.
            // Not the current column, because the diagonal element was updated
            // by the makeGivens function, and the subdiagonal element doesn't
            // have to be set to zero explicitly, it's implicit.
            // Also not the previous columns, because they are already
            // implicitly zero below the diagonal and this rotation wouldn't
            // have any effect there.
            // TODO: can this be sped up by applying it in two blocks instead
            //       of column by column?
            for (size_t cc = r_succ(c); cc != r_idx_end; cc = r_succ(cc))
                R.col(cc).applyOnTheLeft(r, r + 1, G.adjoint());
            // Apply the inverse of the Givens rotation to Q.
            Q.block(0, 0, Q.rows(), q_idx).applyOnTheRight(r, r + 1, G);
            // Advance indices to next diagonal element of R.
            ++r;
            c = r_succ(c);
        }
        // Remove rightmost column of Q, since it corresponds to the bottom row
        // of R, which was set to zero by the Givens rotations
        --q_idx;
        // Remove the first column of R.
        r_idx_start = r_succ(r_idx_start);
    }

    /// Solve the least squares problem Ax = b.
    template <class VecB, class VecX>
    void solve_col(const VecB &b, VecX &x) const {
        // Iterate over the diagonal of R, starting at the bottom right,
        // this is standard back substitution
        // (recall that R is stored in a circular buffer, so R.col(i) is
        // not the mathematical i-th column)
        auto rev_bgn = ring_reverse_iter().begin();
        auto rev_end = ring_reverse_iter().end();
        auto fwd_end = ring_iter().end();
        for (auto it_d = rev_bgn; it_d != rev_end; ++it_d) {
            // Row/column index of diagonal element of R
            auto [rR, cR] = *it_d;
            // (r is the zero-based mathematical index, c is the index in
            // the circular buffer)
            x(rR) = Q.col(rR).transpose() * b; // Compute rhs Qᵀb
            // In the current row of R, iterate over the elements to the
            // right of the diagonal
            // Iterating from left to right seems to give better results
            for (auto it_c = it_d.forwardit; it_c != fwd_end; ++it_c) {
                auto [rX2, cR2] = *it_c;
                x(rR) -= R(rR, cR2) * x(rX2);
            }
            x(rR) /= R(rR, cR); // Divide by diagonal element
        }
    }

    /// Solve the least squares problem AX = B.
    template <class MatB, class MatX>
    void solve(const MatB &B, MatX &X) const {
        assert(B.cols() <= X.cols());
        assert(B.rows() >= Q.rows());
        assert(X.rows() >= Eigen::Index(num_columns()));
        // Each column of the right hand side is solved as an individual system
        for (Eigen::Index cB = 0; cB < B.cols(); ++cB) {
            auto b = B.col(cB);
            auto x = X.col(cB);
            solve_col(b, x);
        }
    }

    template <class Derived>
    using solve_ret_t = std::conditional_t<
        Eigen::internal::traits<Derived>::ColsAtCompileTime == 1, vec, mat>;

    /// Solve the least squares problem AX = B.
    template <class Derived>
    solve_ret_t<Derived> solve(const Eigen::DenseBase<Derived> &B) {
        solve_ret_t<Derived> X(m(), B.cols());
        solve(B, X);
        return X;
    }

    /// Get the full, raw storage for the orthogonal matrix Q.
    const mat &get_raw_Q() const { return Q; }
    /// Get the full, raw storage for the upper triangular matrix R.
    /// The columns of this matrix are permuted because it's stored as a
    /// circular buffer for efficiently appending columns to the end and
    /// popping columns from the front.
    const mat &get_raw_R() const { return R; }

    /// Get the full storage for the upper triangular matrix R but with the
    /// columns in the correct order.
    /// @note   Meant for tests only, creates a permuted copy.
    mat get_full_R() const {
        if (r_idx_start == 0)
            return R;
        // Using a permutation matrix here isn't as efficient as rotating the
        // matrix manually, but this function is only used in tests, so it
        // shouldn't matter.
        Eigen::PermutationMatrix<Eigen::Dynamic> P(R.cols());
        P.setIdentity();
        std::rotate(P.indices().data(), P.indices().data() + r_idx_start,
                    P.indices().data() + P.size());
        return R * P;
    }
    /// Get the matrix R such that Q times R is the original matrix.
    /// @note   Meant for tests only, creates a permuted copy.
    mat get_R() const {
        return get_full_R()
            .block(0, 0, q_idx, q_idx)
            .triangularView<Eigen::Upper>();
    }
    /// Get the matrix Q such that Q times R is the original matrix.
    /// @note   Meant for tests only, creates a copy.
    mat get_Q() const { return Q.block(0, 0, n(), q_idx); }

    /// Multiply the matrix R by a scalar.
    void scale_R(real_t scal) {
        for (auto [i, r_idx] : ring_iter())
            R.col(r_idx).topRows(i + 1) *= scal;
    }

    /// Get the number of MGS reorthogonalizations.
    size_t get_reorth_count() const { return reorth_count; }
    /// Reset the number of MGS reorthogonalizations.
    void clear_reorth_count() { reorth_count = 0; }

    /// Reset all indices, clearing the Q and R matrices.
    void reset() {
        q_idx        = 0;
        r_idx_start  = 0;
        r_idx_end    = 0;
        reorth_count = 0;
    }

    /// Re-allocate storage for a problem with a different size. Causes
    /// a @ref reset.
    void resize(size_t n, size_t m) {
        Q.resize(n, m);
        R.resize(m, m);
        reset();
    }

    /// Get the number of columns that are currently stored.
    size_t num_columns() const { return q_idx; }
    /// Get the head index of the circular buffer (points to the oldest
    /// element).
    size_t ring_head() const { return r_idx_start; }
    /// Get the tail index of the circular buffer (points to one past the most
    /// recent element).
    size_t ring_tail() const { return r_idx_end; }
    /// Get the next index in the circular buffer.
    size_t ring_next(size_t i) const { return r_succ(i); }
    /// Get the previous index in the circular buffer.
    size_t ring_prev(size_t i) const { return r_pred(i); }

    /// Get iterators in the circular buffer.
    CircularRange<size_t> ring_iter() const {
        return {q_idx, r_idx_start, r_idx_end, m()};
    }
    /// Get reverse iterators in the circular buffer.
    ReverseCircularRange<size_t> ring_reverse_iter() const {
        return ring_iter();
    }

  private:
    mat Q; ///< Storage for orthogonal factor Q.
    mat R; ///< Storage for upper triangular factor R.

    size_t q_idx       = 0; ///< Number of columns of Q being stored.
    size_t r_idx_start = 0; ///< Index of the first column of R.
    size_t r_idx_end   = 0; ///< Index of the one-past-last column of R.

    size_t reorth_count = 0; ///< Number of MGS reorthogonalizations.

    /// Get the next index in the circular storage for R.
    size_t r_succ(size_t i) const { return i + 1 < m() ? i + 1 : 0; }
    /// Get the previous index in the circular storage for R.
    size_t r_pred(size_t i) const { return i == 0 ? m() - 1 : i - 1; }
};

} // namespace alpaqa