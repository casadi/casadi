#pragma once

#include <alpaqa/config/config.hpp>
#include <Eigen/Sparse>
#include <memory>
#include <span>
#include <type_traits>
#include <variant>

namespace alpaqa {

template <Config Conf>
struct DenseMatrixStructure {
    USING_ALPAQA_CONFIG(Conf);
    DenseMatrixStructure(length_t rows, length_t cols) : rows{rows}, cols{cols} {}
    [[nodiscard]] auto view(std::span<const real_t> storage) const {
        assert(static_cast<length_t>(storage.size()) >= rows * cols);
        return cmmat{storage.data(), rows, cols};
    }
    length_t rows, cols;
    static constexpr bool is_symmetric = false;
};

template <Config Conf>
struct DiagonalMatrixStructure {
    USING_ALPAQA_CONFIG(Conf);
    DiagonalMatrixStructure(length_t rows, [[maybe_unused]] length_t cols) : rows{rows} {
        assert(rows == cols);
    }
    [[nodiscard]] auto view(std::span<const real_t> storage) const {
        assert(static_cast<length_t>(storage.size()) >= rows);
        return cmvec{storage.data(), rows}.asDiagonal();
    }
    length_t rows;
    static constexpr bool is_symmetric = true;
};

template <Config Conf, Eigen::UpLoType UpLo>
struct SymmetricMatrixStructure {
    USING_ALPAQA_CONFIG(Conf);
    /// @todo   Wastes ~50% of storage, but can be optimized later.
    SymmetricMatrixStructure(length_t rows, [[maybe_unused]] length_t cols) : rows{rows} {
        assert(rows == cols);
    }
    [[nodiscard]] auto view(std::span<const real_t> storage) const {
        return cmmat{storage.data(), rows, rows}.template selfadjointView<UpLo>();
    }
    length_t rows;
    static constexpr bool is_symmetric = true;
};

template <Config Conf, class StorageIndex>
struct SparseMatrixStructure {
    USING_ALPAQA_CONFIG(Conf);
    struct Sparsity {
        length_t rows, cols, nnz;
        const StorageIndex *outer_index_ptr, *inner_index_ptr, *inner_non_zeros_ptr;
    };
    /// Shared pointer allows owning, reference-counted sparsity patterns,
    /// through `shared_ptr`'s aliasing constructor.
    SparseMatrixStructure(std::shared_ptr<Sparsity> sparsity) : sparsity{std::move(sparsity)} {}
    [[nodiscard]] auto view(std::span<const real_t> storage) const {
        assert(static_cast<length_t>(storage.size()) >= sparsity->nnz);
        return Eigen::Map<const Eigen::SparseMatrix<real_t, Eigen::ColMajor, StorageIndex>>{
            sparsity->rows,
            sparsity->cols,
            sparsity->nnz,
            sparsity->outer_index_ptr,
            sparsity->inner_index_ptr,
            storage.data(),
            sparsity->inner_non_zeros_ptr,
        };
    }
    std::shared_ptr<Sparsity> sparsity;
    static constexpr bool is_symmetric = false;
};

template <class Base>
struct TransposedVariantMatrix;

template <Config Conf>
struct VariantMatrix {
    USING_ALPAQA_CONFIG(Conf);

    // clang-format off
    using variant_t = std::variant<
        DenseMatrixStructure<config_t>,
        DiagonalMatrixStructure<config_t>,
        SymmetricMatrixStructure<config_t, Eigen::Upper>,
        SymmetricMatrixStructure<config_t, Eigen::Lower>,
        SparseMatrixStructure<config_t, int32_t>,
        SparseMatrixStructure<config_t, int64_t>
    >;
    // clang-format on
    variant_t value;

    template <class Visitor>
    decltype(auto) visit(Visitor &&vis) {
        return std::visit(std::forward<Visitor>(vis), value);
    }
    template <class Visitor>
    decltype(auto) visit(Visitor &&vis) const {
        return std::visit(std::forward<Visitor>(vis), value);
    }

    TransposedVariantMatrix<VariantMatrix &> transpose() & { return {*this}; }
    TransposedVariantMatrix<const VariantMatrix &> transpose() const & { return {*this}; }
    TransposedVariantMatrix<VariantMatrix> transpose() && { return {std::move(*this)}; }
};

template <class Base>
struct TransposedVariantMatrix {
    Base value;

    template <class Visitor>
    decltype(auto) visit(Visitor &&vis) {
        return value.visit(std::forward<Visitor>(vis));
    }
    template <class Visitor>
    decltype(auto) visit(Visitor &&vis) const {
        return value.visit(std::forward<Visitor>(vis));
    }
};

enum class VariantMatrixMultiplySide { Left, Right };
enum class VariantMatrixMultiplyTranspose { No, Yes };

template <Config Conf, VariantMatrixMultiplySide Side, VariantMatrixMultiplyTranspose Transp>
struct VariantMatrixMultiplyAdd {
    USING_ALPAQA_CONFIG(Conf);

    template <class T>
    void operator()(const T &op_a) {
        using enum VariantMatrixMultiplySide;
        using enum VariantMatrixMultiplyTranspose;
        if constexpr (Side == Left && (Transp == No || T::is_symmetric))
            out.noalias() += op_a.view(storage) * op_b;
        else if constexpr (Side == Left && Transp == Yes)
            out.noalias() += op_a.view(storage).transpose() * op_b;
        else if constexpr (Side == Right && (Transp == No || T::is_symmetric))
            out.noalias() += op_b * op_a.view(storage);
        else if constexpr (Side == Right && Transp == Yes)
            out.noalias() += op_b * op_a.view(storage).transpose();
    }

    std::span<const real_t> storage;
    crmat op_b;
    rmat out;
};

template <Config Conf>
void multiply_add(std::span<const typename Conf::real_t> storage, const VariantMatrix<Conf> &left,
                  typename Conf::crmat right, typename Conf::rmat out) {
    left.visit(VariantMatrixMultiplyAdd<Conf, VariantMatrixMultiplySide::Left,
                                        VariantMatrixMultiplyTranspose::No>{storage, right, out});
}

template <Config Conf>
void multiply_add(std::span<const typename Conf::real_t> storage, typename Conf::crmat left,
                  const VariantMatrix<Conf> &right, typename Conf::rmat out) {
    right.visit(VariantMatrixMultiplyAdd<Conf, VariantMatrixMultiplySide::Right,
                                         VariantMatrixMultiplyTranspose::No>{storage, left, out});
}

template <class MV>
void multiply_add(std::span<const typename std::remove_cvref_t<MV>::config_t::real_t> storage,
                  const TransposedVariantMatrix<MV> &left,
                  typename std::remove_cvref_t<MV>::config_t::crmat right,
                  typename std::remove_cvref_t<MV>::config_t::rmat out) {
    left.visit(VariantMatrixMultiplyAdd<typename std::remove_cvref_t<MV>::config_t,
                                        VariantMatrixMultiplySide::Left,
                                        VariantMatrixMultiplyTranspose::Yes>{storage, right, out});
}

template <class MV>
void multiply_add(std::span<const typename std::remove_cvref_t<MV>::config_t::real_t> storage,
                  typename std::remove_cvref_t<MV>::config_t::crmat left,
                  const TransposedVariantMatrix<MV> &right,
                  typename std::remove_cvref_t<MV>::config_t::rmat out) {
    right.visit(VariantMatrixMultiplyAdd<typename std::remove_cvref_t<MV>::config_t,
                                         VariantMatrixMultiplySide::Right,
                                         VariantMatrixMultiplyTranspose::Yes>{storage, left, out});
}

} // namespace alpaqa
