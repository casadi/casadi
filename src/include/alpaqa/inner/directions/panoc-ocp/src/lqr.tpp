#pragma once

#include <alpaqa/inner/directions/panoc-ocp/lqr.hpp>

namespace alpaqa {

template <Config Conf>
void LQRFactor<Conf>::assign_possibly_diagonal(rmat dest, crmat src) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest = src.reshaped().asDiagonal();
    } else {
        dest = src;
    }
}

template <Config Conf>
void LQRFactor<Conf>::add_possibly_diagonal(rmat dest, crmat src) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest += src.reshaped().asDiagonal();
    } else {
        dest += src;
    }
}

template <Config Conf>
void LQRFactor<Conf>::add_possibly_diagonal_masked(rmat dest, crmat src,
                                                   const auto &V) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest += src.reshaped()(V).asDiagonal();
    } else {
        dest += src(V, V);
    }
}

template <Config Conf>
void LQRFactor<Conf>::mat_vec_possibly_diagonal(rvec dest, crmat M, crvec v) {
    if (M.cols() == 1 || M.rows() == 1) {
        dest = M.reshaped().asDiagonal() * v;
    } else {
        dest = M * v;
    }
}

// --- TODO: DRY

template <Config Conf>
void StatefulLQRFactor<Conf>::assign_possibly_diagonal(rmat dest, crmat src) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest = src.reshaped().asDiagonal();
    } else {
        dest = src;
    }
}

template <Config Conf>
void StatefulLQRFactor<Conf>::add_possibly_diagonal(rmat dest, crmat src) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest += src.reshaped().asDiagonal();
    } else {
        dest += src;
    }
}

template <Config Conf>
void StatefulLQRFactor<Conf>::add_possibly_diagonal_masked(rmat dest, crmat src,
                                                           const auto &V) {
    if (src.cols() == 1 || src.rows() == 1) {
        dest += src.reshaped()(V).asDiagonal();
    } else {
        dest += src(V, V);
    }
}

template <Config Conf>
void StatefulLQRFactor<Conf>::mat_vec_possibly_diagonal(rvec dest, crmat M,
                                                        crvec v) {
    if (M.cols() == 1 || M.rows() == 1) {
        dest = M.reshaped().asDiagonal() * v;
    } else {
        dest = M * v;
    }
}

} // namespace alpaqa