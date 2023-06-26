#include <alpaqa/export.h>
#include <alpaqa/implementation/util/print.tpp>

namespace alpaqa {

using Eigen::MatrixX;
using Eigen::VectorX;
using std::complex;
using std::ostream;

#ifndef DOXYGEN
// clang-format off

template ALPAQA_EXPORT std::string float_to_str(float value, int precision);
template ALPAQA_EXPORT std::string float_to_str(double value, int precision);
template ALPAQA_EXPORT std::string float_to_str(long double value, int precision);

template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<float>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<double>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<long double>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<float>>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<double>>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<long double>>> &, std::string_view, std::string_view, std::string_view);

template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<float>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<long double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<float>>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<double>>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<long double>>> &, std::string_view);

template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<float>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<long double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<float>>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<double>>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<long double>>> &, std::string_view);

#ifdef ALPAQA_WITH_QUAD_PRECISION
template ALPAQA_EXPORT std::string float_to_str(__float128 value, int precision);

template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<__float128>> &, std::string_view, std::string_view, std::string_view);
template ALPAQA_EXPORT ostream &print_csv_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<__float128>>> &, std::string_view, std::string_view, std::string_view);

template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<__float128>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_matlab_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<__float128>>> &, std::string_view);

template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<__float128>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python_impl(ostream &os, const Eigen::Ref<const MatrixX<complex<__float128>>> &, std::string_view);
#endif

// clang-format on
#endif // DOXYGEN

} // namespace alpaqa