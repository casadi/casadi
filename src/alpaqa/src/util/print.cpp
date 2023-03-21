
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

template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<float> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<double> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<long double> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<complex<float>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<complex<double>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<complex<long double>> &);

template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<float> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<double> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<long double> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<complex<float>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<complex<double>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<complex<long double>> &);

template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<float> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<double> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<long double> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<complex<float>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<complex<double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<complex<long double>> &, std::string_view);

template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<float> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<double> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<long double> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<complex<float>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<complex<double>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<complex<long double>> &, std::string_view);

#ifdef ALPAQA_WITH_QUAD_PRECISION
template ALPAQA_EXPORT std::string float_to_str(__float128 value, int precision);

template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<__float128> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<complex<__float128>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<__float128> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<complex<__float128>> &);

template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<__float128> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const MatrixX<complex<__float128>> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<__float128> &, std::string_view);
template ALPAQA_EXPORT ostream &print_python(ostream &os, const VectorX<complex<__float128>> &, std::string_view);
#endif

// clang-format on
#endif // DOXYGEN

} // namespace alpaqa