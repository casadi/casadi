#include <alpaqa/export.h>
#include <alpaqa/util/src/print.tpp>

namespace alpaqa {

using Eigen::MatrixX;
using Eigen::VectorX;
using std::complex;
using std::ostream;

#ifndef DOXYGEN
// clang-format off

template ALPAQA_EXPORT std::string float_to_str(float value);
template ALPAQA_EXPORT std::string float_to_str(double value);
template ALPAQA_EXPORT std::string float_to_str(long double value);

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

#ifdef ALPAQA_WITH_QUAD_PRECISION
template ALPAQA_EXPORT std::string float_to_str(__float128 value);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<__float128> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const MatrixX<complex<__float128>> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<__float128> &);
template ALPAQA_EXPORT ostream &print_matlab(ostream &os, const VectorX<complex<__float128>> &);
#endif

// clang-format on
#endif // DOXYGEN

} // namespace alpaqa