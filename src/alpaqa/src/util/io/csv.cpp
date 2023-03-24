#include <alpaqa/export.h>
#include <alpaqa/implementation/util/io/csv.tpp>

namespace alpaqa::csv {

template void ALPAQA_EXPORT //
read_row_impl(std::istream &, Eigen::Ref<Eigen::VectorX<float>>, char);
template void ALPAQA_EXPORT //
read_row_impl(std::istream &, Eigen::Ref<Eigen::VectorX<double>>, char);
template void ALPAQA_EXPORT //
read_row_impl(std::istream &, Eigen::Ref<Eigen::VectorX<long double>>, char);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void ALPAQA_EXPORT //
read_row_impl(std::istream &, Eigen::Ref<Eigen::VectorX<__float128>>, char);
#endif

template std::vector<float> ALPAQA_EXPORT //
read_row_std_vector(std::istream &, char);
template std::vector<double> ALPAQA_EXPORT //
read_row_std_vector(std::istream &, char);
template std::vector<long double> ALPAQA_EXPORT //
read_row_std_vector(std::istream &, char);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template std::vector<__float128> ALPAQA_EXPORT //
read_row_std_vector(std::istream &, char);
#endif

} // namespace alpaqa::csv