#include <alpaqa/config/config.hpp>

#include <alpaqa/problem/src/problem.tpp>

namespace alpaqa {

ALPAQA_EXPORT_TEMPLATE(class, ProblemBase, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, ProblemBase, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, ProblemBase, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, ProblemBase, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, ProblemBase, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(class, Problem, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, Problem, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, Problem, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, Problem, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, Problem, EigenConfigq);
#endif

ALPAQA_EXPORT_TEMPLATE(class, FunctionalProblem, DefaultConfig);
ALPAQA_EXPORT_TEMPLATE(class, FunctionalProblem, EigenConfigf);
ALPAQA_EXPORT_TEMPLATE(class, FunctionalProblem, EigenConfigd);
ALPAQA_EXPORT_TEMPLATE(class, FunctionalProblem, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_TEMPLATE(class, FunctionalProblem, EigenConfigq);
#endif

} // namespace alpaqa