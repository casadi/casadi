#include <panoc-alm/interop/cutest/CUTEstLoader.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <yaml-cpp/emitter.h>
#include <yaml-cpp/emittermanip.h>
#include <yaml-cpp/yaml.h>

#include <sstream>

inline YAML::Emitter &operator<<(YAML::Emitter &out, const pa::vec &v) {
    out << YAML::Flow;
    out << YAML::BeginSeq;
    for (pa::vec::Index i = 0; i < v.size(); ++i)
        out << v[i];
    out << YAML::EndSeq;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 CUTEstProblem::Report::Status s) {
    return out << enum_name(s);
}

inline YAML::Emitter &operator<<(YAML::Emitter &out,
                                 const CUTEstProblem::Report &r) {
    out << YAML::BeginMap;
    out << YAML::Key << r.name;
    out << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "nvar" << YAML::Value << r.nvar;
    out << YAML::Key << "ncon" << YAML::Value << r.ncon;
    out << YAML::Key << "status" << YAML::Value << r.status;
    out << YAML::Key << "calls" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "objective" << YAML::Value << r.calls.objective
        << YAML::Key << "objective_grad" << YAML::Value
        << r.calls.objective_grad //
        << YAML::Key << "objective_hess" << YAML::Value
        << r.calls.objective_hess //
        << YAML::Key << "hessian_times_vector" << YAML::Value
        << r.calls.hessian_times_vector;
    if (r.ncon > 0)
        out << YAML::Key << "constraints" << YAML::Value
            << r.calls.constraints //
            << YAML::Key << "constraints_grad" << YAML::Value
            << r.calls.constraints_grad //
            << YAML::Key << "constraints_hess" << YAML::Value
            << r.calls.constraints_hess;
    out << YAML::EndMap;
    out << YAML::Key << "time_setup" << YAML::Value << r.time_setup //
        << YAML::Key << "time" << YAML::Value << r.time;
    out << YAML::EndMap << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, pa::EvalCounter ctr) {
    out << YAML::BeginMap;
    out << YAML::Key << "f" << YAML::Value << ctr.f;
    out << YAML::Key << "grad_f" << YAML::Value << ctr.grad_f;
    out << YAML::Key << "g" << YAML::Value << ctr.g;
    out << YAML::Key << "grad_g" << YAML::Value << ctr.grad_g;
    out << YAML::EndMap;
    return out;
}

inline YAML::Emitter &operator<<(YAML::Emitter &out, pa::SolverStatus s) {
    return out << enum_name(s);
}
