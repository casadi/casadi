#include <alpaqa/implementation/params/params.tpp>
#include <alpaqa/ipopt-adapter-export.h>

#include <IpIpoptApplication.hpp>

#include <stdexcept>

namespace alpaqa::params {

template <class T>
static auto possible_keys(const T &tbl) {
    if (tbl.empty())
        return std::string{};
    auto penult       = std::prev(tbl.end());
    auto quote_concat = [](std::string &&a, auto b) {
        return a + "'" + b.first + "', ";
    };
    return std::accumulate(tbl.begin(), penult, std::string{}, quote_concat) +
           "'" + std::string(penult->first) + "'";
}

template <>
void IPOPT_ADAPTER_EXPORT set_param(Ipopt::IpoptApplication &app,
                                    ParamString s) {

    // Split the key to get the option name (val_key is expected to be empty)
    auto [opt_name, val_key] = split_key(s.key);
    ParamString val_param{
        .full_key = s.full_key,
        .key      = val_key,
        .value    = s.value,
    };

    // Search the option name in the list of Ipopt options
    const auto &ipopt_opts = app.RegOptions()->RegisteredOptionsList();
    const auto regops_it   = ipopt_opts.find(std::string(opt_name));
    if (regops_it == ipopt_opts.end())
        throw std::invalid_argument(
            "Invalid key '" + std::string(opt_name) + "' for type '" +
            "IpoptApplication" + "' in '" + std::string(s.full_key) +
            "',\n  possible keys are: " + possible_keys(ipopt_opts));

    // Depending on the type, set the value of the option
    bool success    = false;
    const auto type = regops_it->second->Type();
    switch (type) {
        case Ipopt::OT_Number: {
            double value;
            set_param(value, val_param);
            success = app.Options()->SetNumericValue(std::string(opt_name),
                                                     value, false);
        } break;
        case Ipopt::OT_Integer: {
            Ipopt::Index value;
            set_param(value, val_param);
            success = app.Options()->SetIntegerValue(std::string(opt_name),
                                                     value, false);
        } break;
        case Ipopt::OT_String: {
            success = app.Options()->SetStringValue(
                std::string(opt_name), std::string(val_param.value), false);
        } break;
        case Ipopt::OT_Unknown:
        default: {
            throw std::invalid_argument("Unknown type in '" +
                                        std::string(s.full_key) + "'");
        }
    }
    if (!success)
        throw std::invalid_argument("Invalid option in '" +
                                    std::string(s.full_key) + "'");
}

} // namespace alpaqa::params
