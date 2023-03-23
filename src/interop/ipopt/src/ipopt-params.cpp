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
void IPOPT_ADAPTER_EXPORT set_params(Ipopt::IpoptApplication &app,
                                     std::string_view prefix,
                                     std::span<const std::string_view> options,
                                     std::optional<std::span<bool>> used) {
    const auto ipopt_opts = app.RegOptions()->RegisteredOptionsList();

    size_t index = 0;
    for (const auto &kv : options) {
        auto [key, value]     = split_key(kv, '=');
        auto [pfx, remainder] = split_key(key);
        auto curr_index       = index++;
        if (pfx != prefix)
            continue;
        if (used)
            (*used)[curr_index] = true;

        auto [opt_key, val_key] = split_key(remainder);
        ParamString val_param{.full_key = kv, .key = val_key, .value = value};

        const auto regops_it = ipopt_opts.find(std::string(opt_key));
        if (regops_it == ipopt_opts.end())
            throw std::invalid_argument(
                "Invalid key '" + std::string(key) + "' for type '" +
                "IpoptApplication" + "' in '" + std::string(kv) +
                "',\n  possible keys are: " + possible_keys(ipopt_opts));

        bool success    = false;
        const auto type = regops_it->second->Type();
        switch (type) {
            case Ipopt::OT_Number: {
                double value;
                set_param(value, val_param);
                success = app.Options()->SetNumericValue(std::string(opt_key),
                                                         value, false);
            } break;
            case Ipopt::OT_Integer: {
                Ipopt::Index value;
                set_param(value, val_param);
                success = app.Options()->SetIntegerValue(std::string(opt_key),
                                                         value, false);
            } break;
            case Ipopt::OT_String: {
                success = app.Options()->SetStringValue(
                    std::string(opt_key), std::string(val_param.value), false);
            } break;
            case Ipopt::OT_Unknown:
            default: {
                throw std::invalid_argument("Unknown type in '" +
                                            std::string(kv) + "'");
            }
        }
        if (!success)
            throw std::invalid_argument("Invalid option in '" +
                                        std::string(kv) + "'");
    }
}

} // namespace alpaqa::params
