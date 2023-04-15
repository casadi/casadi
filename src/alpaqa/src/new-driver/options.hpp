#pragma once

#include <alpaqa/params/params.hpp>

#include <memory>
#include <span>
#include <string_view>
#include <vector>

class Options {
  private:
    std::vector<std::string_view> options_storage;
    std::unique_ptr<bool[]> used_storage;

  public:
    Options(int argc, const char *const argv[]) {
        std::copy(argv, argv + argc, std::back_inserter(options_storage));
        used_storage = std::make_unique<bool[]>(options_storage.size());
    }
    [[nodiscard]] std::span<const std::string_view> options() const {
        return options_storage;
    }
    [[nodiscard]] std::span<bool> used() {
        return std::span{used_storage.get(), options_storage.size()};
    }
};

template <class T>
decltype(auto) set_params(T &t, std::string_view prefix, Options &opts) {
    return alpaqa::params::set_params(t, prefix, opts.options(), opts.used());
}
