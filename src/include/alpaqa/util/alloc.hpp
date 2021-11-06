#pragma once

#include <alpaqa/util/problem.hpp>

#include <algorithm>
#include <cassert>
#include <stack>
#include <utility>
#include <vector>

namespace alpaqa {

class vec_allocator {
  private:
    size_t num_vec;
    Eigen::Index n;
    std::vector<real_t> storage;
    struct findstack : std::stack<real_t *, std::vector<real_t *>> {
        using stack::stack;

        auto begin() { return this->c.begin(); }
        auto begin() const { return this->c.begin(); }
        auto end() { return this->c.end(); }
        auto end() const { return this->c.end(); }

        auto cbegin() const { return this->c.cbegin(); }
        auto cend() const { return this->c.cend(); }

        auto rbegin() { return this->c.rbegin(); }
        auto rbegin() const { return this->c.rbegin(); }
        auto rend() { return this->c.rend(); }
        auto rend() const { return this->c.rend(); }

        auto crbegin() const { return this->c.crbegin(); }
        auto crend() const { return this->c.crend(); }
    } stack;
    size_t highwatermark = 0;

  public:
    vec_allocator(size_t num_vec, Eigen::Index n)
        : num_vec(num_vec), n(n), storage(num_vec * n, NaN) {
        for (auto it = storage.begin(); it < storage.end(); it += n)
            stack.push(&*it);
    }

    vec_allocator(const vec_allocator &) = delete;
    vec_allocator(vec_allocator &&)      = delete;
    vec_allocator &operator=(const vec_allocator &) = delete;
    vec_allocator &operator=(vec_allocator &&) = delete;

    struct alloc_raii_wrapper {
        using mvec = Eigen::Map<vec>;
        mvec v;
        vec_allocator *alloc;

        alloc_raii_wrapper(real_t *dptr, Eigen::Index n, vec_allocator *alloc)
            : v{dptr, n}, alloc{alloc} {}
        alloc_raii_wrapper(mvec &&v, vec_allocator *alloc)
            : v{std::move(v)}, alloc{alloc} {}
        ~alloc_raii_wrapper() {
            assert(alloc);
            alloc->free(v);
        }
        alloc_raii_wrapper(const alloc_raii_wrapper &) = delete;
        alloc_raii_wrapper(alloc_raii_wrapper &&o)
            : v{std::exchange(o.v, {nullptr, 0})}, //
              alloc{std::exchange(o.alloc, nullptr)} {}
        alloc_raii_wrapper &operator=(const alloc_raii_wrapper &) = delete;
        alloc_raii_wrapper &operator=(alloc_raii_wrapper &&o) {
            this->v     = std::exchange(o.v, {nullptr, 0});
            this->alloc = std::exchange(o.alloc, nullptr);
            return *this;
        }

        alloc_raii_wrapper &operator=(crvec v) {
            this->v = v;
            return *this;
        }
        operator crvec() const { return v; }
        operator rvec() { return v; }
    };

    auto alloc() {
        if (stack.empty())
            throw std::bad_alloc();
        auto dptr = stack.top();
        stack.pop();
        highwatermark = std::max(used_space(), highwatermark);
        return Eigen::Map<vec>(dptr, n);
    }

    alloc_raii_wrapper alloc_raii() { return {alloc(), this}; }

    void free(rvec v) {
        auto dptr = v.data();
        assert(dptr >= &*storage.begin());
        assert(dptr <= &*storage.end() - n);
        assert(std::find(stack.begin(), stack.end(), dptr) == stack.end() &&
               "double free");
        stack.push(dptr);
    }

    template <class... Vecs>
    void free(rvec first, Vecs &&...vecs) {
        free(first);
        free(vecs...);
    }

    size_t size() const { return stack.size(); }
    size_t used_space() const { return num_vec - size(); }

    Eigen::Index vector_size() const { return n; }

    size_t highwater() const { return highwatermark; }
};

} // namespace alpaqa