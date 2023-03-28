#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace alpaqa {

template <Config Conf     = DefaultConfig,
          class Allocator = std::allocator<std::byte>>
class TestTypeErasedProblem
    : public util::TypeErased<ProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = ProblemVTable<Conf>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using TypeErased::TypeErased;

  public:
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TestTypeErasedProblem make(Args &&...args) {
        return TypeErased::template make<TestTypeErasedProblem, T>(
            std::forward<Args>(args)...);
    }
};

} // namespace alpaqa

struct TestReqProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    TestReqProblem()          = default;
    virtual ~TestReqProblem() = default;
    TestReqProblem(const TestReqProblem &) { throw std::logic_error("copy"); }
    TestReqProblem(TestReqProblem &&) { throw std::logic_error("move"); }

    // clang-format off
    MOCK_METHOD(void, eval_proj_diff_g, (crvec g, rvec e), (const));
    MOCK_METHOD(void, eval_proj_multipliers, (rvec y, real_t M), (const));
    MOCK_METHOD(real_t, eval_prox_grad_step, (real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p), (const));
    MOCK_METHOD(real_t, eval_f, (crvec x), (const));
    MOCK_METHOD(void, eval_grad_f, (crvec x, rvec grad_fx), (const));
    MOCK_METHOD(void, eval_g, (crvec x, rvec gx), (const));
    MOCK_METHOD(void, eval_grad_g_prod, (crvec x, crvec y, rvec grad_gxy), (const));
    MOCK_METHOD(void, check, (), (const));
    // clang-format on

    length_t get_n() const { return 0; }
    length_t get_m() const { return 0; }
};

TEST(TypeErasedProblem, RequiredProblem) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob = alpaqa::TestTypeErasedProblem<>::make<TestReqProblem>();
    vec x;

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_proj_diff_g);
    te_prob.vtable.eval_proj_diff_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_proj_multipliers);
    te_prob.vtable.eval_proj_multipliers(te_prob.self, x, 0);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_prox_grad_step);
    te_prob.vtable.eval_prox_grad_step(te_prob.self, 0, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_f);
    te_prob.vtable.eval_f(te_prob.self, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_grad_f);
    te_prob.vtable.eval_grad_f(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_g);
    te_prob.vtable.eval_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    EXPECT_CALL(te_prob.as<TestReqProblem>(), eval_grad_g_prod);
    te_prob.vtable.eval_grad_g_prod(te_prob.self, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestReqProblem>());

    // No defaults for second-order functions
    EXPECT_EQ(te_prob.vtable.eval_grad_gi, te_prob.vtable.default_eval_grad_gi);
    EXPECT_EQ(te_prob.vtable.eval_hess_L_prod,
              te_prob.vtable.default_eval_hess_L_prod);
    EXPECT_EQ(te_prob.vtable.eval_hess_L, te_prob.vtable.default_eval_hess_L);

    // Defaults for combined evaluations
    EXPECT_EQ(te_prob.vtable.eval_f_grad_f,
              te_prob.vtable.default_eval_f_grad_f);
    EXPECT_EQ(te_prob.vtable.eval_f_g, te_prob.vtable.default_eval_f_g);
    EXPECT_EQ(te_prob.vtable.eval_grad_f_grad_g_prod,
              te_prob.vtable.default_eval_grad_f_grad_g_prod);

    // Defaults for Lagrangians
    EXPECT_EQ(te_prob.vtable.eval_grad_L, te_prob.vtable.default_eval_grad_L);
    EXPECT_EQ(te_prob.vtable.eval_ψ, te_prob.vtable.default_eval_ψ);
    EXPECT_EQ(te_prob.vtable.eval_grad_ψ_from_ŷ,
              te_prob.vtable.default_eval_grad_ψ_from_ŷ);
    EXPECT_EQ(te_prob.vtable.eval_grad_ψ, te_prob.vtable.default_eval_grad_ψ);
    EXPECT_EQ(te_prob.vtable.eval_ψ_grad_ψ,
              te_prob.vtable.default_eval_ψ_grad_ψ);
}

struct TestOptProblem : TestReqProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    TestOptProblem()          = default;
    virtual ~TestOptProblem() = default;
    TestOptProblem(const TestOptProblem &) { throw std::logic_error("copy"); }
    TestOptProblem(TestOptProblem &&) { throw std::logic_error("move"); }

    // clang-format off
    MOCK_METHOD(void, eval_grad_gi, (crvec x, index_t i, rvec grad_gi), (const));
    MOCK_METHOD(void, eval_hess_L_prod, (crvec x, crvec y, real_t scale, crvec v, rvec Hv), (const));
    MOCK_METHOD(void, eval_hess_L, (crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values), (const));
    MOCK_METHOD(real_t, eval_f_grad_f, (crvec x, rvec grad_fx), (const));
    MOCK_METHOD(real_t, eval_f_g, (crvec x, rvec g), (const));
    MOCK_METHOD(void, eval_grad_f_grad_g_prod, (crvec x, crvec y, rvec grad_f, rvec grad_gxy), (const));
    MOCK_METHOD(void, eval_grad_L, (crvec x, crvec y, rvec grad_L, rvec work_n), (const));
    MOCK_METHOD(real_t, eval_ψ, (crvec x, crvec y, crvec Σ, rvec ŷ), (const));
    MOCK_METHOD(void, eval_grad_ψ_from_ŷ, (crvec x, crvec ŷ, rvec grad_ψ, rvec work_n), (const));
    MOCK_METHOD(void, eval_grad_ψ, (crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m), (const));
    MOCK_METHOD(real_t, eval_ψ_grad_ψ, (crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m), (const));
    // clang-format on
};

TEST(TypeErasedProblem, OptionalProblem) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob = alpaqa::TestTypeErasedProblem<>::make<TestOptProblem>();
    vec x;
    indexvec i;

    ASSERT_NE(te_prob.vtable.eval_proj_diff_g, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_proj_diff_g);
    te_prob.vtable.eval_proj_diff_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_proj_multipliers, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_proj_multipliers);
    te_prob.vtable.eval_proj_multipliers(te_prob.self, x, 0);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_prox_grad_step, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_prox_grad_step);
    te_prob.vtable.eval_prox_grad_step(te_prob.self, 0, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_f, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f);
    te_prob.vtable.eval_f(te_prob.self, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_f, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_f);
    te_prob.vtable.eval_grad_f(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_g, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_g);
    te_prob.vtable.eval_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_g_prod, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_g_prod);
    te_prob.vtable.eval_grad_g_prod(te_prob.self, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_gi, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_gi);
    te_prob.vtable.eval_grad_gi(te_prob.self, x, 0, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_hess_L_prod, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_hess_L_prod);
    te_prob.vtable.eval_hess_L_prod(te_prob.self, x, x, 1, x, x,
                                    te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_hess_L, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_hess_L);
    te_prob.vtable.eval_hess_L(te_prob.self, x, x, 1, i, i, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_f_grad_f, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f_grad_f);
    te_prob.vtable.eval_f_grad_f(te_prob.self, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_f_g, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f_g);
    te_prob.vtable.eval_f_g(te_prob.self, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_f_grad_g_prod, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_f_grad_g_prod);
    te_prob.vtable.eval_grad_f_grad_g_prod(te_prob.self, x, x, x, x,
                                           te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_L, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_L);
    te_prob.vtable.eval_grad_L(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_ψ, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_ψ);
    te_prob.vtable.eval_ψ(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_ψ_from_ŷ, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_ψ_from_ŷ);
    te_prob.vtable.eval_grad_ψ_from_ŷ(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_grad_ψ, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_ψ);
    te_prob.vtable.eval_grad_ψ(te_prob.self, x, x, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    ASSERT_NE(te_prob.vtable.eval_ψ_grad_ψ, nullptr);
    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_ψ_grad_ψ);
    te_prob.vtable.eval_ψ_grad_ψ(te_prob.self, x, x, x, x, x, x,
                                 te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());
}

TEST(TypeErasedProblem, CountedOptionalProblem) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    using Problem = alpaqa::ProblemWithCounters<TestOptProblem &>;
    TestOptProblem prob;
    auto te_prob = alpaqa::TestTypeErasedProblem<>::make<Problem>(prob);
    auto &evals  = *te_prob.as<Problem>().evaluations;
    vec x(1);
    indexvec i;

    EXPECT_EQ(evals.proj_diff_g, 0);
    ASSERT_NE(te_prob.vtable.eval_proj_diff_g, nullptr);
    EXPECT_CALL(prob, eval_proj_diff_g);
    te_prob.vtable.eval_proj_diff_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.proj_diff_g, 1);

    EXPECT_EQ(evals.proj_multipliers, 0);
    ASSERT_NE(te_prob.vtable.eval_proj_multipliers, nullptr);
    EXPECT_CALL(prob, eval_proj_multipliers);
    te_prob.vtable.eval_proj_multipliers(te_prob.self, x, 0);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.proj_multipliers, 1);

    EXPECT_EQ(evals.prox_grad_step, 0);
    ASSERT_NE(te_prob.vtable.eval_prox_grad_step, nullptr);
    EXPECT_CALL(prob, eval_prox_grad_step);
    te_prob.vtable.eval_prox_grad_step(te_prob.self, 0, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.prox_grad_step, 1);

    EXPECT_EQ(evals.f, 0);
    ASSERT_NE(te_prob.vtable.eval_f, nullptr);
    EXPECT_CALL(prob, eval_f);
    te_prob.vtable.eval_f(te_prob.self, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.f, 1);

    EXPECT_EQ(evals.grad_f, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_f, nullptr);
    EXPECT_CALL(prob, eval_grad_f);
    te_prob.vtable.eval_grad_f(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_f, 1);

    EXPECT_EQ(evals.g, 0);
    ASSERT_NE(te_prob.vtable.eval_g, nullptr);
    EXPECT_CALL(prob, eval_g);
    te_prob.vtable.eval_g(te_prob.self, x, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.g, 1);

    EXPECT_EQ(evals.grad_g_prod, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_g_prod, nullptr);
    EXPECT_CALL(prob, eval_grad_g_prod);
    te_prob.vtable.eval_grad_g_prod(te_prob.self, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_g_prod, 1);

    EXPECT_EQ(evals.grad_gi, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_gi, nullptr);
    EXPECT_CALL(prob, eval_grad_gi);
    te_prob.vtable.eval_grad_gi(te_prob.self, x, 0, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_gi, 1);

    EXPECT_EQ(evals.hess_L_prod, 0);
    ASSERT_NE(te_prob.vtable.eval_hess_L_prod, nullptr);
    EXPECT_CALL(prob, eval_hess_L_prod);
    te_prob.vtable.eval_hess_L_prod(te_prob.self, x, x, 1, x, x,
                                    te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.hess_L_prod, 1);

    EXPECT_EQ(evals.hess_L, 0);
    ASSERT_NE(te_prob.vtable.eval_hess_L, nullptr);
    EXPECT_CALL(prob, eval_hess_L);
    te_prob.vtable.eval_hess_L(te_prob.self, x, x, 1, i, i, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.hess_L, 1);

    EXPECT_EQ(evals.f_grad_f, 0);
    ASSERT_NE(te_prob.vtable.eval_f_grad_f, nullptr);
    EXPECT_CALL(prob, eval_f_grad_f);
    te_prob.vtable.eval_f_grad_f(te_prob.self, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.f_grad_f, 1);

    EXPECT_EQ(evals.f_g, 0);
    ASSERT_NE(te_prob.vtable.eval_f_g, nullptr);
    EXPECT_CALL(prob, eval_f_g);
    te_prob.vtable.eval_f_g(te_prob.self, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.f_g, 1);

    EXPECT_EQ(evals.grad_f_grad_g_prod, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_f_grad_g_prod, nullptr);
    EXPECT_CALL(prob, eval_grad_f_grad_g_prod);
    te_prob.vtable.eval_grad_f_grad_g_prod(te_prob.self, x, x, x, x,
                                           te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_f_grad_g_prod, 1);

    EXPECT_EQ(evals.grad_L, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_L, nullptr);
    EXPECT_CALL(prob, eval_grad_L);
    te_prob.vtable.eval_grad_L(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_L, 1);

    EXPECT_EQ(evals.ψ, 0);
    ASSERT_NE(te_prob.vtable.eval_ψ, nullptr);
    EXPECT_CALL(prob, eval_ψ);
    te_prob.vtable.eval_ψ(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.ψ, 1);

    EXPECT_EQ(evals.grad_ψ_from_ŷ, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_ψ_from_ŷ, nullptr);
    EXPECT_CALL(prob, eval_grad_ψ_from_ŷ);
    te_prob.vtable.eval_grad_ψ_from_ŷ(te_prob.self, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_ψ_from_ŷ, 1);

    EXPECT_EQ(evals.grad_ψ, 0);
    ASSERT_NE(te_prob.vtable.eval_grad_ψ, nullptr);
    EXPECT_CALL(prob, eval_grad_ψ);
    te_prob.vtable.eval_grad_ψ(te_prob.self, x, x, x, x, x, x, te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.grad_ψ, 1);

    EXPECT_EQ(evals.ψ_grad_ψ, 0);
    ASSERT_NE(te_prob.vtable.eval_ψ_grad_ψ, nullptr);
    EXPECT_CALL(prob, eval_ψ_grad_ψ);
    te_prob.vtable.eval_ψ_grad_ψ(te_prob.self, x, x, x, x, x, x,
                                 te_prob.vtable);
    testing::Mock::VerifyAndClearExpectations(&prob);
    EXPECT_EQ(evals.ψ_grad_ψ, 1);
}

struct TestOptProblemNoHess : TestOptProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    TestOptProblemNoHess()                             = default;
    virtual ~TestOptProblemNoHess()                    = default;
    TestOptProblemNoHess(const TestOptProblemNoHess &) = default;
    TestOptProblemNoHess(TestOptProblemNoHess &&)      = default;

    bool provides_eval_grad_gi() { return true; }
    bool provides_eval_hess_L_prod() { return false; }
    bool provides_eval_hess_L() { return false; }
};

TEST(TypeErasedProblem, providesNoHess) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob =
        alpaqa::TestTypeErasedProblem<>::make<TestOptProblemNoHess>();

    EXPECT_NE(te_prob.vtable.eval_grad_gi, te_prob.vtable.default_eval_grad_gi);
    EXPECT_EQ(te_prob.vtable.eval_hess_L_prod,
              te_prob.vtable.default_eval_hess_L_prod);
    EXPECT_EQ(te_prob.vtable.eval_hess_L, te_prob.vtable.default_eval_hess_L);
}

struct TestOptProblemNoPsi : TestOptProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    TestOptProblemNoPsi()                            = default;
    virtual ~TestOptProblemNoPsi()                   = default;
    TestOptProblemNoPsi(const TestOptProblemNoPsi &) = default;
    TestOptProblemNoPsi(TestOptProblemNoPsi &&)      = default;

    bool provides_eval_grad_L() { return true; }
    bool provides_eval_ψ() { return false; }
    bool provides_eval_ψ_grad_ψ() { return false; }
};

TEST(TypeErasedProblem, providesNoPsi) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob = alpaqa::TestTypeErasedProblem<>::make<TestOptProblemNoPsi>();

    EXPECT_NE(te_prob.vtable.eval_grad_L, te_prob.vtable.default_eval_grad_L);
    EXPECT_EQ(te_prob.vtable.eval_ψ, te_prob.vtable.default_eval_ψ);
    EXPECT_EQ(te_prob.vtable.eval_ψ_grad_ψ,
              te_prob.vtable.default_eval_ψ_grad_ψ);
}

TEST(TypeErasedProblem, TEOptionalProblem) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob = alpaqa::TypeErasedProblem<>::make<TestOptProblem>();
    vec x;
    indexvec i;

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_proj_diff_g);
    te_prob.eval_proj_diff_g(x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_proj_multipliers);
    te_prob.eval_proj_multipliers(x, 0);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_prox_grad_step);
    te_prob.eval_prox_grad_step(0, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f);
    te_prob.eval_f(x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_f);
    te_prob.eval_grad_f(x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_g);
    te_prob.eval_g(x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_g_prod);
    te_prob.eval_grad_g_prod(x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_gi);
    te_prob.eval_grad_gi(x, 0, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_hess_L_prod);
    te_prob.eval_hess_L_prod(x, x, 1, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_hess_L);
    te_prob.eval_hess_L(x, x, 1, i, i, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f_grad_f);
    te_prob.eval_f_grad_f(x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_f_g);
    te_prob.eval_f_g(x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_f_grad_g_prod);
    te_prob.eval_grad_f_grad_g_prod(x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_L);
    te_prob.eval_grad_L(x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_ψ);
    te_prob.eval_ψ(x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_ψ_from_ŷ);
    te_prob.eval_grad_ψ_from_ŷ(x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_grad_ψ);
    te_prob.eval_grad_ψ(x, x, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());

    EXPECT_CALL(te_prob.as<TestOptProblem>(), eval_ψ_grad_ψ);
    te_prob.eval_ψ_grad_ψ(x, x, x, x, x, x);
    testing::Mock::VerifyAndClearExpectations(&te_prob.as<TestOptProblem>());
}

TEST(TypeErasedProblem, TEprovidesNoHess) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    auto te_prob = alpaqa::TypeErasedProblem<>::make<TestOptProblemNoHess>();
    vec x;
    indexvec i;

    EXPECT_TRUE(te_prob.provides_eval_grad_gi());
    EXPECT_FALSE(te_prob.provides_eval_hess_L_prod());
    EXPECT_FALSE(te_prob.provides_eval_hess_L());

    EXPECT_CALL(te_prob.as<TestOptProblemNoHess>(), eval_grad_gi);
    te_prob.eval_grad_gi(x, 0, x);
    testing::Mock::VerifyAndClearExpectations(
        &te_prob.as<TestOptProblemNoHess>());

    EXPECT_THROW(te_prob.eval_hess_L_prod(x, x, 1, x, x),
                 alpaqa::not_implemented_error);

    EXPECT_THROW(te_prob.eval_hess_L(x, x, 1, i, i, x),
                 alpaqa::not_implemented_error);
}
