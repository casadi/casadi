#pragma once

#include <memory>
#include <panoc-alm/problem.hpp>

class CUTEstProblem {
  private:
    CUTEstProblem();

  public:
    CUTEstProblem(CUTEstProblem &&);
    CUTEstProblem &operator=(CUTEstProblem &&);
    ~CUTEstProblem();

  public:
    pa::Problem problem;
    pa::vec x0;
    pa::vec y0;

  private:
    std::unique_ptr<class CUTEstLoader> implementation;

  public:
    friend CUTEstProblem load_CUTEst(const char *so_fname,
                                     const char *outsdif_fname);
};

CUTEstProblem load_CUTEst(const char *so_fname, const char *outsdif_fname);
