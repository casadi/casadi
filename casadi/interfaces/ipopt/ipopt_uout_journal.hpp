/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_IPOPT_UOUT_JOURNAL_HPP
#define CASADI_IPOPT_UOUT_JOURNAL_HPP

#include "casadi/core/casadi_logger.hpp"
#include <IpJournalist.hpp>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

/// \cond INTERNAL
namespace casadi {

  /** \brief Ipopt journal printing through casadi::uout()

      uout() is thread_local, so the stream is looked up at every write rather than stored:
      a StreamJournal holding &uout() keeps pointing at the stream of the thread that
      created the solver memory, which is destroyed when that thread exits (e.g. a
      map(..., "thread") worker). The next solve on another thread then crashes in
      Journalist::FlushBuffer, even with printing disabled.

      Shared by the ipopt and bonmin plugins (header-only).
  */
  class UoutJournal : public Ipopt::Journal {
  public:
    UoutJournal(const std::string& name, Ipopt::EJournalLevel default_level)
      : Ipopt::Journal(name, default_level) {}

  protected:
    void PrintImpl(Ipopt::EJournalCategory, Ipopt::EJournalLevel, const char* str) override {
      uout() << str;
    }

    void PrintfImpl(Ipopt::EJournalCategory, Ipopt::EJournalLevel,
        const char* pformat, va_list ap) override {
      va_list ap2;
      va_copy(ap2, ap);
      int n = std::vsnprintf(nullptr, 0, pformat, ap2);
      va_end(ap2);
      if (n < 0) return;
      std::vector<char> buf(n + 1);
      std::vsnprintf(buf.data(), buf.size(), pformat, ap);
      uout() << buf.data();
    }

    void FlushBufferImpl() override {
      uout().flush();
    }
  };

} // namespace casadi
/// \endcond

#endif // CASADI_IPOPT_UOUT_JOURNAL_HPP
