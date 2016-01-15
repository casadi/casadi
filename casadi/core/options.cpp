/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#include "options.hpp"
#include <algorithm>

using namespace std;

namespace casadi {

  const Options::Entry* Options::find(const std::string& name) const {
    // Check if in one of the bases
    for (auto&& b : bases) {
      const Options::Entry* entry = b->find(name);
      if (entry) return entry;
    }

    // Lookup in this class
    auto it = entries.find(name);
    if (it!=entries.end()) {
      return &it->second;
    } else {
      return 0;
    }
  }

  void Options::Entry::print(const std::string& name, std::ostream &stream) const {
    stream << "> \"" << name << "\"          ["
           << GenericType::get_type_description(this->type)
           << "] ";

    // Print out the description on a new line.
    stream << "     \"" << this->description << "\""<< std::endl;
  }

  void Options::print(std::ostream &stream) const {
    // Print bases
    for (auto&& b : bases) {
      b->print(stream);
    }

    // Print all entries
    for (auto&& e : entries) {
      e.second.print(e.first, stream);
    }
  }

  double Options::word_distance(const std::string &a, const std::string &b) {
    /// Levenshtein edit distance
    if (a == b) return 0;
    int na = a.size();
    int nb = b.size();
    if (na == 0) return nb;
    if (nb == 0) return na;

    vector<int> v0(nb+1, 0);
    vector<int> v1(nb+1, 0);

    for (int i=0;i<nb+1;++i)
      v0[i] = i;

    char s;
    char t;
    std::locale loc;
    for (int i=0;i<na;i++) {
      v1[0] = i + 1;
      for (int j=0; j<nb; j++) {
        s = std::tolower(a[i], loc);
        t = std::tolower(b[j], loc);
        int cost = 0;
        if (s != t)
          cost = 1;

        v1[j+1] = min(min(v1[j] + 1, v0[j+1] + 1), v0[j] + cost);
      }

      for (int j=0; j<nb+1; j++)
        v0[j] = v1[j];
    }

    return v1[nb];
  }

  vector<string> Options::suggestions(const string& word, int amount) const {
    // Best distances so far
    const double inf = numeric_limits<double>::infinity();
    vector<pair<double, string> > best(amount, {inf, ""});

    // Iterate over elements
    best_matches(word, best);

    // Sort the elements in ascending order
    stable_sort(best.begin(), best.end());

    // Collect the values that are non-infinite
    vector<string> ret;
    ret.reserve(amount);
    for (auto&& e : best) {
      if (e.first!=inf) {
        ret.push_back(e.second);
      }
    }
    return ret;
  }

  void Options::best_matches(const std::string& word,
                             vector<pair<double, string> >& best) const {
    // Iterate over bases
    for (auto&& b : bases) {
      b->best_matches(word, best);
    }

    // Worst match so far
    auto worst = max_element(best.begin(), best.end());

    // Loop over entries
    for (auto&& e : entries) {
      // Get word distance
      double d = word_distance(e.first, word);

      // Keep if better than the worst amongst the suggestions
      if (d < worst->first) {
        worst->first = d;
        worst->second = e.first;
        worst = max_element(best.begin(), best.end());
      }
    }
  }

} // namespace casadi
