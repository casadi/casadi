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


#include "casadi_file.hpp"
#include "exception.hpp"

using namespace std;

namespace casadi {

  ParsedFile::ParsedFile(const std::string& fname) {
    parse(fname);
  }

  ParsedFile::ParsedFile(const std::vector<std::string>& lines, int offset) {
    parse(lines, offset);
  }

  void ParsedFile::parse(const std::string& fname) {
    // Lines to be extracted
    std::vector<std::string> lines;

    // Open file and read line-by-line
    ifstream file(fname);
    casadi_assert_message(file.good(), "File \"" + fname + "\" cannot be opened.");
    string line;
    while (true) {
      // Read line
      getline(file, line);
      if (file.eof()) break;
      casadi_assert(file.good());

      // Save to class
      lines.push_back(line);
    }

    // Parse the extracted lines
    parse(lines, 0);
  }

  void ParsedFile::parse(const std::vector<std::string>& lines, int offset) {
    // Loop over the lines
    auto line_it = lines.cbegin();
    while (line_it!=lines.cend()) {
      // If comment or empty line, skip
      if (line_it->empty() || line_it->at(0)=='#') {
        line_it++;
        continue;
      }

      // Current line number
      int line_no = line_it - lines.cbegin() + 1 + offset;

      // Get command string
      casadi_assert_message(line_it->at(0)==':',
                            "Syntax error: " + *line_it + " is not a command string");
      string cmd = line_it->substr(1, line_it->find(' ')-1);

      // New entry
      stringstream ss;

      // Collect the meta data
      size_t start = cmd.size()+2;
      while (true) {
        // Find the backslash, if any
        size_t stop = line_it->find('\\');

        // Add to entry
        ss << line_it->substr(start, stop-start);

        // Break if no more lines or not multiline
        if (++line_it==lines.cend() || stop == string::npos) break;

        // Multiline entry
        if (start!=stop) ss << std::endl;

        // No offset for subsequent lines
        start = 0;
      }

      // Insert new element in map
      auto new_el = commands.insert(make_pair(cmd, make_pair(line_no, ss.str())));
      casadi_assert_message(new_el.second, "Duplicate entry: \"" + cmd + "\"");
    }
  }

  void ParsedFile::print(std::ostream &stream) const {
    // Print all the commands
    for (auto&& c : commands) {
      stream << c.first << " (line " << c.second.first << "):" << endl;
      stream << c.second.second;
    }
  }

  std::string ParsedFile::to_text(const std::string& cmd, int ind) const {
    if (ind>=0) return to_text(indexed(cmd, ind));
    casadi_assert_message(has(cmd), "No such command: " + cmd);
    return commands.at(cmd).second;
  }

  bool ParsedFile::has(const std::string& cmd, int ind) const {
    if (ind>=0) return has(indexed(cmd, ind));
    return commands.find(cmd) != commands.end();
  }

} // namespace casadi
