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

  File::File(const std::string& fname) {
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
  }

  ParsedFile::ParsedFile(const File& file) {
    // Loop over the lines
    auto line_it = file.lines.cbegin();
    while (line_it!=file.lines.cend()) {
      // Current line number and command string
      const string& cmd = *line_it++;
      int line_no = line_it - file.lines.cbegin();

      // If comment or empty line, skip
      if (cmd.empty() || cmd.at(0)=='#') continue;

      // Make sure command
      casadi_assert(cmd.at(0)==':');

      // New entry
      stringstream ss;
      while (true) {
        casadi_assert_message(line_it!=file.lines.cend(),
                              "End of file reached looking for " + cmd);
        const string& line = *line_it++;

        // End of entry?
        if (line==cmd) break;

        // Add to entry
        ss << line << endl;
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

  std::string ParsedFile::to_text(const std::string& cmd) const {
    casadi_assert_message(has(cmd), "No such command: " + cmd);
    return commands.at(cmd).second;
  }

  bool ParsedFile::has(const std::string& cmd) const {
    return commands.find(cmd) != commands.end();
  }

} // namespace casadi
