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


#include "function.hpp"
#include <iomanip>

using namespace casadi;

int eval_dump(const std::string& name) {
    // Load function
    Function f = Function::load(name+".casadi");
    f.change_option("dump_in", false);
    f.change_option("dump_out", false);

    // Helper to format filenames 
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(6);

    std::vector<DM> inputs;
    // Loop over all inputs
    for (int i=0;i<1000000;++i) {
        ss << i;
        try {
            inputs = f.generate_in(name+"."+ss.str()+ ".in.txt");
        } catch (CasadiException& ex) {
            if (i==0) {
                casadi_warning(ex.what());
                casadi_assert(i>0, "Could not find a single input file "
                                 "with file name " + name + ".<dddddd>.in.txt");
            }
            // No more input files
            break;
        }
        // Run function
        std::vector<DM> res = f(inputs);
        // Generate output file
        f.generate_out(name+"."+ss.str()+ ".out.txt", res);
        ss.clear();
    }
    return 0;
}

int eval_dump_parse(const std::vector<std::string>& args) {
    casadi_assert(args.size()>0, "Name is missing in $ casadi-cli eval_dump name.");
    std::string name = args[0];
    return eval_dump(name);
}

int main(int argc, char* argv[]) {
    // Retrieve all arguments
    std::vector<std::string> args(argv + 1, argv + argc);

    // Branch on 'command' (first argument)
    std::set<std::string> commands = {"eval_dump"};
    casadi_assert(args.size()>0, "Must provide a command. Use one of: " + str(commands) + ".");
    std::string cmd = args[0];
    if (cmd=="eval_dump") {
        return eval_dump_parse(std::vector<std::string>(args.begin()+1, args.end()));
    } else {
        casadi_assert(commands.find(cmd)!=commands.end(),
            "Unrecognised command '" + cmd + "'. Use one of: " + str(commands) + ".");
    }
}