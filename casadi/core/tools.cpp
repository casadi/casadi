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


#include "tools.hpp"

#include "casadi_misc.hpp"
#include "importer.hpp"
#include "serializer.hpp"
#include <casadi/config.h>
#include "casadi_os.hpp"
#include "casadi_meta.hpp"

namespace casadi {

void callback_stdout(const char* s) {
    uout() << s << std::endl;
}

void callback_stderr(const char* s) {
    uerr() << s << std::endl;
}

Function external_transform(const std::string& name,
                    const std::string& op,
                    const Function& f,
                    const Dict& opts) {
    std::string signature = "f";
    Importer li(name + SHARED_LIBRARY_SUFFIX, "dll");
    std::string op_full = op + "__" + signature;
    external_transform_t t = (external_transform_t) li.get_function(op_full);
    if (!t) {
        casadi_error("Transformation '" + op + "' not defined in library '" + name +
                     "' for the following signature: " + signature);
    }

    StringSerializer ss;
    ss.pack(f);
    ss.pack(opts);

    char api_version = 0;
    std::string casadi_version = CasadiMeta::version();

    std::string s = ss.encode();
    const char* out = t(api_version, casadi_version.c_str(), s.c_str(),
        callback_stdout, callback_stderr);
    if (!out) {
        casadi_error("Transformation '" + op + "' defined in library '" + name +
                     "' failed.");
    }

    StringDeserializer sd(out);
    Function r = sd.unpack_function();

    return r;
}

} // namespace casadi

const char* external_transform_test_success__f(char api_version, const char* casadi_version,
        const char* in,
        casadi::external_print_callback_t stdout, casadi::external_print_callback_t stderr) {
    if (api_version != 0) {
        stderr("version mismatch");
        return 0;
    }
    casadi::StringDeserializer sd(in);
    casadi::Function f = sd.unpack_function();
    casadi::Dict opts = sd.unpack_generictype();

    std::string msg = "passed options: " + str(opts);

    stdout(msg.c_str());

    stdout("Doing a lot of stuff...");
    stderr("Warning here...");

    casadi::StringSerializer ss;
    ss.pack(f);
    static std::string s = ss.encode();

    return s.c_str();
}

const char* external_transform_test_fail__f(char api_version, const char* casadi_version,
        const char* in,
        casadi::external_print_callback_t stdout, casadi::external_print_callback_t stderr) {
    stdout("This is going to fail");
    stderr("Fatal error");
    return 0;
}
