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


      #include "shell_compiler.hpp"
      #include <string>

      const std::string casadi::ShellCompiler::meta_doc=
      "\n"
"\n"
"\n"
"Interface to the JIT compiler SHELL\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_22w \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+----------------------+-----------------+---------------------------------+\n"
"|          Id          |      Type       |           Description           |\n"
"+======================+=================+=================================+\n"
"| cleanup              | OT_BOOL         | Cleanup temporary files when    |\n"
"|                      |                 | unloading. Default: true        |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| compiler             | OT_STRING       | Compiler command                |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| compiler_flags       | OT_STRINGVECTOR | Alias for 'compiler_flags'      |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| compiler_output_flag | OT_STRING       | Compiler flag to denote object  |\n"
"|                      |                 | output. Default: '-o '          |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| compiler_setup       | OT_STRING       | Compiler setup command.         |\n"
"|                      |                 | Intended to be fixed. The       |\n"
"|                      |                 | 'flag' option is the prefered   |\n"
"|                      |                 | way to set custom flags.        |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| directory            | OT_STRING       | Directory to put temporary      |\n"
"|                      |                 | objects in. Must end with a     |\n"
"|                      |                 | file separator.                 |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| extra_suffixes       | OT_STRINGVECTOR | List of suffixes for extra      |\n"
"|                      |                 | files that the compiler may     |\n"
"|                      |                 | generate. Default: None         |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| flags                | OT_STRINGVECTOR | Compile flags for the JIT       |\n"
"|                      |                 | compiler. Default: None         |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| linker               | OT_STRING       | Linker command                  |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| linker_flags         | OT_STRINGVECTOR | Linker flags for the JIT        |\n"
"|                      |                 | compiler. Default: None         |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| linker_output_flag   | OT_STRING       | Linker flag to denote shared    |\n"
"|                      |                 | library output. Default: '-o '  |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| linker_setup         | OT_STRING       | Linker setup command. Intended  |\n"
"|                      |                 | to be fixed. The 'flag' option  |\n"
"|                      |                 | is the prefered way to set      |\n"
"|                      |                 | custom flags.                   |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| name                 | OT_STRING       | The file name used to write out |\n"
"|                      |                 | compiled objects/libraries. The |\n"
"|                      |                 | actual file names used depend   |\n"
"|                      |                 | on 'temp_suffix' and include    |\n"
"|                      |                 | extensions. Default:            |\n"
"|                      |                 | 'tmp_casadi_compiler_shell'     |\n"
"+----------------------+-----------------+---------------------------------+\n"
"| temp_suffix          | OT_BOOL         | Use a temporary (seemingly      |\n"
"|                      |                 | random) filename suffix for     |\n"
"|                      |                 | file names. This is desired for |\n"
"|                      |                 | thread-safety. This behaviour   |\n"
"|                      |                 | may defeat caching compiler     |\n"
"|                      |                 | wrappers. Default: true         |\n"
"+----------------------+-----------------+---------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
