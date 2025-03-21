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


#include "libzip.hpp"
#include "casadi/core/filesystem_impl.hpp"
#include <zip.h>

namespace casadi {

    zip_t* open_zip_from_istream(std::istream& stream) {
        // Read stream content into a string
        std::string buffer((std::istreambuf_iterator<char>(stream)),
            std::istreambuf_iterator<char>());

        // Open zip archive from memory buffer
        zip_error_t errorp;
        zip_source_t* src = zip_source_buffer_create(buffer.data(), buffer.size(), 0, &errorp);
        if (!src) {
            casadi_error("Failed to create zip source: " +
                std::string(zip_error_strerror(&errorp)) + "\n");
            return nullptr;
        }

        zip_t* archive = zip_open_from_source(src, 0, &errorp);
        if (!archive) {
            zip_source_free(src);
            casadi_error("Failed to open zip from source: " +
                std::string(zip_error_strerror(&errorp)) + "\n");
        }
        return archive;
    }

    bool extract_zip_internal(zip_t* za, const std::string& output_dir); // Declare before use

    bool extract_zip_from_stream(std::istream& src, const std::string& output_dir) {
        return extract_zip_internal(open_zip_from_istream(src), output_dir);
    }

    bool extract_zip_from_path(const std::string& zip_path, const std::string& output_dir) {
        int err;
        zip_t* za = zip_open(zip_path.c_str(), ZIP_RDONLY, &err);
        if (!za) {
            casadi_error("Cannot open ZIP file: " + zip_path);
            return false;
        }
        return extract_zip_internal(za, output_dir);
    }

    bool extract_zip_internal(zip_t* za, const std::string& output_dir) {
        Filesystem::assert_enabled();
        auto filesystem = Filesystem::getPlugin("ghc");

        zip_int64_t num_entries = zip_get_num_entries(za, 0);
        if (num_entries < 0) {
            zip_close(za);
            casadi_error("Cannot read ZIP contents.");
            return false;
        }

        for (zip_uint64_t i = 0; i < static_cast<zip_uint64_t>(num_entries); ++i) {
            const char* name = zip_get_name(za, i, 0);
            if (!name) {
                uerr() << "Error: Cannot get file name for entry " << i << std::endl;
                continue;
            }

            std::string full_path = output_dir + filesep() + name;
            std::string filesep_str = filesep();
            char filesep_char = filesep_str[0];

            if (full_path.back() == filesep_char) {  // Directory entry
                filesystem.exposed.create_directories(full_path);
            } else {  // File entry
                std::string dir_path = full_path.substr(0, full_path.find_last_of(filesep_char));
                filesystem.exposed.create_directories(dir_path);

                zip_file_t* zf = zip_fopen_index(za, i, 0);
                if (!zf) {
                    uerr() << "Error: Cannot open file in ZIP: " << name << std::endl;
                    continue;
                }

                std::ofstream out_file(full_path, std::ios::binary);
                if (!out_file) {
                    uerr() << "Error: Cannot write file: " << full_path << std::endl;
                    zip_fclose(zf);
                    continue;
                }

                char buffer[8192];
                zip_int64_t bytes_read;
                while ((bytes_read = zip_fread(zf, buffer, sizeof(buffer))) > 0) {
                    out_file.write(buffer, bytes_read);
                }

                if (bytes_read < 0) {
                    uerr() << "Error: Read failed for file in ZIP: " << name << std::endl;
                }

                out_file.close();
                zip_fclose(zf);
            }
        }

        zip_close(za);
        return true;
    }

   extern "C"
   int CASADI_ARCHIVER_LIBZIP_EXPORT
   casadi_register_archiver_libzip(Archiver::Plugin* plugin) {
     plugin->name = "libzip";
     plugin->doc = Libzip::meta_doc.c_str();
     plugin->version = CASADI_VERSION;
     plugin->exposed.unpack = &extract_zip_from_path;
     plugin->exposed.unpack_from_stream = &extract_zip_from_stream;
     return 0;
   }

   extern "C"
   void CASADI_ARCHIVER_LIBZIP_EXPORT casadi_load_archiver_libzip() {
     Archiver::registerPlugin(casadi_register_archiver_libzip);
   }


} // namespace casadi
