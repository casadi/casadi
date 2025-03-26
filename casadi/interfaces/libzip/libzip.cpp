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
#include <cstring>
namespace casadi {

    bool extract_zip_internal(zip_t* za, const std::string& output_dir); // Declare before use

    bool extract_zip_from_stringstream(std::stringstream& src, const std::string& output_dir) {
        src.clear();
        src.seekg(0, std::ios::beg);
        const std::string& s = src.str();
        // Open zip archive from memory buffer
        zip_error_t errorp;
        zip_source_t* zip_src = zip_source_buffer_create(s.data(), s.size(), 0, &errorp);

        if (!zip_src) {
            casadi_error("Failed to create zip source: " +
                std::string(zip_error_strerror(&errorp)) + "\n");
            return false;
        }

        zip_t* archive = zip_open_from_source(zip_src, 0, &errorp);
        if (!archive) {
            zip_source_free(zip_src);
            casadi_error("Failed to open zip from source: " +
                std::string(zip_error_strerror(&errorp)) + "\n");
        }
        return extract_zip_internal(archive, output_dir);
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

            std::string full_path = output_dir + "/" + name;
            if (full_path.back() == '/') {  // Directory entry
                filesystem.exposed.create_directories(full_path);
            } else {  // File entry
                std::string dir_path = full_path.substr(0, full_path.find_last_of('/'));
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

    bool add_file_to_zip(zip_t* archive, const std::string& file_path,
            const std::string& archive_name) {
        std::ifstream file(file_path, std::ios::binary | std::ios::ate);
        if (!file) {
            uerr() << "Error: Cannot open file: " << file_path << std::endl;
            return false;
        }

        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);

        char* data = static_cast<char*>(malloc(size));  // ← use malloc
        if (!data) {
            uerr() << "Error: Memory allocation failed for file: " << file_path << std::endl;
            return false;
        }

        if (!file.read(data, size)) {
            uerr() << "Error: Cannot read file: " << file_path << std::endl;
            free(data);
            return false;
        }

        zip_error_t ziperr;
        zip_source_t* source = zip_source_buffer_create(data, size, 1, &ziperr); // ← freep = 1
        if (!source) {
            uerr() << "Error: Cannot create zip source for file: " << file_path
                   << ": " << zip_error_strerror(&ziperr) << std::endl;
            free(data); // not strictly needed, but safe
            zip_error_fini(&ziperr);
            return false;
        }

        zip_int64_t idx = zip_file_add(archive, archive_name.c_str(), source, ZIP_FL_ENC_UTF_8);
        if (idx < 0) {
            zip_source_free(source);  // Only needed if not added
            uerr() << "Error: Cannot add file to archive: " << archive_name << std::endl;
            return false;
        }

        return true;
    }

    void add_directory_recursive(zip_t* archive,
        const std::string& base_dir,
        const std::string& current_dir,
        const std::string& rel_prefix) {
        auto filesystem = Filesystem::getPlugin("ghc");
        std::vector<std::string> entries = filesystem.exposed.iterate_directory_names(current_dir);

        for (const std::string& full_path : entries) {
            std::string rel_path = full_path.substr(base_dir.size() + 1);

            if (filesystem.exposed.is_directory(full_path)) {
                zip_dir_add(archive, (rel_path + "/").c_str(), ZIP_FL_ENC_UTF_8);
                add_directory_recursive(archive, base_dir, full_path, rel_path);
            } else {
                add_file_to_zip(archive, full_path, rel_path);
            }
        }
    }
    bool zip_to_stream(const std::string& dir, std::ostream& output) {
        zip_error_t error;
        zip_error_init(&error);

        zip_source_t* src = zip_source_buffer_create(nullptr, 0, 0, &error);
        if (!src) {
            uerr() << "Failed to create zip source buffer: "
                   << zip_error_strerror(&error) << std::endl;
            zip_error_fini(&error);
            return false;
        }

        // Prevent zip_close from destroying the source
        zip_source_keep(src);

        zip_t* archive = zip_open_from_source(src, ZIP_TRUNCATE, &error);
        if (!archive) {
            uerr() << "Failed to open zip archive from source: "
                   << zip_error_strerror(&error) << std::endl;
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        try {
            add_directory_recursive(archive, dir, dir, "");
        } catch (const std::exception& e) {
            uerr() << "Exception while zipping directory: " << e.what() << std::endl;
            zip_discard(archive);  // also frees src
            zip_error_fini(&error);
            return false;
        }

        if (zip_close(archive) != 0) {
            uerr() << "Failed to finalize zip archive: "
                   << zip_error_strerror(&error) << std::endl;
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        // At this point, src contains the archive in memory.
        if (zip_source_open(src) < 0) {
            uerr() << "Failed to open zip source for reading." << std::endl;
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        // Seek to end to get size
        if (zip_source_seek(src, 0, SEEK_END) < 0) {
            uerr() << "Failed to seek to end of zip source." << std::endl;
            zip_source_close(src);
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        zip_int64_t size = zip_source_tell(src);
        if (size < 0) {
            uerr() << "Failed to get size of zip source." << std::endl;
            zip_source_close(src);
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        if (zip_source_seek(src, 0, SEEK_SET) < 0) {
            uerr() << "Failed to rewind zip source." << std::endl;
            zip_source_close(src);
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        if (zip_source_seek(src, 0, SEEK_SET) < 0) {
            uerr() << "Failed to rewind zip source." << std::endl;
            zip_source_close(src);
            zip_source_free(src);
            zip_error_fini(&error);
            return false;
        }

        // Efficient streaming read/write
        char buf[8192];
        zip_int64_t bytes_read;

        while ((bytes_read = zip_source_read(src, buf, sizeof(buf))) > 0) {
            output.write(buf, bytes_read);
            if (!output) {
                uerr() << "Write error while streaming zip data to output." << std::endl;
                zip_source_close(src);
                zip_source_free(src);
                zip_error_fini(&error);
                return false;
            }
        }

        zip_source_close(src);
        zip_source_free(src);
        zip_error_fini(&error);

        if (bytes_read < 0) {
            uerr() << "Error reading from zip source." << std::endl;
            return false;
        }

        return true;
    }


    bool zip_to_path(const std::string& dir_path, const std::string& zip_path) {
        std::ofstream ofs(zip_path, std::ios::binary);
        if (!ofs) {
            uerr() << "Failed to open output file: " << zip_path << std::endl;
            return false;
        }

        return zip_to_stream(dir_path, ofs);
    }

    bool zip_to_path2(const std::string& dir_path, const std::string& zip_path) {
        int errorp;
        zip_t* archive = zip_open(zip_path.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &errorp);
        if (!archive) {
            zip_error_t ziperror;
            zip_error_init_with_code(&ziperror, errorp);
            uerr() << "Error: Cannot open zip archive " << zip_path << ": "
                   << zip_error_strerror(&ziperror) << std::endl;
            zip_error_fini(&ziperror);
            return false;
        }

        try {
            add_directory_recursive(archive, dir_path, dir_path, "");
        } catch (const std::exception& e) {
            uerr() << "Exception while zipping directory: " << e.what() << std::endl;
            zip_discard(archive);
            return false;
        }

        if (zip_close(archive) < 0) {
            uerr() << "Error: Cannot finalize zip archive: " << zip_strerror(archive) << std::endl;
            zip_discard(archive);
            return false;
        }

        return true;
    }

   extern "C"
   int CASADI_ARCHIVER_LIBZIP_EXPORT
   casadi_register_archiver_libzip(Archiver::Plugin* plugin) {
     plugin->name = "libzip";
     plugin->doc = Libzip::meta_doc.c_str();
     plugin->version = CASADI_VERSION;
     plugin->exposed.unpack = &extract_zip_from_path;
     plugin->exposed.unpack_from_stringstream = &extract_zip_from_stringstream;
     plugin->exposed.pack = &zip_to_path;
     plugin->exposed.pack_to_stream = &zip_to_stream;
     return 0;
   }

   extern "C"
   void CASADI_ARCHIVER_LIBZIP_EXPORT casadi_load_archiver_libzip() {
     Archiver::registerPlugin(casadi_register_archiver_libzip);
   }


} // namespace casadi
