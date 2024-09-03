# Search for fatrop via find_package(fatrop)
include(CMakeFindDependencyMacro)
find_dependency(fatrop NO_MODULE)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FATROP DEFAULT_MSG fatrop_FOUND)
