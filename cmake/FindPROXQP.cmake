# Search for proxsuite via find_package(tinyxml2)
include(CMakeFindDependencyMacro)
find_dependency(proxsuite NO_MODULE)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PROXQP DEFAULT_MSG proxsuite_FOUND)
