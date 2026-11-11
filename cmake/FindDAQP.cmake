# DAQP installs daqpConfig.cmake (with a lower-case package name).
find_package(daqp CONFIG QUIET)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DAQP DEFAULT_MSG daqp_FOUND)
