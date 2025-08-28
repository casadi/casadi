[conf]
tools.cmake.cmake_layout:build_folder_vars=['const.matlab', 'settings.build_type']

[options]
alpaqa/*:with_ipopt=False
alpaqa/*:with_external_casadi=True
alpaqa/*:with_qpalm=False
alpaqa/*:with_cutest=False
alpaqa/*:with_json=True
alpaqa/*:with_matlab=True
alpaqa/*:with_drivers=False

[settings]
casadi/*:build_type=Release
