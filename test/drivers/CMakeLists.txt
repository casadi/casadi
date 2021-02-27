find_package(yaml-cpp)
if (yaml-cpp_FOUND)
    add_library(yaml-encoder src/YAMLEncoder.cpp)
    target_include_directories(yaml-encoder
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
        PRIVATE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    )
    target_link_libraries(yaml-encoder PUBLIC yaml-cpp cutest-loader panoc-alm)

    # pa::ALMParams almparams;
    # almparams.max_iter        = 1000;
    # almparams.preconditioning = true;
    # pa::PANOCParams panocparams;
    # panocparams.max_iter = 1000;
    # panocparams.experimental.specialized_lbfgs = true;
    # panocparams.experimental.update_lipschitz_in_linesearch = true;
    set(CUTEst_PROBLEM_LIST
        "ROSENBR"
        "HS111LNP"
        "HS90"
        "HATFLDH"
        "HALDMADS"
        "HAIFAM"
        "GOFFIN"
        "GMNCASE1"
        "GMNCASE2"
        "GMNCASE3"
        "GMNCASE4"
        "GIGOMEZ1"
        "GIGOMEZ2"
        "GIGOMEZ3"
        "FLETCHER"
        "FEEDLOC"
        "FCCU"
        "EXPFITA"
        "EXPFITB"
        "EXPFITC"
        "ERRINBAR"
        "EQC"
        # "ELATTAR"
        "EG2"
        "DUALC1"
        "DUALC2"
        "DUALC5"
        "DUALC8"
        "DUAL1"
        "DUAL2"
        "DUAL3"
        "DUAL4"
        "DNIEPER"
        "DIXCHLNG"
        "DISCS"
        "DISC2"
        "DIPIGRI"
        "DEMYMALO"
        "DEMBO7"
        "DEGENLPA"
        "DEGENLPB"
        "DECONVNE"
        "DECONVU"
        "DECONVC" # slow, many failed L-BFGS steps
        "DECONVBNE"
        "DANWOOD"
        "DANIWOOD"
        "DALLASS"
        "DALLASM"
        # "DALLASL" # PANOC huge error
        "CUBENE"
        "CRESC50" # ???
        "C-RELOAD"
        "CORE2" # PANOC fail
        "COOLHANS" # outer iterations
        "COOLHANSLS"
        "CONGIGMZ"
        "CONCON"
        "COATINGNE" # slow, many failed L-BGFS
        "CLUSTER"
        "CLEUVEN7" # slow (3m45s)
        # "CLEUVEN2" # slow (5m45s) PANOC huge error
        "CHACONN2"
        "CHACONN1"
        "CBS"
        "CANTILVR"
        "BYRDSPHR"
        # "BURKEHAN" # strange ... ε = 0 but δ = 1
        "BT13"
        "BT12"
        "BT11"
        "BT10"
        "BT9"
        "BT8"
        "BT7"
        "BT6"
        "BT5"
        "BT4"
        "BT3"
        "BT2"
        "BT1"
        "BROWNBSNE"
        "BRITGAS"
        # "BRIDGEND" # Looks hard, PANOC makes little progress, "Line search failed" in first (?) iteration
        "BOXBOD"
        "BOX3NE"
        "BOOTH"
        "BIGGSC4"
        "BIGGS6NE"
        "BEALENE"
        "BATCH" # What causes ∇ψₖ to suddenly become inf?
        "BARDNE" # Many yᵀs <= ε failures
        "BA-L1SP"
        "BA-L1"
        "A4X12" # PANOC struggling
        "AVION2"
        "AVGASB"
        "AVGASA"
        "ANTWERP"
        "ALSOTAME" # ε = 0 but δ = 0.6
        "ALLINITC"
        "ALLINITA"
        "AIRPORT"
        "AIRCRFTA"
        # "AGG" # Many yᵀs <= ε failures, tiny steps, little progress
        # "ACOPR300"
        # "ACOPR118" # constraint violation is really small, but penalty is 1e16 and PANOC barely converges. ε = 0 but might be a numerical error because of tiny step size and huge penalty
        # "ACOPR57" # idem
        # "ACOPR30" # idem
        "ACOPR14"
        # "ACOPP300" # same as ACOPR118
        # "ACOPP118" # idem
        # "ACOPP57" # idem
        "ACOPP30" # idem
        "ACOPP14"
        "A5NSSNSM" # large but fast
        "A5NSSSSL" # really large and quite fast for its size (2m45s)
        "A5NSDSDM"
        "A5ESSNDL" # really large and quite fast for its size (2m15s)
        "3PK"

        "OSORIO" # 10201 var, 202 con, pretty fast (1m05s)
        "PENTAGON"
        "PFIT1"
        "PFIT2"
        "PFIT3"
        "PFIT4"
        # "POLAK1" # PANOC didn't converge
        "POLAK2"
        "POLAK3"
        "POLAK4"
        "POLAK5"
        "POLAK6"

        "PORTFL1"
        "PORTFL2"
        "PORTFL3"
        "PORTFL4"
        "PORTFL6"

        "POWELLBS" # small 2×2 problem with inf ∇
        "POWELLSQ"

        "PRICE3NE"
        "PRICE4NE"

        "PRIMAL1"
        "PRIMAL2"
        "PRIMAL3"
        "PRIMAL4"
        "PRIMALC1"
        "PRIMALC2"
        "PRIMALC5"
        "PRIMALC8"

        "PRODPL0"
        "PRODPL1"
        "QC"
        "QCNEW" # many yᵀs <= ε
        "QPCBLEND"
        # "QPCBOEI1" # little progress, why does the panoc error go up even if Σ stays the same?
        # "QPCBOEI2" # idem
        "QPCSTAIR" # ε = 0
        "QPNBLEND"
        # "QPNBOEI1" # little progress
        # "QPNBOEI2" # idem
        # "QPNSTAIR" # many line search failures, why does PANOC error go up?
        "RAT42" # failed
        "RAT43" # many ‖s‖² <= ε failures
        # "READING6" # δ = 1e20
        "READING7" # ψ gets very large and negative very quickly
        "RES"
        "RK23"
        "ROBOT"
        "ROSENMMX" # many ‖s‖² <= ε failures, huge ψ and ∇ψ
        # "ROSZMAN1"
        "ROTDISC" # might get there eventually if enough PANOC iter
        "S268"
        "S365"
        "S365MOD"
        "S368"
        # "SANTA"
        # "SARO" # really big and really slow
        "SMBANK"
        "SMMPSF"
        "SNAKE"
        "SPANHYD"
    )
    string(JOIN "\", \"" CUTEst_PROBLEM_LIST_STRINGS ${CUTEst_PROBLEM_LIST})
    string(CONCAT CUTEst_PROBLEM_LIST_STRINGS "\"" "${CUTEst_PROBLEM_LIST_STRINGS}" "\"")
    configure_file(include/drivers/CUTEstProblemList.h.in
                   include/drivers/CUTEstProblemList.h @ONLY)

    set(TESTNAME "baseline")
    set(TESTRESULT_DIR ${CMAKE_BINARY_DIR}/testresults/${TESTNAME}/CUTEst)
    file(MAKE_DIRECTORY ${TESTRESULT_DIR})
    foreach(P ${CUTEst_PROBLEM_LIST})
        cutest_sif_problem(${P})
    endforeach()

    add_executable(driver src/driver.cpp)
    target_link_libraries(driver PRIVATE yaml-encoder)

    foreach(P ${CUTEst_PROBLEM_LIST})
        add_custom_command(OUTPUT ${TESTRESULT_DIR}/${P}.yaml
                           COMMAND driver ${P} ${TESTRESULT_DIR}
                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                           DEPENDS CUTEst::problem-${P} driver)
        add_custom_target(cutest-result-${P}-${TESTNAME}
                          DEPENDS ${TESTRESULT_DIR}/${P}.yaml)
        list(APPEND CUTEst_PROBLEMS_RESULT_TARGETS 
                    cutest-result-${P}-${TESTNAME})
    endforeach()
    add_custom_target(cutest-results-${TESTNAME}
                      DEPENDS ${CUTEst_PROBLEMS_RESULT_TARGETS})
endif()