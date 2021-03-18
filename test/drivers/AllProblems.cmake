set(CUTEst_PROBLEM_LIST
    "ROSENBR"
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

    "EXPLIN" # Very weird, no progress at all (unconstrained)
    "EXPLIN2"
    "EXPQUAD"
    "EXTROSNBNE" # Faster than Ipopt
    "FERRISDC" # Much faster than Ipopt
    # "FLOSP2HH" # Blows up
    "FLOSP2TM" # Blows up
    "FREURONE" # Blows up
    # "GASOIL" # Slow
    "GAUSSELM"
    "GENROSENE" # Violation remains high, blows up
    "GILBERT" # Faster than Ipopt
    # "GLIDER" # Blows up

    "GOULDQP1" # A bit faster than Ipopt
    "GOULDQP2" # A bit faster than Ipopt
    "GOULDQP3" # A bit faster than Ipopt

    "GRIDGENA" # Inf
    "GRIDNETA" # PANOC stalls
    "GRIDNETH"
    "GRIDNETI" # PANOC stalls
    "GROUPING"
    "GROWTH" # Blows up
    "HADAMARD" # Blows up
    # "HAGER1" # Blows up
    # "HAGER2" # Blows up
    # "HAGER3" # Blows up
    # "HAGER4" # Blows up
    "HAIFAL"
    "HAIFAS"
    "HAHN1" # Blows up
    "HANGING"
    "HARKERP2"
    "HATFLDDNE"
    "HATFLDENE" # Blows up
    "HATFLDF"
    "HATFLDFLNE"
    "HATFLDG"
    "HEART6"
    "HEART8"
    "HELIXNE"
    "HELSBY"
    "HET-Z"
    "HIE1327D" # PANOC stalls
    "HIE1372D" # PANOC stalls

    "HILBERTA"
    "HILBERTB"
    "HIMMELBA"
    "HIMMELBB"
    "HIMMELBC"
    "HIMMELBCLS"
    "HIMMELBD"
    "HIMMELBE"
    "HIMMELBF"
    "HIMMELBFNE"
    "HIMMELBG"
    "HIMMELBH"
    "HIMMELBI"
    "HIMMELBJ"
    "HIMMELBK"
    "HIMMELP1"
    "HIMMELP2"
    "HIMMELP3"
    "HIMMELP4"
    "HIMMELP5"
    "HIMMELP6"
    "HOLMES"
    "HONG"

    "HS10"
    "HS100"
    "HS100LNP"
    "HS100MOD"
    "HS101"
    "HS102"
    "HS103"
    "HS104"
    "HS105"
    "HS106"
    "HS107"
    "HS108"
    "HS109"
    "HS11"
    "HS110"
    "HS111"
    "HS111LNP"
    "HS112"
    "HS113"
    "HS114"
    "HS116"
    "HS117"
    "HS118"
    "HS119"
    "HS12"
    "HS13"
    "HS14"
    "HS15"
    "HS16"
    "HS17"
    "HS18"
    "HS19"
    "HS2"
    "HS2NE"
    "HS20"
    "HS21"
    "HS21MOD"
    "HS22"
    "HS23"
    "HS24"
    "HS25"
    "HS25NE"
    "HS26"
    "HS268"
    "HS27"
    "HS28"
    "HS29"
    "HS3"
    "HS30"
    "HS31"
    "HS32"
    "HS33"
    "HS34"
    "HS35"
    "HS35I"
    "HS35MOD"
    "HS36"
    "HS37"
    "HS38"
    "HS39"
    "HS3MOD"
    "HS4"
    "HS40"
    "HS41"
    "HS42"
    "HS43"
    "HS44"
    "HS44NEW"
    "HS45"
    "HS46"
    "HS47"
    "HS48"
    "HS49"
    "HS5"
    "HS50"
    "HS51"
    "HS52"
    "HS53"
    "HS54"
    "HS55"
    "HS56"
    "HS57"
    "HS59"
    "HS6"
    "HS60"
    "HS61"
    "HS62"
    "HS63"
    "HS64"
    "HS65"
    "HS66"
    "HS67"
    "HS68"
    "HS69"
    "HS7"
    "HS70"
    "HS71"
    "HS72"
    "HS73"
    "HS74"
    "HS75"
    "HS76"
    "HS76I"
    "HS77"
    "HS78"
    "HS79"
    "HS8"
    "HS80"
    "HS81"
    "HS83"
    "HS84"
    "HS85"
    "HS86"
    "HS87"
    "HS88"
    "HS89"
    "HS9"
    "HS90"
    "HS91"
    "HS92"
    "HS93"
    "HS95"
    "HS96"
    "HS97"
    "HS98"
    "HS99"
    "HS99EXP" # Blows up
    "HUBFIT"
    "HUESTIS" # No progress
    "HUES-MOD" # Slower than Ipopt
    "HVYCRASH" # No progress
    "HYDCAR20" # Blows up
    "HYDCAR6" # Blows up
    "HYDROELL" # About as fast as Ipopt
    "HYDROELM" # Faster than Ipopt
    "HYDROELS" # Faster than Ipopt
    "HYPCIR"
    "INTEGREQ" # Faster than Ipopt
    "INTEQNE"
    "JANNSON3" # PANOC stalls
    "JANNSON4" # PANOC stalls
    "JENSMPNE" # Maxiter
    "JJTABEL3"
    "JUDGENE" # Blows up
    "JUNKTURN" # PANOC stalls
    "KIRBY2" # Blows up
    "KISSING" # Blows up
    "KISSING2" 
    "KIWCRESC"
    "KOWOSBNE" # Blows up
    "KSIP"
    "KTMODEL" # NaN
    "LAKES" # NaN
    "LANCZOS1" # Blows up
    "LANCZOS2"
    "LANCZOS3" # Blows up
    "LAUNCH" # Blows up
    "LCH" # Huge error (1e23)
    "LEAKNET" # PANOC stalls, blows up
    "LEUVEN7" # PANOC stalls, blows up
    "LEWISPOL" # PANOC stalls GOOD TEST PROBLEM (6×9)
    "LIARWHDNE"
    "LIN"
    "LINCONT" # Blows up
    "LINSPANH"
    "LISWET12" # PANOC stalls, then blows up
    "LISWET9" # Much slower than Ipopt
    "LOADBAL"
    # "LOBSTERZ" # Huge and slow
    "LOOTSMA" # Blows up GOOD TEST PROBLEM (3×2)
    "LOTSCHD"
    "LSC1" # Blows up
    "LSC2" # Blows up
    "LSNNODOC"
    "LSQFIT"

    "LUKVLE18"

    "LUKVLI1"
    "LUKVLI2" # Huge error (1e46) ?!
    "LUKVLI3" # About as fast as Ipopt
    "LUKVLI4" # Inf (Ipopt: restoration failed)
    "LUKVLI5" # Slower than Ipopt (10s vs 2s)
    "LUKVLI6" # Much slower than Ipopt
    "LUKVLI7" # Much slower than Ipopt
    "LUKVLI8" # Inf, Ipopt is fast
    "LUKVLI9" # PANOC stalls, Ipopt is much faster
    "LUKVLI10" # MaxTime
    "LUKVLI11" # Slower than Ipopt
    "LUKVLI12" # Faster than Ipopt
    "LUKVLI13" # Much slower than Ipopt
    "LUKVLI14" # PANOC stalls
    "LUKVLI15" # Slow
    "LUKVLI16" # Twice as slow as Ipopt
    "LUKVLI17" # Slower than Ipopt
    "LUKVLI18" # Twice as slow as Ipopt
    "MADSEN"
    "MADSSCHJ" # PANOC doesn't converge at all (ε = 1) until suddenly it does, very strange
    "MAKELA1"
    "MAKELA2"
    "MAKELA3"
    "MAKELA4"
    "MANCINONE"
    "MANNE" # Big but fast, much faster than Ipopt (0.65s vs >20 minutes, I didn't wait for Ipopt result, convergence was very close, but Ipopt was converging to same value for f(x))
    "MARINE" # Huge, slow, PANOC doesn't converge, blows up
    "MATRIX2"
    "MCONCON" # Faster than Ipopt
    "MESH" # Blows up
    "METHANL8" # Blows up
    # "METHANOL" # Huge, slow progress, Ipopt is much faster
    "MEYER3NE" # Inf grad
    "MGH09LS"
    "MGH09" # Blows up
    "MGH10LS"
    "MGH10" # Inf grad
    "MGH10SLS"
    "MGH10S" # Blows up (3 var, 16 eq constr)
    "MGH17LS"
    "MGH17" # Inf grad
    "MGH17SLS"
    "MGH17S" # Inf grad
    "MIFFLIN1"
    "MIFFLIN2"
    "MINC44"
    "MINMAXBD"
    "MINMAXRB"
    "MINPERM" # Much faster than Ipopt
    "MINSURFO" # Faster than Ipopt
    "MISRA1A"
    "MISRA1B"
    "MISRA1C" # NaN gradient (Ipopt can't solve either)
    "MISRA1D"
    "MISTAKE" # Faster than Ipopt
    # "MNISTS5" # Big and slow
    "MODEL" # Blows up (Ipopt reaches local infeasibility)
    "MPC16" # PANOC doesn't converge, large violation, eventually blows up
    "MRIBASIS" # Faster than Ipopt 
    "MSQRTA" # PANOC slow
    "MSQRTB"
    "MSS1"
    "MSS2" # Large but fast
    "MSS3" # Large but fast
    "MWRIGHT"
    "MUONSINE" # Violation stays high, then blows up
    "NASH" # Blows up immediately, why?
    # "NCVXQP1" # Error increases again, very big and slow
    "NELSON" # Inf gradient
    # "NET3" # NaN on first iteration (irregular)

    "NGONE" # Works with 10000 inner iterations
    "NINE12"
    "NINE5D" # Slow and steady at first, then error increases again and PANOC doesn't improve. Why are there sudden jumps in violation?
    "NINENEW" # PANOC gets close but takes time to converge from 2e-5 to 9e-6
    "NUFFIELD" # PANOC doesn't converge and blows up (without penalty blowing up)
    "NYSTROM5"
    "ODFITS"
    "ORTHREGB" 
    "ORTHREGC" # PANOC doesn't converge
    "ORTHRGDM" # PANOC doesn't converge
    "OSBORNE1" # Blows up
    "OSBORNE2" # Blows up

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

    # "SAROMM" # Huge and slow
    "SCOND1LS" # Faster than Ipopt (unconstrained)
    "STNQP2" # Gets close to 5e-5, doesn't make progress
    "STREG"
    "STREGNE"
    "STRTCHDVNE"
    "SUPERSIM"
    # "SVANBERG" # Gets close to 1e-4, but doesn't make progress afterwards
    "SWOPF" # PANOC doesn't converge

    # "SYNPOP24" # Huge
    "SYNTHES1" # Faster than Ipopt
    "SYNTHES2" # Faster than Ipopt
    "SYNTHES3" # Faster than Ipopt

    "TABLE1" # PANOC doesn't converge
    # "TARGUS" # Converges initially, but then doesn't get closer to the tolerance :(
    # "TOYSARAH" # PANOC doesn't converge, even with 100000 iterations :(
    # "TRIMLOSS" # PANOC not converging
    # "TRO11X3" # Huge constraint violation causes huge penalty → PANOC fails as well
    # "TRO3X3" # Idem (Ipopt is super fast, 0.2s)
    "TRUSPYR1" # 
    "TRUSPYR2" # converges really quickly with lbfgs mem = 100 (and Lipschitz update in ls)
    "TRY-B" # clear difference between update lipschitz in ls
    "TWIRISM1"
    # "TWIRIMD1" # PANOC doesn't converge
    # "TWIRIBG1" # (large) PANOC doesn't converge
    "TWOBARS"
    "VIBRBEAMNE" # Penalty blows up, works with ε₀ = 1e-5, Δ = 2, update_lipschitz_in_linesearch, blows up if update_lipschitz_in_linesearch = false
    # "VESUVIA" # really slow, worse with specialized_lbfgs
    "WACHBIEG"
    "WATER" # Faster than Ipopt
    "WOMFLET"
    "YFITNE" # Ipopt can't solve (not enough degrees of freedom)
    "YORKNET" # Penalty (and PANOC) blow up, converges with Δ=2, but slower than Ipopt and worse solution
    "ZAMB2-8" # Ipopt is faster
    "ZAMB2-9" # PANOC doesn't converge
    "ZAMB2-10" # PANOC doesn't converge
    "ZAMB2-11" # PANOC doesn't converge →→→→→ FIND A WAY TO DETECT WHEN THIS HAPPENS
    "ZANGWIL3"
    "ZECEVIC2"
    "ZECEVIC3"
    "ZECEVIC4"
    "ZIGZAG" # Penalty blows up 
    "ZY2"
)