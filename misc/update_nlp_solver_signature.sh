#!/bin/bash

# Input scheme
perl -pi -e 's/NLP_X_INIT/NLP_SOLVER_X0/g' `find *`
perl -pi -e 's/NLP_LBX/NLP_SOLVER_LBX/g' `find *`
perl -pi -e 's/NLP_UBX/NLP_SOLVER_UBX/g' `find *`
perl -pi -e 's/NLP_LBG/NLP_SOLVER_LBG/g' `find *`
perl -pi -e 's/NLP_UBG/NLP_SOLVER_UBG/g' `find *`
perl -pi -e 's/NLP_LAMBDA_INIT/NLP_SOLVER_LAM_G0/g' `find *`
perl -pi -e 's/NLP_P/NLP_SOLVER_P/g' `find *`
perl -pi -e 's/NLP_NUM_IN/NLP_SOLVER_NUM_IN/g' `find *`

#output scheme
perl -pi -e 's/NLP_X_OPT/NLP_SOLVER_X/g' `find *`
perl -pi -e 's/NLP_COST/NLP_SOLVER_F/g' `find *`
perl -pi -e 's/NLP_LAMBDA_G/NLP_SOLVER_LAM_G/g' `find *`
perl -pi -e 's/NLP_LAMBDA_X/NLP_SOLVER_LAM_X/g' `find *`
perl -pi -e 's/NLP_LAMBDA_P/NLP_SOLVER_LAM_P/g' `find *`
perl -pi -e 's/NLP_G/NLP_SOLVER_G/g' `find *`
perl -pi -e 's/NLP_NUM_OUT/NLP_SOLVER_NUM_OUT/g' `find *`

