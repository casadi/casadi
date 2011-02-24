cmake_minimum_required(VERSION 2.0)

include_directories(Include)

set(CSPARSE_SRC
        Include/cs.h
        Source/cs_add Source/cs_amd Source/cs_chol Source/cs_cholsol Source/cs_counts Source/cs_cumsum 
        Source/cs_droptol Source/cs_dropzeros Source/cs_dupl Source/cs_entry 
        Source/cs_etree Source/cs_fkeep Source/cs_gaxpy Source/cs_happly Source/cs_house Source/cs_ipvec 
        Source/cs_lsolve Source/cs_ltsolve Source/cs_lu Source/cs_lusol Source/cs_util Source/cs_multiply 
        Source/cs_permute Source/cs_pinv Source/cs_post Source/cs_pvec Source/cs_qr Source/cs_qrsol 
        Source/cs_scatter Source/cs_schol Source/cs_sqr Source/cs_symperm Source/cs_tdfs Source/cs_malloc 
        Source/cs_transpose Source/cs_compress Source/cs_usolve Source/cs_utsolve Source/cs_scc 
        Source/cs_maxtrans Source/cs_dmperm Source/cs_updown Source/cs_print Source/cs_norm Source/cs_load 
        Source/cs_dfs Source/cs_reach Source/cs_spsolve Source/cs_ereach Source/cs_leaf Source/cs_randperm
)

add_library(csparse STATIC ${CSPARSE_SRC})
if(ENABLE_SHARED)
add_library(csparse_shared SHARED ${CSPARSE_SRC})
endif(ENABLE_SHARED)
