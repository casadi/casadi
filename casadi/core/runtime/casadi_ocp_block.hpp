//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//


// SYMBOL "ocp_block"
struct casadi_ocp_block {
    casadi_int offset_r;
    casadi_int offset_c;
    casadi_int rows;
    casadi_int cols;
};
// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"

// SYMBOL "unpack_ocp_blocks"
template<typename T1>
void casadi_unpack_ocp_blocks(casadi_ocp_block* blocks, const casadi_int* packed) {
    casadi_int i;
    casadi_int N = *packed++;
    for (i=0;i<N;++i) {
        blocks[i].offset_r = *packed++;
        blocks[i].offset_c = *packed++;
        blocks[i].rows = *packed++;
        blocks[i].cols = *packed++;
    }
}

// SYMBOL "ptr_ocp_block"
template<typename T1>
void casadi_ptr_ocp_block(casadi_int N, T1** vs, T1* v, const casadi_ocp_block* blocks, int eye) {
    casadi_int k, offset = 0;
    for (k=0;k<N;++k) {
        vs[k] = v+offset;
        if (eye) {
        offset += blocks[k].rows;
        } else {
        offset += blocks[k].rows*blocks[k].cols;
        }
    }
}
