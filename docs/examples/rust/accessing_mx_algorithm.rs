//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// Demonstration on how the algorithm of an MX function can be accessed and its operations can be transversed.

use casadi as ca;

// The OP_* operation codes are only reachable through the raw layer
fn op(c: fn(Vec<ca::rt::RV>) -> Result<Vec<ca::rt::RV>, ca::rt::Error>) -> i64 {
    ca::rt::out1(c(vec![]).unwrap())
}

fn main() {
    let op_const = op(ca::raw::constant_OP_CONST);
    let op_input = op(ca::raw::constant_OP_INPUT);
    let op_output = op(ca::raw::constant_OP_OUTPUT);
    let op_add = op(ca::raw::constant_OP_ADD);
    let op_mul = op(ca::raw::constant_OP_MUL);
    let op_mtimes = op(ca::raw::constant_OP_MTIMES);

    // Create a function
    let a = ca::MX::sym("a");
    let b = ca::MX::sym_with("b", 2, None);
    let c = ca::MX::sym_with("c", 2, 2);
    let f = ca::Function::named(
        "f",
        vec![&a, &b, &c],
        vec![3.0 * ca::mtimes(&c, &b) * &a + &b],
        vec!["a", "b", "c"],
        vec!["r"],
    );

    // Input values of the same dimensions as the above
    let input_val = vec![
        ca::DM::new(2.0),
        ca::DM::new(vec![3.0, 4.0]),
        ca::DM::new(vec![vec![5.0, 1.0], vec![8.0, 4.0]]),
    ];

    // Output values to be calculated of the same dimensions as the above
    let mut output_val = vec![ca::DM::zeros_with(2, 1)];

    // Work vector
    let mut work: Vec<ca::DM> = vec![ca::DM::default(); f.sz_w() as usize];

    // Loop over the algorithm
    for k in 0..f.n_instructions() {
        // Get the atomic operation
        let op = f.instruction_id(k);

        let o: Vec<usize> = f.instruction_output(k).iter().map(|&e| e as usize).collect();
        let i: Vec<usize> = f.instruction_input(k).iter().map(|&e| e as usize).collect();

        if op == op_const {
            let v = f.instruction_MX(k).to_DM();
            assert!(v.is_dense());
            println!("work[{}] = {}", o[0], v);
            work[o[0]] = v;
        } else if op == op_input {
            work[o[0]] = input_val[i[0]].clone();
            println!("work[{}] = input[{}]            ---> {}", o[0], i[0], work[o[0]]);
        } else if op == op_output {
            output_val[o[0]] = work[i[0]].clone();
            println!("output[{}] = work[{}]             ---> {}", o[0], i[0], output_val[o[0]]);
        } else if op == op_add {
            work[o[0]] = &work[i[0]] + &work[i[1]];
            println!("work[{}] = work[{}] + work[{}]      ---> {}", o[0], i[0], i[1], work[o[0]]);
        } else if op == op_mul {
            work[o[0]] = &work[i[0]] * &work[i[1]];
            println!("work[{}] = work[{}] * work[{}]        ---> {}", o[0], i[0], i[1], work[o[0]]);
        } else if op == op_mtimes {
            work[o[0]] = ca::mtimes(&work[i[1]], &work[i[2]]) + &work[i[0]];
            println!("work[{}] = work[{}] @ work[{}] + work[{}]        ---> {}", o[0], i[1], i[2], i[0], work[o[0]]);
        } else {
            let disp_in: Vec<String> = i.iter().map(|a| format!("work[{}]", a)).collect();
            let debug_str = ca::print_operator(&f.instruction_MX(k), &disp_in);
            panic!("Unknown operation: {} -- {}", op, debug_str);
            // When evaluating with MX, you may use an identity transformation:
            //let node = f.instruction_MX(k);
            //let args: Vec<ca::MX> = i.iter().map(|&ik| work[ik].clone()).collect();
            //let res = node.eval_mx(args);
            //for (ok, e) in o.iter().zip(res) { work[*ok] = e; }
        }
    }

    println!("------");
    println!("Evaluated {}", f);
    let expected: ca::DM = f.eval(input_val);
    println!("Expected:  [{}]", expected);
    println!("Got:       [{}]", output_val[0]);
}
