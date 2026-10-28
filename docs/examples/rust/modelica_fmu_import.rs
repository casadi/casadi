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

use casadi as ca;
use std::process::Command;

fn main() {
    // Consider the following simple Modelica file
    println!("{}", std::fs::read_to_string("rocket.mo").expect("rocket.mo not found in the working directory"));

    // Different tools exist to compile Modelica models like these into FMUs.
    // The following code shows how we could do it using the open-source tool
    // OpenModelica, here by running a script with its compiler omc,
    // cf. the OpenModelica User's Guide.
    let mos = r#"if not loadFile("rocket.mo") then
  print("Modelica compilation failed: " + getErrorString());
  exit(1);
end if;
setDebugFlags("-disableDirectionalDerivatives");
fmu_file := buildModelFMU(rocket, version="2.0", fmuType="me");
print("translateModelFMU warnings:\n" + getErrorString() + "\n");
print(fmu_file);
"#;
    std::fs::write("rocket_build.mos", mos).unwrap();
    let out = Command::new("omc").arg("rocket_build.mos").output().expect("omc (OpenModelica) not found");
    let out = String::from_utf8_lossy(&out.stdout);
    if !out.trim_end().ends_with(".fmu") {
        panic!("FMU generation failed: {}", out);
    }
    println!("{}", out);

    // Regardless of how the FMU was obtained, we must unzip it before CasADi can process it
    // If there is an existing directory with the same name it is better to explicitly
    // delete it to prevent potential memory corruption on overwriting.
    let unzipped_path = std::env::current_dir().unwrap().join("rocket_fmu_unzipped");
    if unzipped_path.is_dir() {
        std::fs::remove_dir_all(&unzipped_path).unwrap();
    }
    let ok = Command::new("unzip").args(["-q", "rocket.fmu", "-d"]).arg(&unzipped_path)
        .status().expect("unzip failed");
    assert!(ok.success(), "unzip failed");

    // Create a CasADi/DaeBuilder instance from the unzipped FMU
    let dae = ca::DaeBuilder::new_with("rocket", unzipped_path.to_str().unwrap(), None);
    dae.disp_with(true);

    // Get state vector, initial conditions, bounds
    let x = dae.x();
    let lbx: Vec<f64> = dae.min(&x);
    let ubx: Vec<f64> = dae.max(&x);
    let x0 = dae.start(&x);
    println!("x:  {:?}", x);
    println!("lbx:  {:?}", lbx);
    println!("ubx:  {:?}", ubx);
    println!("x0:  {:?}", x0);

    // Get control vector, initial conditions, bounds
    let u = dae.u();
    let lbu: Vec<f64> = dae.min(&u);
    let ubu: Vec<f64> = dae.max(&u);
    let u0 = dae.start(&u);
    println!("u:  {:?}", u);
    println!("lbu:  {:?}", lbu);
    println!("ubu:  {:?}", ubu);
    println!("u0:  {:?}", u0);

    // Let's create a CasADi function for evaluating the ODE right-hand-side.
    // We only need to specify the expressions that are varying:
    let f = dae.create_io("f", vec!["x", "u"], vec!["ode"]);

    // This is a standard CasADi function that can be embedded into CasADi expressions
    println!("{}", f);

    // We can evaluate it and calculate derivatives as with other CasADi functions
    let ode0: ca::DM = f.eval((&x0, &u0));
    println!("ode0:  {}", ode0);

    // Analytic first derivatives with sparsities are also available
    let jac = f.factory("jac_f", vec!["x", "u"], vec!["jac:ode:x", "jac:ode:u"]);
    let (jac_ode_x0, jac_ode_u0): (ca::DM, ca::DM) = jac.eval((&x0, &u0));
    println!("jac_ode_x0:  {}", jac_ode_x0);
    println!("jac_ode_u0:  {}", jac_ode_u0);

    // Code generation for the functions are currently work in progress
    f.generate_with("fmu_codegen", None);
}
