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

// Example how to export FMUs from CasADi
// Joel Andersson, joel@jaeandersson.com

use casadi as ca;
use casadi::dict;
use std::path::Path;
use std::process::Command;

fn main() {
    // Start with an empty DaeBuilder instance
    let dae = ca::DaeBuilder::new("vdp");

    // States
    let _t = dae.add_causality("t", "independent");
    let x1 = dae.add_causality_with("x1", "output", dict! {"start" => 1});
    let x2 = dae.add_causality_with("x2", "output", dict! {"start" => 0});
    let q = dae.add("q");

    // Input
    let u = dae.add_causality("u", "input");

    // Set ODE right-hand-sides
    dae.eq(&dae.der_with(&x1), &((1.0 - &x2 * &x2) * &x1 - &x2 + &u));
    dae.eq(&dae.der_with(&x2), &x1);
    dae.eq(&dae.der_with(&q), &(ca::power(&x1, 2.0) + ca::power(&x2, 2.0) + ca::power(&u, 2.0)));

    // Add bounds
    dae.set_min("u", -0.75);
    dae.set_max("u", 1.0);
    dae.set_start("u", 0.5);

    // Print DAE
    dae.disp_with(true);

    // Export FMU
    let mut fmu_files = dae.export_fmu();
    println!("Generated files: {:?}", fmu_files);

    // Compile DLL
    let fmi_headers = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../external_packages/FMI-Standard-3.0/headers");
    let cfiles: Vec<&String> = fmu_files.keys().filter(|f| f.ends_with(".c")).collect();
    let sofile = dae.name() + ".so";
    Command::new("gcc")
        .args(["--shared", "-fPIC"])
        .arg(format!("-I{}", fmi_headers.display()))
        .args(&cfiles)
        .args(["-o", &sofile])
        .status()
        .expect("gcc failed");
    println!("Compiled {}", sofile);
    fmu_files.insert(sofile.clone(), ca::rt::Value::Str(format!("binaries/x86_64-linux/{}", sofile)));

    // Package into an FMU: stage the files under their archive names, then zip
    let fmuname = dae.name() + ".fmu";
    let stage = std::env::temp_dir().join(format!("{}_fmu", dae.name()));
    let _ = std::fs::remove_dir_all(&stage);
    for (f, arcname) in &fmu_files {
        let ca::rt::Value::Str(arcname) = arcname else { panic!("archive name expected") };
        let dest = stage.join(arcname);
        std::fs::create_dir_all(dest.parent().unwrap()).unwrap();
        std::fs::copy(f, &dest).unwrap();
        std::fs::remove_file(f).unwrap();
    }
    let fmupath = std::env::current_dir().unwrap().join(&fmuname);
    let _ = std::fs::remove_file(&fmupath);
    Command::new("zip").args(["-q", "-r"]).arg(&fmupath).arg(".").current_dir(&stage)
        .status().expect("zip failed");
    std::fs::remove_dir_all(&stage).unwrap();
    println!("Created FMU: {}", fmuname);

    // Load the FMU in FMPy: Python only
    println!("FMPy not available from Rust. Skipping FMU simulation.");
}
