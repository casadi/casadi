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

// Simulating a bouncing ball with DaeBuilder and event handling
// Joel Andersson, 2025

use casadi as ca;
use casadi::dict;
use std::path::Path;
use std::process::Command;

fn main() {
    // Start with an empty DaeBuilder instance
    let dae = ca::DaeBuilder::new("bouncing_ball");

    // Model variables
    let _t = dae.add_3("t", "independent", dict! {});
    let h = dae.add_3("h", "output", dict! {"start" => 5, "initial" => "exact"});
    let v = dae.add_3("v", "output", dict! {"start" => 0, "initial" => "exact"});

    // Dynamic equations
    dae.eq(&dae.der_with(&h), &v);
    dae.eq(&dae.der_with(&v), -9.81);
    dae.disp_with(true);

    // Event dynamics: When h < 0, reinitialize v to -0.8*v
    let pre_v: ca::MX = dae.pre(&v);
    dae.when(&ca::lt(&h, 0.0), vec![dae.reinit("v", &(-0.8 * &pre_v))]);
    dae.disp_with(true);

    // Default experiment
    dae.set_start_time(0.0);
    dae.set_stop_time(7.0);

    // Simulate
    let tgrid: ca::DM = ca::linspace(dae.start_time(), dae.stop_time(), 100);
    let tgrid = tgrid.nonzeros();
    let sim = ca::integrator_with("sim", "cvodes", &dae.create(), 0.0, &tgrid,
                                  dict! {"transition" => dae.transition()});
    let simres = sim.call(&dict! {"x0" => dae.start(dae.x())}).unwrap();

    // Visualize the solution: print the height instead
    let hsol = simres["xf"].at((0, ..)).nonzeros();
    for k in (0..tgrid.len()).step_by(11) {
        println!("t = {:.4}: h = {:.6}", tgrid[k], hsol[k]);
    }

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
