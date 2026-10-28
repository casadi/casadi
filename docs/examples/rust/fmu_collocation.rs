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

// Example on how to formulate and solve an optimal control problem using direct collocation.
// The model is a semi-explicit DAE.
//
// Based on collocation example by Mario Zanon, Sebastien Gross 2012
// Reformulated to a semi-explict formulation and modern syntax by Joel Andersson, 2026

use casadi as ca;
use casadi::dict;
use std::path::Path;
use std::process::Command;

// Numeric value of a DaeBuilder attribute
fn to_f64(v: &ca::rt::Value) -> f64 {
    match v {
        ca::rt::Value::Float(x) => *x,
        ca::rt::Value::Int(i) => *i as f64,
        ca::rt::Value::Bool(b) => *b as i64 as f64,
        _ => panic!("not a number: {:?}", v),
    }
}

// Generate a CasADi FMU from a DaeBuilder instance
fn export_fmu(dae: &ca::DaeBuilder) {
    // Generate FMU
    let mut fmu_files = dae.export_fmu();
    println!("Generated files: {:?}", fmu_files);

    // Compile DLL
    let fmi_headers = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../external_packages/FMI-Standard-3.0/headers");
    let cfiles: Vec<&String> = fmu_files.keys().filter(|f| f.ends_with(".c")).collect();
    // Compile sources (in Linux)
    if cfg!(unix) {
        let sofile = dae.name() + ".so";
        let status = Command::new("gcc")
            .args(["--shared", "-fPIC"])
            .arg(format!("-I{}", fmi_headers.display()))
            .args(&cfiles)
            .args(["-o", &sofile])
            .status();
        match status {
            Ok(_) => {
                println!("Compiled {}", sofile);
                fmu_files.insert(sofile.clone(), ca::rt::Value::Str(format!("binaries/x86_64-linux/{}", sofile)));
            }
            Err(e) => println!("Failed to create FMU binary: {}", e),
        }
    }

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
}

fn main() {
    // Precompiled FMU?
    let fmu: Option<&str> = None; // Generate CasADi FMU
    // let fmu = Some("casadi_symbolic_fmu"); // FMU from CasADi, use serialized CasADi expressions (default)
    // let fmu = Some("casadi_blackbox_fmu"); // FMU from CasADi, do not use serialized CasADi expressions
    // let fmu = Some("dymola_fmu"); // FMU from Dymola, compiled from Modelica

    let dae = match fmu {
        None => {
            // Start with an empty DaeBuilder instance
            let dae = ca::DaeBuilder::new("Crane");

            // Hard coded constants
            let m_l = 1.0; // Mass of the load
            let m_c = 1.0; // Mass of the crane
            let g = 9.81; // Gravity
            let length = 1.0; // Rope/pole length
            let inverted = true; // Inverted pendulum on a cart rather than crane with hanging load

            // States
            let x_l = dae.add_3("xL", "local", dict! {"initial" => "exact", "start" => 0,
                "description" => "Horizontal position of the load"});
            let y = dae.add_3("y", "local", dict! {"initial" => "exact",
                "start" => if inverted { length } else { -length },
                "description" => "Vertical position of the load"});
            let x_c = dae.add_3("xC", "local", dict! {"initial" => "exact", "start" => 0,
                "description" => "Horizontal position of the crane"});
            let dx_l = dae.add_3("dxL", "local", dict! {"initial" => "exact", "start" => 0,
                "description" => "Horizontal velocity of the load"});
            let dy = dae.add_3("dy", "local", dict! {"initial" => "exact", "start" => 0,
                "description" => "Vertical velocity of the load"});
            let dx_c = dae.add_3("dxC", "local", dict! {"initial" => "exact", "start" => 0,
                "description" => "Horizontal velocity of the crane"});

            // Input
            let u = dae.add_3("u", "input", dict! {"min" => -2, "max" => 2, "start" => 0,
                "description" => "Horizontal force applied to the crane"});

            // Algebraic variables
            let t = dae.add_3("T", "local", dict! {"initial" => "approx", "start" => if inverted { g } else { -g },
                "description" => "Tension in the rope/pole"});

            // Ordinary differential equations
            let ddx_l = (&x_l - &x_c) / length * &t / m_l;
            let ddy = g - &y / length * &t / m_l;
            let ddx_c = ((&x_l - &x_c) / length * &t - &u) / m_c;
            dae.eq(&dae.der_with(&x_l), &dx_l);
            dae.eq(&dae.der_with(&y), &dy);
            dae.eq(&dae.der_with(&x_c), &dx_c);
            dae.eq(&dae.der_with(&dx_l), &ddx_l);
            dae.eq(&dae.der_with(&dy), &ddy);
            dae.eq(&dae.der_with(&dx_c), &ddx_c);

            // Original algebraic equation
            let res = ca::power(&y, 2.0) + ca::power(&(&x_l - &x_c), 2.0) - length * length;
            println!("Original residual equation: {}", res);

            // Index reduction: Differentiate with respect to time
            // In CasADi, the most efficient way to do index reduction it to use forward
            // mode AD, e.g. us using Jacobian-times-vector products (jtimes):
            // Using the chain rule we have: d(res)/dt = d(res)/dx * dx/dt
            let x_in = ca::vertcat(dae.inputs("x"));
            let ode = ca::vertcat(dae.outputs_with("ode"));
            let res1: ca::MX = ca::jtimes(&res, &x_in, &ode);
            let res2: ca::MX = ca::jtimes(&res1, &x_in, &ode);
            println!("Residual equation after index reduction: {}", res2);

            // Add algebraic equations
            dae.eq(0.0, &res2);

            // Simulation duration
            dae.set_start_time(0.0);
            dae.set_stop_time(5.0);

            // Generate FMU
            export_fmu(&dae);
            dae
        }
        // Use FMU generated from CasADi using the script above: Standard FMU + CasADi expressions + LS-DAE
        Some("casadi_symbolic_fmu") => ca::DaeBuilder::new_with("Crane", "./Crane.fmu", dict! {}),
        // Use FMU generated from CasADi using the script above: Standard FMU + LS-DAE
        Some("casadi_blackbox_fmu") => ca::DaeBuilder::new_with("Crane", "./Crane.fmu",
                                                                dict! {"enable_ls_serialization" => false}),
        // Use FMU compiled from Modelica using Dymola
        Some("dymola_fmu") => ca::DaeBuilder::new_with("crane", "../../../test/data/Crane.fmu", dict! {}),
        Some(other) => panic!("Unknown FMU option: {}", other),
    };

    // Display the model
    dae.disp_with(true);

    // DAE right-hand-side function
    let ffcn = dae.create_4("ffcn", vec!["x", "z", "u"], vec!["ode", "alg", "y"], dict! {});

    // Variable names
    let u = dae.u();
    let x = dae.x();
    let z = dae.z();
    let y = dae.y();

    // Refence location for crane and load
    let xref = 0.1;

    // Weight/reference vector: value for the listed names, 0 otherwise
    let pick = |names: &[String], sel: &[&str], val: f64| -> ca::DM {
        ca::DM::new(names.iter().map(|n| if sel.contains(&n.as_str()) { val } else { 0.0 }).collect::<Vec<f64>>())
    };

    // Quadratic cost at the beginning of control intervals
    let x_ref = pick(&x, &["xL", "xC"], xref);
    let y_ref = pick(&y, &["xL", "xC"], xref);
    let x_weights = pick(&x, &["xL", "xC"], 1.0);
    let y_weights = pick(&y, &["xL", "xC"], 1.0);

    // Quadratic cost at final time
    let xf_ref = pick(&x, &["xL", "xC"], xref);
    let yf_ref = pick(&y, &["xL", "xC"], xref);
    let xf_weights = pick(&x, &["xL", "xC", "dxL", "dy"], 1.0);
    let yf_weights = pick(&y, &["xL", "xC", "dxL", "dy"], 1.0);

    // Discretization
    let nk = 50; // Number of control intervals
    let deg = 4; // Degree of interpolating polynomial

    // Size of the finite elements
    let tf = 5.0;
    let h = tf / nk as f64;

    // We will use collocation discretization using Legendre roots, cf.
    // Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
    // by Lorenz Biegler (2010).
    // The roots can be queried from CasADi or looked up in the above textbook
    let tau_root = ca::collocation_points_with(deg, "radau");

    // We can query the coefficient for the interpolating polynomials using a helper function
    // C: Coefficients for calculating state derivatives
    // D: Coefficients for calculating state at the end of the interval
    // B: Coefficients for calculating a quadrature at the end of the interval
    let (c, d, _b) = ca::collocation_coeff(&tau_root);

    // We will construct an NLP incrementally, starting with no decision variables,
    // and constraints, and zero objective
    let mut w: Vec<ca::MX> = vec![];
    let mut g: Vec<ca::MX> = vec![];
    let mut j = ca::MX::new(0.0);

    // Solution trajectories
    let mut u_all: Vec<ca::MX> = vec![];
    let mut x_all: Vec<ca::MX> = vec![];
    let mut z_all: Vec<ca::MX> = vec![];
    let mut y_all: Vec<ca::MX> = vec![];

    // Time grids, for solution
    let mut ugrid: Vec<f64> = vec![];
    let mut xgrid: Vec<f64> = vec![];
    let mut zgrid: Vec<f64> = vec![];
    let mut ygrid: Vec<f64> = vec![];

    // Current time
    let mut tk = 0.0;

    // Add a variable for the initial state
    let mut xk = ca::MX::sym_with("X0", x.len() as i64, 1);
    x_all.push(xk.clone());
    w.push(xk.clone());
    xgrid.push(tk);
    let mut yk: Option<ca::MX> = None;
    // Loop over all control intervals, constructing the NLP and x, u trajectories
    for k in 0..nk {
        // Add contribution to the cost
        let dx = &xk - &x_ref;
        j += ca::dot(&dx, &(&x_weights * &dx));
        if k > 0 && !y.is_empty() {
            // Add contribution to the cost
            let dy = yk.as_ref().unwrap() - &y_ref;
            j += ca::dot(&dy, &(&y_weights * &dy));
        }
        // Symbolic expression for the control for the interval
        let uk = ca::MX::sym_with(format!("u{}", k), u.len() as i64, 1);
        u_all.push(uk.clone());
        w.push(uk.clone());
        ugrid.push(tk);
        // Declare variables at each collocation point
        let mut xc: Vec<ca::MX> = vec![];
        let mut zc: Vec<ca::MX> = vec![];
        let mut yc: Vec<ca::MX> = vec![];
        for jj in 1..=deg {
            // Differential state at collocation points
            let xkj = ca::MX::sym_with(format!("X{}_{}", k, jj), x.len() as i64, 1);
            xc.push(xkj.clone());
            x_all.push(xkj.clone());
            w.push(xkj);
            xgrid.push(tk + tau_root[jj - 1] * h);
            // Algebraic variables at collocation points
            let zkj = ca::MX::sym_with(format!("Z{}_{}", k, jj), z.len() as i64, 1);
            zc.push(zkj.clone());
            z_all.push(zkj.clone());
            w.push(zkj);
            zgrid.push(tk + tau_root[jj - 1] * h);
            // Output variables at collocation points
            let ykj = ca::MX::sym_with(format!("Y{}_{}", k, jj), y.len() as i64, 1);
            yc.push(ykj.clone());
            y_all.push(ykj.clone());
            w.push(ykj);
            ygrid.push(tk + tau_root[jj - 1] * h);
        }
        // Collect the collocation equations
        for jj in 0..deg {
            // Expression for the state derivative at the collocation point
            let mut xdot_j = c.at((0, jj)) * &xk;
            for r in 0..deg {
                xdot_j += c.at((r + 1, jj)) * &xc[r];
            }
            // Call the DAE right-hand-side function
            let (ode_j, alg_j, y_j): (ca::MX, ca::MX, ca::MX) = ffcn.eval((&xc[jj], &zc[jj], &uk));
            // Append collocation equations
            g.push(h * &ode_j - &xdot_j);
            // Add algebraic equations
            g.push(alg_j);
            // Add output equations
            g.push(&y_j - &yc[jj]);
        }
        // State at the end of the interval
        let mut xk_end = d.at(0) * &xk;
        for jj in 0..deg {
            xk_end = &xk_end + d.at(jj + 1) * &xc[jj];
        }
        // New variable for the state at the end of the interval
        tk += h;
        xk = ca::MX::sym_with(format!("X{}", k + 1), x.len() as i64, 1);
        x_all.push(xk.clone());
        w.push(xk.clone());
        xgrid.push(tk);
        // Enforce continuity
        g.push(&xk_end - &xk);
        // Output, algebraic variable at the end of the interval (only correct for Radau)
        let _zk = zc.last().unwrap().clone();
        yk = Some(yc.last().unwrap().clone());
    }

    // Include end of control interval
    ugrid.push(tk);

    // Add cost contribution from the final state
    let dx = &xk - &xf_ref;
    j += ca::dot(&dx, &(&xf_weights * &dx));
    if !y.is_empty() {
        // Add contribution to the cost
        let dy = yk.as_ref().unwrap() - &yf_ref;
        j += ca::dot(&dy, &(&yf_weights * &dy));
    }

    // Concatenate vectors
    let w = ca::vertcat(w);
    let g = ca::vertcat(g);

    // Allocate an NLP solver
    let nlp = dict! {"x" => &w, "f" => &j, "g" => &g};
    let opts = dict! {"ipopt.max_iter" => 100, "ipopt.tol" => 1e-4};
    let solver = ca::nlpsol_with("solver", "ipopt", &nlp, opts);

    // Create mappings from (X, Z, U) -> (W) and back
    let x_all = ca::horzcat(x_all);
    let z_all = ca::horzcat(z_all);
    let u_all = ca::horzcat(u_all);
    let y_all = ca::horzcat(y_all);
    let to_w = ca::Function::named("to_W", vec![&x_all, &z_all, &y_all, &u_all], vec![&w],
                                   vec!["X", "Z", "Y", "U"], vec!["W"]);
    let from_w = ca::Function::named("from_W", vec![&w], vec![&x_all, &z_all, &y_all, &u_all],
                                     vec!["W"], vec!["X", "Z", "Y", "U"]);

    // Variable bounds, initial guess: Get data structures of correct size
    type Traj = (ca::DM, ca::DM, ca::DM, ca::DM);
    let (mut x0, mut z0, mut y0, mut u0): Traj = from_w.eval([0.0]);
    let (mut lb_x, mut lb_z, mut lb_y, mut lb_u): Traj = from_w.eval([f64::NEG_INFINITY]);
    let (mut ub_x, mut ub_z, mut ub_y, mut ub_u): Traj = from_w.eval([f64::INFINITY]);

    // Initial guess, bounds for control trajectory
    let uinit = dae.get(&u);
    let umin: Vec<f64> = dae.min(&u);
    let umax: Vec<f64> = dae.max(&u);
    for (i, v0) in uinit.iter().map(to_f64).enumerate() {
        u0.set_at((i, ..), v0);
        lb_u.set_at((i, ..), umin[i]);
        ub_u.set_at((i, ..), umax[i]);
    }

    // Initial guess, bounds for state trajectory
    let xinit = dae.get(&x);
    let xmin: Vec<f64> = dae.min(&x);
    let xmax: Vec<f64> = dae.max(&x);
    for (i, v0) in xinit.iter().map(to_f64).enumerate() {
        x0.set_at((i, ..), v0);
        lb_x.set_at((i, ..), xmin[i]);
        lb_x.set_at((i, 0), v0);
        ub_x.set_at((i, ..), xmax[i]);
        ub_x.set_at((i, 0), v0);
    }

    // Initial guess for algebraic variable trajectory
    let zinit = dae.get(&z);
    let zmin: Vec<f64> = dae.min(&z);
    let zmax: Vec<f64> = dae.max(&z);
    for (i, v0) in zinit.iter().map(to_f64).enumerate() {
        z0.set_at((i, ..), v0);
        lb_z.set_at((i, ..), zmin[i]);
        ub_z.set_at((i, ..), zmax[i]);
    }

    // Initial guess for output trajectory
    let yinit = dae.get(&y);
    let ymin: Vec<f64> = dae.min(&y);
    let ymax: Vec<f64> = dae.max(&y);
    for (i, v0) in yinit.iter().map(to_f64).enumerate() {
        y0.set_at((i, ..), v0);
        lb_y.set_at((i, ..), ymin[i]);
        ub_y.set_at((i, ..), ymax[i]);
    }

    // Solve the NLP
    let w0: ca::DM = to_w.eval((&x0, &z0, &y0, &u0));
    let lb_w: ca::DM = to_w.eval((&lb_x, &lb_z, &lb_y, &lb_u));
    let ub_w: ca::DM = to_w.eval((&ub_x, &ub_z, &ub_y, &ub_u));
    let sol = solver.call(&dict! {"x0" => &w0, "lbx" => &lb_w, "ubx" => &ub_w, "lbg" => 0, "ubg" => 0}).unwrap();
    let (xopt, zopt, yopt, uopt): Traj = from_w.eval([&sol["x"]]);

    // Variables to plot: print the trajectories at a few grid points instead
    println!("Optimal cost: {}", sol["f"]);
    let show = |label: &str, i: usize, n: &str, grid: &[f64], m: &ca::DM| {
        let v = m.at((i, ..)).nonzeros();
        let ks = [0, v.len() / 2, v.len() - 1];
        let pts: Vec<String> = ks.iter().map(|&k| format!("t={:.4}: {:.6}", grid[k], v[k])).collect();
        println!("{}[{}]: {}: {}", label, i, n, pts.join(", "));
    };
    for (i, n) in x.iter().enumerate() { show("x", i, n, &xgrid, &xopt); }
    for (i, n) in z.iter().enumerate() { show("z", i, n, &zgrid, &zopt); }
    for (i, n) in y.iter().enumerate() { show("y", i, n, &ygrid, &yopt); }
    for (i, n) in u.iter().enumerate() { show("u", i, n, &ugrid[..ugrid.len() - 1], &uopt); }
}
