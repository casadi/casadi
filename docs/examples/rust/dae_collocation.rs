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
// @author: Mario Zanon and Sebastien Gross, KU Leuven 2012

use casadi as ca;
use casadi::dict;

const INF: f64 = f64::INFINITY;

fn main() {
    // -----------------------------------------------------------------------------
    // Collocation setup
    // -----------------------------------------------------------------------------
    let nicp = 1; // Number of (intermediate) collocation points per control interval

    let xref = 0.1; // chariot reference

    let l = 1.0; //- -> crane, + -> pendulum
    let m = 1.0;
    let mm = 1.0;
    let g = 9.81;
    let tf = 5.0;
    let nk = 50;
    let ndstate = 6;
    let nastate = 1;
    let _ninput = 1;

    // Degree of interpolating polynomial
    let deg = 4;
    // Radau collocation points
    let cp = "radau";
    // Size of the finite elements
    let h = tf / nk as f64 / nicp as f64;

    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; deg + 1]; deg + 1];
    // Coefficients of the continuity equation
    let mut d = vec![0.0; deg + 1];

    // Collocation point
    let tau = ca::SX::sym("tau");

    // All collocation time points
    let mut tau_root = vec![0.0];
    tau_root.extend(ca::collocation_points_with(deg as i64, cp));

    let mut t = vec![vec![0.0; deg + 1]; nk];
    for i in 0..nk {
        for j in 0..=deg {
            t[i][j] = h * (i as f64 + tau_root[j]);
        }
    }

    // For all collocation points: eq 10.4 or 10.17 in Biegler's book
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    for j in 0..=deg {
        let mut ll = ca::SX::new(1.0);
        for j2 in 0..=deg {
            if j2 != j {
                ll *= (&tau - tau_root[j2]) / (tau_root[j] - tau_root[j2]);
            }
        }

        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        let lfcn = ca::Function::new("lfcn", vec![&tau], vec![&ll]);
        d[j] = lfcn.eval::<ca::DM>([1.0]).to_double();

        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let tfcn = ca::Function::new("tfcn", vec![&tau], vec![ca::tangent(&ll, &tau)]);
        for j2 in 0..=deg {
            c[j][j2] = tfcn.eval::<ca::DM>([tau_root[j2]]).to_double();
        }
    }

    // -----------------------------------------------------------------------------
    // Model setup
    // -----------------------------------------------------------------------------
    // Declare variables (use scalar graph)
    let tt = ca::SX::sym("t"); // time
    let u = ca::SX::sym("u"); // control
    let mut xd = ca::SX::sym_with("xd", ndstate, None); // differential state
    let xa = ca::SX::sym_with("xa", nastate, None); // algebraic state
    let xddot = ca::SX::sym_with("xdot", ndstate, None); // differential state time derivative
    let p = ca::SX::sym_with("p", 0, None); // parameters

    let x = ca::SX::sym("x");
    let y = ca::SX::sym("y");
    let w = ca::SX::sym("w");

    let dx = ca::SX::sym("dx");
    let dy = ca::SX::sym("dy");
    let dw = ca::SX::sym("dw");

    let res = ca::vertcat(vec![
        xddot.at(0) - &dx,
        xddot.at(1) - &dy,
        xddot.at(2) - &dw,
        m * xddot.at(3) + (&x - &w) * &xa,
        m * xddot.at(4) + &y * &xa - g * m,
        mm * xddot.at(5) + (&w - &x) * &xa + &u,
        (&x - &w) * (xddot.at(3) - xddot.at(5)) + &y * xddot.at(4) + &dy * &dy + (&dx - &dw) * (&dx - &dw),
    ]);

    xd.set_at(0, &x);
    xd.set_at(1, &y);
    xd.set_at(2, &w);
    xd.set_at(3, &dx);
    xd.set_at(4, &dy);
    xd.set_at(5, &dw);

    // System dynamics (implicit formulation)
    let ffcn = ca::Function::new("ffcn", vec![&tt, &xddot, &xd, &xa, &u, &p], vec![&res]);

    // Objective function
    let mayer_term = ca::Function::new(
        "mayer",
        vec![&tt, &xd, &xa, &u, &p],
        vec![(&x - xref) * (&x - xref) + (&w - xref) * (&w - xref) + &dx * &dx + &dy * &dy],
    );
    let lagrange_term = ca::Function::new(
        "lagrange",
        vec![&tt, &xd, &xa, &u, &p],
        vec![(&x - xref) * (&x - xref) + (&w - xref) * (&w - xref)],
    );

    // Control bounds
    let u_min = vec![-2.0];
    let u_max = vec![2.0];
    let u_init = vec![vec![0.0]; nk * nicp * (deg + 1)]; // needs to be specified for every time interval (even though it stays constant)

    // Differential state bounds
    //Path bounds
    let xd_min = vec![-INF; 6];
    let xd_max = vec![INF; 6];
    //Initial bounds
    let xdi_min = vec![0.0, l, 0.0, 0.0, 0.0, 0.0];
    let xdi_max = vec![0.0, l, 0.0, 0.0, 0.0, 0.0];
    //Final bounds
    let xdf_min = vec![-INF; 6];
    let xdf_max = vec![INF; 6];

    //Initial guess for differential states
    let xd_init = vec![vec![0.0, l, 0.0, 0.0, 0.0, 0.0]; nk * nicp * (deg + 1)]; // needs to be specified for every time interval

    // Algebraic state bounds and initial guess
    let xa_min = vec![-INF];
    let xa_max = vec![INF];
    let _xai_min = vec![-INF];
    let _xai_max = vec![INF];
    let _xaf_min = vec![-INF];
    let _xaf_max = vec![INF];
    let xa_init = vec![vec![ca::sign(&ca::DM::new(l)).to_double() * 9.81]; nk * nicp * (deg + 1)];

    // Parameter bounds and initial guess
    let p_min: Vec<f64> = vec![];
    let p_max: Vec<f64> = vec![];
    let p_init: Vec<f64> = vec![];

    // -----------------------------------------------------------------------------
    // Constraints setup
    // -----------------------------------------------------------------------------
    // Initial constraint
    let ic_min: Vec<f64> = vec![];
    let ic_max: Vec<f64> = vec![];
    let ic = ca::SX::default();
    let icfcn = ca::Function::new("icfcn", vec![&tt, &xd, &xa, &u, &p], vec![&ic]);
    // Path constraint
    let pc_min: Vec<f64> = vec![];
    let pc_max: Vec<f64> = vec![];
    let pc = ca::SX::default();
    let pcfcn = ca::Function::new("pcfcn", vec![&tt, &xd, &xa, &u, &p], vec![&pc]);
    // Final constraint
    let fc_min: Vec<f64> = vec![];
    let fc_max: Vec<f64> = vec![];
    let fc = ca::SX::default();
    let fcfcn = ca::Function::new("fcfcn", vec![&tt, &xd, &xa, &u, &p], vec![&fc]);

    // -----------------------------------------------------------------------------
    // NLP setup
    // -----------------------------------------------------------------------------
    // Dimensions of the problem
    let nx = (xd.nnz() + xa.nnz()) as usize; // total number of states
    let ndiff = xd.nnz() as usize; // number of differential states
    let nalg = xa.nnz() as usize; // number of algebraic states
    let nu = u.nnz() as usize; // number of controls
    let np = p.nnz() as usize; // number of parameters

    // Total number of variables
    let nxd = nicp * nk * (deg + 1) * ndiff; // Collocated differential states
    let nxa = nicp * nk * deg * nalg; // Collocated algebraic states
    let nuu = nk * nu; // Parametrized controls
    let nxf = ndiff; // Final state (only the differential states)
    let nv = nxd + nxa + nuu + nxf + np;

    // NLP variable vector
    let v = ca::MX::sym_with("V", nv as i64, None);

    // All variables with bounds and initial guess
    let mut vars_lb = vec![0.0; nv];
    let mut vars_ub = vec![0.0; nv];
    let mut vars_init = vec![0.0; nv];
    let mut offset = 0;

    // Get the parameters
    let pp = v.at(offset..offset + np);
    vars_init[offset..offset + np].copy_from_slice(&p_init);
    vars_lb[offset..offset + np].copy_from_slice(&p_min);
    vars_ub[offset..offset + np].copy_from_slice(&p_max);
    offset += np;

    // Get collocated states and parametrized control
    let mut xdv = vec![vec![vec![ca::MX::default(); deg + 1]; nicp]; nk + 1];
    let mut xav = vec![vec![vec![ca::MX::default(); deg]; nicp]; nk];
    let mut uv = vec![ca::MX::default(); nk];
    let mut index = 0;
    for k in 0..nk {
        // Collocated states
        for i in 0..nicp {
            for j in 0..=deg {
                // Get the expression for the state vector
                xdv[k][i][j] = v.at(offset..offset + ndiff);
                if j != 0 {
                    xav[k][i][j - 1] = v.at(offset + ndiff..offset + ndiff + nalg);
                }
                // Add the initial condition
                index = (deg + 1) * (nicp * k + i) + j;
                if k == 0 && j == 0 && i == 0 {
                    vars_init[offset..offset + ndiff].copy_from_slice(&xd_init[index]);

                    vars_lb[offset..offset + ndiff].copy_from_slice(&xdi_min);
                    vars_ub[offset..offset + ndiff].copy_from_slice(&xdi_max);
                    offset += ndiff;
                } else if j != 0 {
                    vars_init[offset..offset + nx]
                        .copy_from_slice(&[&xd_init[index][..], &xa_init[index][..]].concat());

                    vars_lb[offset..offset + nx].copy_from_slice(&[&xd_min[..], &xa_min[..]].concat());
                    vars_ub[offset..offset + nx].copy_from_slice(&[&xd_max[..], &xa_max[..]].concat());
                    offset += nx;
                } else {
                    vars_init[offset..offset + ndiff].copy_from_slice(&xd_init[index]);

                    vars_lb[offset..offset + ndiff].copy_from_slice(&xd_min);
                    vars_ub[offset..offset + ndiff].copy_from_slice(&xd_max);
                    offset += ndiff;
                }
            }
        }

        // Parametrized controls
        uv[k] = v.at(offset..offset + nu);
        vars_lb[offset..offset + nu].copy_from_slice(&u_min);
        vars_ub[offset..offset + nu].copy_from_slice(&u_max);
        vars_init[offset..offset + nu].copy_from_slice(&u_init[index]);
        offset += nu;
    }

    // State at end time
    xdv[nk][0][0] = v.at(offset..offset + ndiff);
    vars_lb[offset..offset + ndiff].copy_from_slice(&xdf_min);
    vars_ub[offset..offset + ndiff].copy_from_slice(&xdf_max);
    vars_init[offset..offset + ndiff].copy_from_slice(xd_init.last().unwrap());
    offset += ndiff;
    assert_eq!(offset, nv);

    // Constraint function for the NLP
    let mut gv: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // Initial constraints
    let ick: ca::MX = icfcn.eval((0.0, &xdv[0][0][0], &xav[0][0][0], &uv[0], &pp));
    gv.push(ick);
    lbg.extend(&ic_min);
    ubg.extend(&ic_max);

    // For all finite elements
    let (mut k, mut i, mut j) = (0, 0, 0);
    for k_ in 0..nk {
        k = k_;
        for i_ in 0..nicp {
            i = i_;
            // For all collocation points
            for j_ in 1..=deg {
                j = j_;
                // Get an expression for the state derivative at the collocation point
                let mut xp_jk = ca::MX::new(0.0);
                for j2 in 0..=deg {
                    xp_jk += c[j2][j] * &xdv[k][i][j2]; // get the time derivative of the differential states (eq 10.19b)
                }

                // Add collocation equations to the NLP
                let fk: ca::MX = ffcn.eval((0.0, &xp_jk / h, &xdv[k][i][j], &xav[k][i][j - 1], &uv[k], &pp));
                gv.push(fk.at(..ndiff)); // impose system dynamics (for the differential states (eq 10.19b))
                lbg.extend(vec![0.0; ndiff]); // equality constraints
                ubg.extend(vec![0.0; ndiff]); // equality constraints
                gv.push(fk.at(ndiff..)); // impose system dynamics (for the algebraic states (eq 10.19b))
                lbg.extend(vec![0.0; nalg]); // equality constraints
                ubg.extend(vec![0.0; nalg]); // equality constraints

                //  Evaluate the path constraint function
                let pck: ca::MX = pcfcn.eval((0.0, &xdv[k][i][j], &xav[k][i][j - 1], &uv[k], &pp));

                gv.push(pck);
                lbg.extend(&pc_min);
                ubg.extend(&pc_max);
            }

            // Get an expression for the state at the end of the finite element
            let mut xf_k = ca::MX::new(0.0);
            for j in 0..=deg {
                xf_k += d[j] * &xdv[k][i][j];
            }

            // Add continuity equation to NLP
            if i == nicp - 1 {
                gv.push(&xdv[k + 1][0][0] - &xf_k);
            } else {
                gv.push(&xdv[k][i + 1][0] - &xf_k);
            }

            lbg.extend(vec![0.0; ndiff]);
            ubg.extend(vec![0.0; ndiff]);
        }
    }

    // Periodicity constraints
    //   none

    // Final constraints (Const, dConst, ConstQ)
    let fck: ca::MX = fcfcn.eval((0.0, &xdv[k][i][j], &xav[k][i][j - 1], &uv[k], &pp));
    gv.push(fck);
    lbg.extend(&fc_min);
    ubg.extend(&fc_max);

    // Objective function of the NLP
    //Implement Mayer term
    let mut obj_total = ca::MX::new(0.0);
    let obj: ca::MX = mayer_term.eval((0.0, &xdv[k][i][j], &xav[k][i][j - 1], &uv[k], &pp));
    obj_total += obj;

    // Implement Lagrange term
    let l_dot_at_tau_root = ca::DM::new(c.clone()).T();
    let l_at_one = ca::DM::new(d.clone());

    let ld_inv = ca::inv(&l_dot_at_tau_root.at((1.., 1..)));
    let _ld0 = l_dot_at_tau_root.at((1.., 0));
    let mut lagrange = ca::MX::new(0.0);
    for k in 0..nk {
        for i in 0..nicp {
            let terms: Vec<ca::MX> =
                (1..=deg).map(|j| lagrange_term.eval((0.0, &xdv[k][i][j], &xav[k][i][j - 1], &uv[k], &pp))).collect();
            let dqs = h * ca::veccat(&terms);
            let qs = ca::mtimes(&ld_inv, &dqs);
            let m = ca::mtimes(&qs.T(), &l_at_one.at(1..));
            lagrange += m;
        }
    }

    obj_total += lagrange;

    // NLP
    let nlp = dict! {"x" => &v, "f" => &obj_total, "g" => ca::vertcat(&gv)};

    // ----
    // SOLVE THE NLP
    // ----

    // NLP solver options
    let opts = dict! {
        "expand" => true,
        "ipopt.max_iter" => 1000,
        "ipopt.tol" => 1e-4,
        "ipopt.linear_solver" => "ma27",
    };

    // Allocate an NLP solver
    let solver = ca::nlpsol_with("solver", "ipopt", &nlp, &opts);
    let arg = dict! {
        // Initial condition
        "x0" => &vars_init,

        // Bounds on x
        "lbx" => &vars_lb,
        "ubx" => &vars_ub,

        // Bounds on g
        "lbg" => &lbg,
        "ubg" => &ubg,
    };

    // Solve the problem
    let res = solver.call(&arg).unwrap();

    // Print the optimal cost
    println!("optimal cost:  {}", res["f"].to_double());

    // Retrieve the solution
    let v_opt = res["x"].nonzeros();

    // ----
    // RETRIEVE THE SOLUTION
    // ----
    let npts = (deg + 1) * nicp * nk + 1;
    let mut xd_opt = vec![vec![0.0; npts]; ndiff];
    let mut xa_opt = vec![vec![0.0; deg * nicp * nk]; nalg];
    let mut u_opt = vec![vec![0.0; npts]; nu];
    let mut offset = 0;
    let mut offset2 = 0;
    let mut offset3 = 0;
    let mut offset4 = 0;

    for _k in 0..nk {
        for _i in 0..nicp {
            for j in 0..=deg {
                for r in 0..ndiff {
                    xd_opt[r][offset2] = v_opt[offset + r];
                }
                offset2 += 1;
                offset += ndiff;
                if j != 0 {
                    for r in 0..nalg {
                        xa_opt[r][offset4] = v_opt[offset + r];
                    }
                    offset4 += 1;
                    offset += nalg;
                }
            }
        }
        let utemp = &v_opt[offset..offset + nu];
        for _i in 0..nicp {
            for _j in 0..=deg {
                for r in 0..nu {
                    u_opt[r][offset3] = utemp[r];
                }
                offset3 += 1;
            }
        }
        offset += nu;
    }

    for r in 0..ndiff {
        xd_opt[r][npts - 1] = v_opt[offset + r];
    }

    // The algebraic states are not defined at the first collocation point of the finite elements:
    // with the polynomials we compute them at that point
    let mut da = vec![0.0; deg];
    for j in 1..=deg {
        // Lagrange polynomials for the algebraic states: exclude the first point
        let mut la = ca::SX::new(1.0);
        for j2 in 1..=deg {
            if j2 != j {
                la *= (&tau - tau_root[j2]) / (tau_root[j] - tau_root[j2]);
            }
        }
        let lafcn = ca::Function::new("lafcn", vec![&tau], vec![&la]);
        da[j - 1] = lafcn.eval::<ca::DM>([tau_root[0]]).to_double();
    }

    let mut xa_plt = vec![vec![0.0; npts]; nalg];
    let mut offset4 = 0;
    let mut offset5 = 0;
    for _k in 0..nk {
        for _i in 0..nicp {
            for j in 0..=deg {
                if j != 0 {
                    for r in 0..nalg {
                        xa_plt[r][offset5] = xa_opt[r][offset4];
                    }
                    offset4 += 1;
                    offset5 += 1;
                } else {
                    for r in 0..nalg {
                        xa_plt[r][offset5] = (0..deg).map(|j| da[j] * xa_opt[r][offset4 + j]).sum();
                    }
                    offset5 += 1;
                }
            }
        }
    }

    for r in 0..nalg {
        xa_plt[r][npts - 1] = xa_plt[r][npts - 2];
    }

    let tg: Vec<f64> = tau_root.iter().map(|tr| tr * h).collect();
    let mut tgrid: Vec<f64> = vec![];
    for k in 0..nk * nicp {
        let t0 = if k == 0 { 0.0 } else { *tgrid.last().unwrap() };
        tgrid.extend(tg.iter().map(|x| t0 + x));
    }
    tgrid.push(*tgrid.last().unwrap());

    // Print the results
    let every = |v: &[f64]| v.iter().step_by(deg + 1).map(|x| format!("{:.6}", x)).collect::<Vec<_>>().join(", ");
    println!("t:      [{}]", every(&tgrid));
    println!("x:      [{}]", every(&xd_opt[0]));
    println!("y:      [{}]", every(&xd_opt[1]));
    println!("w:      [{}]", every(&xd_opt[2]));
    println!("u:      [{}]", every(&u_opt[0]));
    println!("lambda: [{}]", every(&xa_plt[0]));
}
