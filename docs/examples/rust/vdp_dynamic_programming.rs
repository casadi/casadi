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
// The Python original uses numpy only (no CasADi); this port uses plain Rust vectors.

const INF: f64 = f64::INFINITY;

// System dynamics, can be called with matrices (here: elementwise on a grid)
fn f(x1: f64, x2: f64, u: f64) -> (f64, f64, f64) {
    let x1_dot = (1.0 - x2 * x2) * x1 - x2 + u;
    let x2_dot = x1;
    let q_dot = x1 * x1 + x2 * x2 + u * u;
    (x1_dot, x2_dot, q_dot)
}

fn linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n).map(|i| a + (b - a) * i as f64 / (n - 1) as f64).collect()
}

fn main() {
    // End time
    let t_end = 10.0;

    // Number of control intervals
    let n = 20;

    // Number of Runge-Kutta 4 steps per interval and step size
    let nk = 20;
    let dt = t_end / (n * nk) as f64;

    // Number of discrete control values
    let nu = 101;

    // Number of discrete state values
    let nx = 101;

    // Control enumeration
    let u_grid = linspace(-1.0, 1.0, nu);

    // State space enumeration (meshgrid: element [i2*nx + i1] holds (x1[i1], x2[i2]))
    let x1 = linspace(-1.0, 1.0, nx);
    let x2 = linspace(-1.0, 1.0, nx);
    let idx = |i2: usize, i1: usize| i2 * nx + i1;

    // For each control action and state, precalculate next state and stage cost
    let mut stage_j: Vec<Vec<f64>> = vec![];
    let mut next_x1: Vec<Vec<usize>> = vec![];
    let mut next_x2: Vec<Vec<usize>> = vec![];
    for &u in &u_grid {
        let mut q_k = vec![0.0; nx * nx];
        let mut nx1 = vec![0usize; nx * nx];
        let mut nx2 = vec![0usize; nx * nx];
        for i2 in 0..nx {
            for i1 in 0..nx {
                // Take number of integration steps
                let (mut x1_k, mut x2_k, mut q) = (x1[i1], x2[i2], 0.0);
                for _ in 0..nk {
                    // RK4 integration for x1, x2 and q
                    let (k1_x1, k1_x2, k1_q) = f(x1_k, x2_k, u);
                    let (k2_x1, k2_x2, k2_q) = f(x1_k + dt / 2.0 * k1_x1, x2_k + dt / 2.0 * k1_x2, u);
                    let (k3_x1, k3_x2, k3_q) = f(x1_k + dt / 2.0 * k2_x1, x2_k + dt / 2.0 * k2_x2, u);
                    let (k4_x1, k4_x2, k4_q) = f(x1_k + dt * k3_x1, x2_k + dt * k3_x2, u);
                    x1_k += dt / 6.0 * (k1_x1 + 2.0 * k2_x1 + 2.0 * k3_x1 + k4_x1);
                    // NB: the Python original updates X1_k in place while k1_x2 aliases it,
                    // so its x2 update uses the new x1 (minimal cost 2.817428...); this is plain RK4
                    x2_k += dt / 6.0 * (k1_x2 + 2.0 * k2_x2 + 2.0 * k3_x2 + k4_x2);
                    q += dt / 6.0 * (k1_q + 2.0 * k2_q + 2.0 * k3_q + k4_q);
                }

                // Find out which state comes next (index)
                let mut j1 = ((x1_k + 1.0) / 2.0 * (nx - 1) as f64).round_ties_even() as i64;
                let mut j2 = ((x2_k + 1.0) / 2.0 * (nx - 1) as f64).round_ties_even() as i64;

                // Infinite cost if state gets out-of-bounds
                if j1 < 0 { q = INF; j1 = 0; }
                if j2 < 0 { q = INF; j2 = 0; }
                if j1 >= nx as i64 { q = INF; j1 = 0; }
                if j2 >= nx as i64 { q = INF; j2 = 0; }

                q_k[idx(i2, i1)] = q;
                nx1[idx(i2, i1)] = j1 as usize;
                nx2[idx(i2, i1)] = j2 as usize;
            }
        }

        // Save the stage cost and next state
        next_x1.push(nx1);
        next_x2.push(nx2);
        stage_j.push(q_k);
    }

    // Calculate cost-to-go (no end cost) and optimal control
    let mut j = vec![0.0; nx * nx];
    let mut u_opt_tab: Vec<Vec<usize>> = vec![];
    for _k in (0..n).rev() {
        // Cost to go for the previous step, optimal control action
        let mut j_prev = vec![INF; nx * nx];
        let mut u_prev = vec![usize::MAX; nx * nx];

        // Test all control actions
        for uind in 0..nu {
            for e in 0..nx * nx {
                let j_prev_test = j[idx(next_x2[uind][e], next_x1[uind][e])] + stage_j[uind][e];
                if j_prev_test < j_prev[e] {
                    u_prev[e] = uind;
                    j_prev[e] = j_prev_test;
                }
            }
        }

        // Update cost-to-go and save optimal control
        j = j_prev;
        u_opt_tab.push(u_prev);
    }

    // Reorder U_opt by stage
    u_opt_tab.reverse();

    // Find optimal control starting at x1=0, x2=1
    let mut i1 = nx / 2;
    let mut i2 = nx - 1;
    let mut u_opt = vec![];
    let mut x1_opt = vec![x1[i1]];
    let mut x2_opt = vec![x2[i2]];
    let mut cost = 0.0;
    for k in 0..n {
        // Get the optimal control and go to next step
        let u_ind = u_opt_tab[k][idx(i2, i1)];
        cost += stage_j[u_ind][idx(i2, i1)];
        (i1, i2) = (next_x1[u_ind][idx(i2, i1)], next_x2[u_ind][idx(i2, i1)]);

        // Save the trajectories
        u_opt.push(u_grid[u_ind]);
        x1_opt.push(x1[i1]);
        x2_opt.push(x2[i2]);
    }

    // Optimal cost
    println!("Minimal cost:  {}", cost);
    assert!((cost - j[idx(nx - 1, nx / 2)]).abs() < 1e-8); // Consistency check

    // Print the trajectories
    println!("x1 trajectory: {:?}", x1_opt);
    println!("x2 trajectory: {:?}", x2_opt);
    println!("u trajectory: {:?}", u_opt);
}
