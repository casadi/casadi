# Test conclusions

## Tests

All QPs and all tests without box constraints on x were removed from the test
collection. That left a total of 219 test problems.

## Experiments

Four main experiments were carried out: 
1. ALM+PANOC with standard L-BFGS directions
2. ALM+PANOC with specialized L-BFGS directions without flushing when γ changes
3. ALM+LBFGS++ with all box constraints moved to the outer ALM problem
4. ALM+LBFGSB++

## Results

|                 | ALM+PANOC | ALM+PANOC (+) | ALM+LBFGS++ | ALM+LBFGSB++ |
|:----------------|:---------:|:-------------:|:-----------:|:------------:|
| Converged       | 138 / 219 | 143 / 219     | 116 / 219   | 150 / 219    |
| Total time      | 1860.50   | 1880.36       | 916.25      | 1135.14      |
| Conv. time      | 164.39    | 227.10        | 214.27      | 136.72       |
| All conv. time* | 31.73     | 51.75         | 11.55       | 6.93         |
| Objective evaluations* | 409942 | 373863 | 43169 | 21544 |
| Gradient evaluations* | 363342 | 316104 | 42492 | 20891 |
| Outer iterations* | 770 | 712 | 677 | 653 |
| Inner iterations* | 160001 | 130252 | 33657 | 17583 |
| Linesearch failures* | 4287 | 5979 | - | - |
| L-BFGS failures* | 0 | 0 | - | - |
| L-BFGS rejected* | 3109 | 2094 | - | - |

(*) Only for tests that converged for all four solvers.

The bottom three statistics are PANOC-specific, so they have no value for the
LBFGS++ columns.

An L-BFGS failure is a non-finite L-BFGS step (infinite or NaN). An L-BFGS 
rejection is a failed curvature or C-BFGS condition, resulting in the update
being rejected. 
