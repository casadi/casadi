/** \class CasADi::Sundials::KinsolInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::KinsolSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunctionInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>symbolic_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>generate jacobian symbolically by source code transformation</td><td>CasADi::SXFunctionInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunction
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>symbolic_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>generate jacobian symbolically by source code transformation</td><td>CasADi::SXFunctionInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration (gives an initial value for INTEGRATOR_T0, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration (gives an initial value for INTEGRATOR_TF, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesIntegrator
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration (gives an initial value for INTEGRATOR_T0, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration (gives an initial value for INTEGRATOR_TF, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
</table>
*/
/** \class CasADi::JacobianInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>finite_differences</td><td>OT_BOOLEAN</td><td>false</td><td>Using finite differences instead of automatic differentiation</td><td>CasADi::JacobianInternal</td></tr>
</table>
*/
/** \class CasADi::Jacobian
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>finite_differences</td><td>OT_BOOLEAN</td><td>false</td><td>Using finite differences instead of automatic differentiation</td><td>CasADi::JacobianInternal</td></tr>
</table>
*/
/** \class CasADi::CplexInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
</table>
*/
/** \class CasADi::CplexSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoInterface
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoInternal</td></tr>
</table>
*/
/** \class CasADi::IntegratorInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration (gives an initial value for INTEGRATOR_T0, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration (gives an initial value for INTEGRATOR_TF, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
</table>
*/
/** \class CasADi::Integrator
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration (gives an initial value for INTEGRATOR_T0, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration (gives an initial value for INTEGRATOR_TF, which will be removed)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLUInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLU
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
</table>
*/
/** \class CasADi::ParallelizerInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
</table>
*/
/** \class CasADi::Parallelizer
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
</table>
*/
/** \class CasADi::IntegratorJacobianInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>derivative_index</td><td>OT_INTEGER</td><td>INTEGRATOR_P</td><td>integrate with respect to what?</td><td>CasADi::IntegratorJacobianInternal</td></tr>
</table>
*/
/** \class CasADi::IntegratorJacobian
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>derivative_index</td><td>OT_INTEGER</td><td>INTEGRATOR_P</td><td>integrate with respect to what?</td><td>CasADi::IntegratorJacobianInternal</td></tr>
</table>
*/
