function metrics = build_metrics(Opt, info, pinfo, prec_time, total_time)

eff = info.efficiency;

metrics = struct();

metrics.method = "Davidson";
metrics.total_time = total_time;
metrics.preconditioner_setup_time = prec_time;

metrics.initialization_time = pinfo.total_time + info.initial_subspace_time;
metrics.davidson_solver_time = max(0, ...
    eff.total_time - info.initial_subspace_time);

metrics.outer_iterations = info.outer_iterations;
metrics.restart_count = info.restart_count;

metrics.final_R_max = info.final_R_max;
metrics.converged = logical(info.converged);
metrics.success = logical(info.success);

metrics.fallback_used = logical(info.fallback_used);
metrics.inner_failure_count = info.inner_failure_count;

metrics.solid_shifted_solve_calls = eff.linear.solid.calls;
metrics.fluid_shifted_solve_calls = eff.linear.fluid.calls;
metrics.total_shifted_solve_calls = eff.total_shifted_solve_calls;

metrics.solid_krylov_iterations = eff.linear.solid.iter_total;
metrics.fluid_krylov_iterations = eff.linear.fluid.iter_total;
metrics.total_krylov_iterations = eff.total_krylov_iterations;

metrics.solid_operator_applications = eff.matvec.solid.count;
metrics.fluid_operator_applications = eff.matvec.fluid.count;
metrics.total_operator_applications = eff.total_operator_applications;

metrics.solid_matvec_time = eff.matvec.solid.time;
metrics.fluid_matvec_time = eff.matvec.fluid.time;
metrics.total_matvec_time = eff.total_matvec_time;

metrics.total_linear_solve_time = eff.total_linear_solve_time;

metrics.mass_norm_error = info.mass_norm_error;
metrics.qep_orth_error = info.qep_orth_error;

metrics.init_solid = string(Opt.init.solid);
metrics.init_fluid = string(Opt.init.fluid);
metrics.preconditioner_enabled = logical(Opt.preconditioner.enable);

metrics.target_tolerance = Opt.eig_tol;
metrics.inner_tolerance = Opt.linear.tol;

end