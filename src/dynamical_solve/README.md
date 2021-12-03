The structure of dynamical_solve is shown below:   
```
src/dynamical_solve
├── Dynamical_Solver
│   ├── Dynamical_Solver.cxx               # Resembles SDP_Solver
│   │ 
│   ├── run_dynamical.cxx                  # calls dynamical_step till termination conditions are satisfied
│   │ 
│   ├── dynamical_navigator
│   │   ├── dynamical_step.cxx             # update (dx, dy, dX, dY) and external_step 
│   │   ├── mixed_hess.cxx                 # Compute H_xp and H^-1_xx H_xp 
│   │   ├── external_grad_hess.cxx         # Compute H_pp and Del_p L
│   │   ├── internal_search_direction.cxx  # Compute H^-1_xx Del_x L  (sign?)
│   │   ├── compute_update_sdp.cxx         # Decide whether to carry out the steps above (called by run_dynamical). TODO
│   │   ├── imports.hxx                    # Functions imported from sdp_solve 
│   │   ├── compute_lag.cxx                # Currently not used. Call Approx_Objective instead 
│   │   └── compute_search_direction.cxx.  # Currently not used.
│   │ 
│   ├── load_checkpoint
│   │   ├── load_binary_checkpoint.cxx
│   │   ├── load_checkpoint.cxx
│   │   └── load_text_checkpoint.cxx
│   ├── run
│   │   ├── compute_feasible_and_termination.cxx
│   │   ├── print_header_dynamical.cxx
│   │   └── print_iteration.cxx
│   └── save_checkpoint.cxx
│   │ 
├── Dynamical_Solver.hxx
├── Dynamical_Solver_Parameters
│   ├── Dynamical_Solver_Parameters.cxx
│   ├── ostream.cxx
│   └── to_property_tree.cxx
│   │ 
├── Dynamical_Solver_Parameters.hxx
├── Dynamical_Solver_Terminate_Reason
│   └── ostream.cxx
│   │ 
└── Dynamical_Solver_Terminate_Reason.hxx

```

# Notes
1. There's potentially a **sign error** with internal_search_direction, as I'm not sure if it computes `- H^-1_xx Del_x L` or `+ H^-1_xx Del_x L`.
2. The external_step is scaled by dual_step_length at the very end of dynamical_step.cxx.
3. More detailed comments can be found within the files. 
