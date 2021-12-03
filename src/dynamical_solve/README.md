
Workflow: 

  Dynamical_step: compute a step in internal variables (x,y,X,Y) and external parameter step using Eq(12) and (13)
  
     Input: central SDP, surrounding SDPs 
     
     Output: external_step in external parameter space
     
  run_dynamical: use Dynamical_step as a subroutine, run multiple iterations of Dynamical_step until terminating criteria are met 
  
```
├── Dynamical_Solver
│   ├── dynamical_navigator
│   │   ├── dynamical_step.cxx             # update (dx, dy, dX, dY) and external_step 
│   │   ├── mixed_hess.cxx                 # Compute H_xp and H^-1_xx H_xp 
│   │   ├── external_grad_hess.cxx         # Compute H_pp and Del_p L
│   │   ├── internal_search_direction.cxx  # Compute H^-1_xx Del_x L  (sign?)
│   │   ├── compute_update_sdp.cxx         # Decide whether to carry out the steps above (called by run_dynamical). TODO
│   │   ├── imports.hxx
│   │   ├── compute_lag.cxx
│   │   └── compute_search_direction.cxx
│   ├── run_dynamical.cxx                  # calls dynamical_step till termination conditions are satisfied
│   ├── Dynamical_Solver.cxx               # Resembles SDP_Solver
│   ├── load_checkpoint
│   │   ├── load_binary_checkpoint.cxx
│   │   ├── load_checkpoint.cxx
│   │   └── load_text_checkpoint.cxx
│   ├── run
│   │   ├── compute_feasible_and_termination.cxx
│   │   ├── print_header_dynamical.cxx
│   │   └── print_iteration.cxx
│   └── save_checkpoint.cxx
├── Dynamical_Solver.hxx
├── Dynamical_Solver_Parameters
│   ├── Dynamical_Solver_Parameters.cxx
│   ├── ostream.cxx
│   └── to_property_tree.cxx
├── Dynamical_Solver_Parameters.hxx
├── Dynamical_Solver_Terminate_Reason
│   └── ostream.cxx
└── Dynamical_Solver_Terminate_Reason.hxx

```
     
