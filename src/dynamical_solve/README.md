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
│   │   ├── mixed_hess.cxx                 # Compute H_xp and (- H^-1_xx H_xp) 
│   │   ├── external_grad_hess.cxx         # Compute H_pp and Del_p L
│   │   ├── internal_search_direction.cxx  # Compute (- H^-1_xx Del_x L)  
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
1. The external_step is scaled by dual_step_length at the very end of dynamical_step.cxx.
2. The hessians and gradients are printed out around line 205 in dynamical_step.cxx. Running on multiple cores will gives several copies of these lines and can make the log file rather messy. 
3. More detailed comments can be found within the files. 

# Updates
1. By choosing the a higher precision when approximating the OPE vector in the haskell input, we solved the resolution problems with the hessians and gradients. For now, we require that ```theta``` is given as a BigFloat of 512 bits and the OPE vector is approximated by a closeby rational with the resolution ```e^(-200)```

```
approxRational 1e-200 <$> toV (cos theta, sin theta)
```
With this setup, we found that choosing alpha step size anywhere from e^-20 to e^-70 gives good approximations to the derivatives. 

2. After correcting some sign errors, we found a primal_dual optimal solution at the minimum in the external parameter (theta in this example) after ~300 iterations. 
```
primalObjective = -0.00119767921900...
```
The BFGS algorithm that started from the same initial conditions took 782 iterations in total to get to 
```
primalObjective = -0.00119614134740... 
```
and 855 iterations in total to
```
primalObjective = -0.00119770961688...
```
(Whether the paramters I used was optimal for the BFGS algorithm is not clear to me, but I believe that the algorithm has hot-starting implemented) 

3. I'm now working on finding the boundary of the island rather than the minimum of the navigator function. These lines are commented out in the lastest commit. 



