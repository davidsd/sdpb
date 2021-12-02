
Workflow: 
  Dynamical_step: compute a step in internal variables (x,y,X,Y) and external parameter step using Eq(12) and (13)
     Input: central SDP, surrounding SDPs 
     Output: external_step in external parameter space
  run_dynamical: use Dynamical_step as a subroutine, run multiple iterations of Dynamical_step until terminating criteria are met 
  
     
     
