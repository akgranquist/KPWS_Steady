The Project_summary.tex file is a LaTeX document with a summary of all of the analytical and numerical results. 
The MATLAB files integrate the Boussinesq equation using a spectral IFRK4 method. They are as follows:
  Boussinesq_RK4: This integrates the Boussinesq equation 
        u_tt - u_xx + a(u^2)_xx + g u_xxxx = 0, adapted to the system: 

        u_t = v_x
        v_t = -gam*u_xxx + u_x - alpha*(u^2)_xx

        on [-L,L] with initial conditions u(x,0) and v(x,0).
        This method is not used in the main numerical analysis, because it cannot integrate Riemann boundary conditions. It is included for completenes' sake, but not recommended for use.
        
  Boussinesq_RK4_Riemann: 
      Integrates Boussinesq equation
      u_tt - u_xx + a(u^2)_xx + g u_xxxx = 0, adapted to the system: 
      
      w_t = v_x
      v_t = w_x - gam*w_xxx - 2*alpha*(w^2 + antideriv(w)*w_x

      where w = u_x and v_x = w_t (ie v = u_t). 
      This method is used in the main numerical analysis. It should be able to integrate boundary conditions that 
        a) decay at the edges of the domain, e.g. a soliton; or
        b) approach constant values at the edges of the domain, e.g. Riemann boundary conditions.

  Boussinesq_comparison_right_moving:
       Integrates Riemann + soliton boundary conditions. This evolves to the interaction between a soliton and a rarefaction wave. 
       The conditions are set such that the soliton starts to the right of the jump and moves to the right. The solution is compared at each time step to the analytical solution as described in the Project Summary.

  Boussinesq_comparison_left_moving:
       Integrates Riemann + soliton boundary conditions. This evolves to the interaction between a soliton and a rarefaction wave. 
       The conditions are set such that the soliton starts to the right of the jump and moves to the left. The solution is compared at each time step to the analytical solution as described in the Project Summary.


