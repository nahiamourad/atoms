 PROGRAM Atoms_and_perturbation
 implicit none

 integer::Ng,Na,Ls,lrho,ms,Zi,Zf,wVxc,logic_oda,grid_div
 integer::norm,Rayleih,my_method,basis,continuation,pert
 integer::diff_NI,diff_R,diff_beta,threshold
 integer::ninter,ninter0,step_NI,final_NI,exceed,max_iter
 double precision::scal0,steplength
 double precision::epsE,eps,vc,Length,Length0,beta
 double precision::step_R,final_R,step_beta,final_beta,final_eps,final_epsE
!%%%%%%%%%%%%%%%%%%%%%%%
!%%integration parameter
 Ng=15           !%Ng: number of Gauss points for the radial variable
 Na=30           !%Na: number of Gauss points for the \theta variable
!%%%%%%%%%%%%%%%%%%%%%%%
!%%charge and spherical harmonics
 Zi=25            !%Zi:initial atomic number
 Zf=25            !%Zf:final atomic number
 wVxc=0          !<------ logic 0 or 1: 0 means rHF and 1 means X\alpha
 Ls=2            !%Ls: maximum angular momentum
 ms=2            !%ms:maximum magnetic number
!%%%%%%%%%%%%%%%%%%%%%%%
!%%grid
 scal0=0.9d0     !%scale0: 0<scale0\le 1, the scale used to build the grid (h_{k-1} = scal0 * h_k), if equal 1 then uniform grid is used
 ninter=90       !%ninter: number of intervals
 Length=200.d0    !%Length: length of the simulation interval
 grid_div=1      !<------ logic 0 or 1: 1 means that the radial grid is divided into two parts,the first one depends on scal0 and the second one is uniform 
                 !its length is is "Length0+(Length-Length0)" and has "ninter0+(ninter-ninter0)" intervals
 ninter0=80       !%if grid_dev=1, ninterv0 is the number of intervals for the first part of the grid (the non-uniform one)
 Length0=80.d0    !%if grid_dev=1, Length0 is the Radius of the first part of the grid (the non-uniform one)
 diff_NI=0       !<------ logic 0 or 1: 1 means run for different number of grid points
 step_NI=0       !%if diff_NI=1, step_NI is the difference of two consecutive number of intervals (can be positive or negative)
 final_NI=0      !%if diff_NI=1, final_NI is the maximum (or minimum) desired number of intervals
 diff_R=0        !<------ logic 0 or 1: 1 means run for different grid Length
 step_R=0.d0     !%if diff_R=1, step_R is the difference of two consecutive Radii (can be positive or negative)
 final_R=0.d0    !%if diff_R=1, final_R is the largest (or smallest) desired Radius
!%%%%%%%%%%%%%%%%%%%%%%%
!%%threshold
 epsE=1.d-10     !%epsE: is the eps used for the convergence of the energy
 eps=1.d-7       !%eps: is the eps used for the convergence of the H1_norm of the eigenfunctions if turned on
 vc=1.d-4        !%vc: parameter used for the difference bt the eigenvalues in the accidental degeneracy case
 threshold=0     !<------ logic 0 or 1: 1 means run for different eps and epsE
 final_eps=0.d0  !%if threshold=1, final_eps is the smallest desired epsilon, where the initial is eps and the step is *10{-1}
 final_epsE=0.d0 !%if threshold=1, final_epsE is the smallest desired epsilon, where the initial is epsE and the step is *10{-1}
!%%%%%%%%%%%%%%%%%%%%%%%
!%%electric field
 beta=0.d-2      !%beta: \beta
 diff_beta=0     !<------ logic 0 or 1: 1 means run for different beta
 step_beta=0.d0  !%if diff_beta=1, step_beta is the difference between two consecutive betas (can be positive or negative)
 final_beta=0.d0 !%if diff_beta=1, final_beta is the maximum (or. minimum) desired beta
!%%%%%%%%%%%%%%%%%%%%%%%
!%%ODA
 logic_oda=1    !<------ logic 0 or 1: 1 means ODA is used
 steplength=0.7d0!%steplength: is the step length used if no ODA (logic_oda=0)  
!%%%%%%%%%%%%%%%%%%%%%%%
!%%optional
 basis=2        !<------- 0,1 or 2: 1 means Robin boundary condition turned on, 2 means the base 1/r outside the ball of simulation is added,0 else wise
                !N.B. 1 and 2 are defined just for the rHF case 
 continuation=0 !<------ logic 0 or 1: 1 means the continuation principle is used i.e the initial point is the convergent point of the previous calculations
 norm=0         !<------ logic 0 or 1: 1 means the H1 norm of the eigenfunctions is calculated
 Rayleih=0      !<------ logic 0 or 1: 1 means that the converged eigenvalues are calculated by Rayleigh quotient
 my_method=0    !<------ logic 0 or 1: 1 means that the occupation number is calculated by minimizing the energy at each iteration (special for z=5 or 6)
 pert=0         !<------ logic 0 or 1: 1 means skip all after calculating the matrices and go to perturbation.f95 file
 lrho=0!2*Ls      !%lrho: is the maximum l such that rho_l is not zero.... 
!%%%%%%%%%%%%%%%%%%%%%%%
!%force to stop
exceed=1        !<------ logic 0 or 1: 1 means if the number of iteration exceeds "max_iter" then stop
max_iter=500    !%if exceed=1, max_iter is the maximum number of iteration desired
!%%%%%%%%%%%%%%%%%%%%%%%
 call main(Zi,Zf,wVxc,beta,Ng,Na,Ls,ms,lrho,basis,steplength,logic_oda,grid_div,scal0,Length,Length0,ninter,ninter0,&
           epsE,eps,vc,continuation,norm,Rayleih,my_method,pert,diff_NI,step_NI,final_NI,&
           diff_R,step_R,final_R,diff_beta,step_beta,final_beta,threshold,final_eps,final_epsE,exceed,max_iter)
End program Atoms_and_perturbation
