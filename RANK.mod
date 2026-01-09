//------------------------------------------------------------------------
//
// RANK.mod
//
// January 2026
//
// Replication codes: 
// " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
//
// (c) Beatriz Gonzalez, Galo Nuño and Dominik Thaler
// 
// Dynare mod file solving for the RANK model 
//
//------------------------------------------------------------------------


@#include "RANK_DeclareVars.mod"

var r Yn;

@#if OMP==1
    @#include "RANK_DeclareMult.mod"
@#endif


varexo TPS MPS CPS; 

parameters phiK nocost rho rho_hh delta gamma alpha  Deltat Kss  psi Upsilon  epsilon phiTR pibar thetaP eta_q share_d rhoTR zeta N_period  zlb SSsubsidy tau Zhat wedge;

  
model_parameters

RANK_normalization

rho      = params(1);
rho_hh   = params(2);
delta    = params(3);
gamma    = params(4);
alpha    = params(5);
psi      = params(6);
sigma    = params(7);
theta    = params(8);
Deltat   = params(9);
Deltaz   = params(10);
Upsilon  = params(11);
thetaP   = params(12);
epsilon  = params(13);
phiTR    = params(14); 
pibar    = params(15);
eta_q    = params(16);
timeless = params(17);
share_d  = params(18);
rhoTR    = params(19);
phiK     = params(20);
N_period = params(22);
I_zlb    = params(23);
zlb      = params(24);
zeta     = params(26);
SSsubsidy= params(30);
tau      = params(35);
nocost   = params(36);

wedge    = 1; // Wedge to reproduce HANK path
Zhat     = 1; // Exogenous productivity in RANK = Productivity in HANK


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MODEL BLOCK
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


model;

@#include "RANK_DynamicsEquCond.mod"

//---------------------------------------
// Planner follows a Taylor rule
//---------------------------------------
@#if OMP == 0
    i     = (MPS==0)*((1-rhoTR*Deltat)*i(-1) + rhoTR*Deltat*(rho_hh+phiTR*(Pi-pibar)+pibar)) + MPS; // Taylor rule with inertia
@#endif

//---------------------------------------
// Planner solves Ramsey problem
//---------------------------------------

@#if OMP==1
     @#include "RANK_DynamicsRamsey.mod"
@#endif

//---------------------------------------
// Exogenous path of pi=0
//---------------------------------------

@#if OMP == 2
    Pi=steady_state(Pi);
@#endif

//---------------------------------------
// Define auxiliy variables (optional)
//---------------------------------------

r=((q(+1)-q)/Deltat + R -delta*q(@{If}))/q(@{If});
steady_state(m)/(1-tau)= (Upsilon/(1 - alpha)) * Yn^((psi + alpha)/(1 - alpha)) * (Zhat*Q^alpha)^(-(psi + 1)/(1 - alpha)) * K^(-alpha*(psi + 1)/(1 - alpha)) * (Yn - delta*K)^zeta;


end;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// STEADY STATE
//
// Provide SS of model, knowing that Pi=0 in Ramsey SS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


steady_state_model; // provide SS as function of Pi without optimal policy

@#if OMP != 1
    Pi       = pibar;
@#else
    Pi=0;
@#endif

//----------------------------------------------------------------
// Get the steady state conditional on Pi, normalizing q=1 
//----------------------------------------------------------------

[L,r,m,Kss,w,i,Y,T,C,rho_hh_t,Q,Welf,q,R,mbar]=RANK_SSfile(Pi,0);

tau_t=tau;
K = Kss;
Yn = Y;

//----------------------------------------
// SS multipliers of Ramsey problem
//----------------------------------------
@#if OMP == 1
    @#include "RANK_SSmultfilecall.mod"
@#endif
end;

//----------------------------------------
// Find SS of model with optimal policy
//----------------------------------------

steady(maxit=20,tolf=1e-9);
check;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TRANSITIONAL DYNAMICS
//
// Given the SS and the model equations, get transitional dynamics
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//----------------------------------------
// Preliminaries
//----------------------------------------

options_.periods = round(N_period);
start = oo_.steady_state;

//Solve with SS as initial guess
perfect_foresight_setup;
 
//---------------------------------------- 
// TIMELESS vs TIME 0
//----------------------------------------

if timeless==1
    //by default dynare does timeless
else
    oo_.endo_simul(strmatch('MULT',M_.endo_names(:,:)),1)=0;
end

//----------------------------------------
// DEFINE SHOCK
//----------------------------------------

// Monetary policy shock
//----------------------------------------


if shocktype==1 // Monetary policy shock
    oo_.exo_simul([2],strmatch('MPS',M_.exo_names(:,:)))= rho_hh+pibar+shocksize./10000;
end


// Time preference shock
//----------------------------------------

if shocktype==2 
    oo_.exo_simul([2], strmatch('TPS',M_.exo_names(:,:),'exact'))= shocksize/10000;
end


// Cost push shock
//----------------------------------------

if shocktype==3 
     oo_.exo_simul([2], strmatch('CPS',M_.exo_names(:,:),'exact'))= shocksize./10000;
end


// TFP shock
//----------------------------------------

if shocktype==4
    oo_.endo_simul(strmatch('Q',M_.endo_names(:,:)),1)=start(strmatch('Q',M_.endo_names(:,:)))*(1+shocksize/10000);  // TFP  shock
end



if linearize==1
   perfect_foresight_solver(maxit=10,tolx=1e-13,tolf=1e-13, linear_approximation, stack_solve_algo=0);
elseif linearize==0
   perfect_foresight_solver(maxit=10,tolx=1e-10,tolf=1e-10, stack_solve_algo=0);
elseif linearize==2
   perfect_foresight_solver(maxit=10,tolx=1e-13,tolf=1e-13, linear_approximation, stack_solve_algo=0);
   options_.linear_approximation = false;
   perfect_foresight_solver(maxit=10,tolx=1e-10,tolf=1e-10, stack_solve_algo=0);
end


//----------------------------------------
// FIGURES
//----------------------------------------

Igraphs_all=0; //Indicator for graph folder we are computing RANK

// Variables not available in RANK needed for graphs
clear X zstar D A 
zstar=Q*0+1;
D=Q*0+1;
A=D;

postsimulation;


