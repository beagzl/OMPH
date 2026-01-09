//------------------------------------------------------------------------
//
// HANK.mod
//
// January 2026
//
// Replication codes: 
// " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
//
// (c) Beatriz Gonzalez, Galo Nuño and Dominik Thaler
// 
// Dynare mod file solving for the HANK model with firm-level OU shocks
//
//------------------------------------------------------------------------


@#include "HANK_DeclareVars.mod"
// declare variables that are not in the automatically generated code 
var logY logZ r Yn;
@#if OMP==1
    @#include "HANK_DeclareMult.mod"
@#endif


varexo TPS MPS CPS; 

parameters Ibirth Kss rho rho_hh gamma alpha sigma theta Deltat  psi Upsilon  epsilon phiTR pibar thetaP eta_q share_d rhoTR phiK N_period I_zlb zlb IP zeta eta varpsi shocktype tau nocost delta;

@#for i in 1:N
    parameters z@{i} ;
    parameters varrho@{i} ;
    parameters xi@{i} ;
    parameters bettaA@{i} ;
    parameters Deltazc@{i} ;
    parameters f@{i} ;
@#endfor

@#for i in 1:N-1
    parameters Deltaz@{i} ;
@#endfor

model_parameters

HANK_normalization

rho      = params(1);
rho_hh   = params(2);
delta    = params(3);
gamma    = params(4);
alpha    = params(5);
psi      = params(6);
sigma    = params(7);
theta    = params(8);
Deltat   = params(9);
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
IP       = params(25);
zeta     = params(26);
eta      = params(28);
varpsi   = params(29);
SSsubsidy= params(30);
shocktype= params(31);
tau      = params(35);
nocost   = params(38);
Ibirth   = params(39);
Kss      = params(34);


@#for  i in 1:N
    z@{i}      = params(end-7*@{N}+@{i});
    varrho@{i} = params(end-6*@{N}+@{i});
    xi@{i}     = params(end-5*@{N}+@{i});
    bettaA@{i} = params(end-4*@{N}+@{i});
    Deltazc@{i}= params(end-3*@{N}+@{i});
    f@{i}      = params(end-2*@{N}+@{i});
@#endfor

@#for i in 1:N-1
    Deltaz@{i}= params(end-1*@{N}+@{i});
@#endfor

WS=0;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MODEL BLOCK
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

model;

@#include "HANK_DynamicsEquCond.mod"

//-------------------------------------------------
// Planner follows a Taylor rule or exogenous rule
//-------------------------------------------------

@#if OMP == 0
   i     = (MPS==0)*((1-rhoTR*Deltat)*i(-1) + rhoTR*Deltat*(rho_hh+phiTR*(Pi-pibar)+pibar)) + MPS; // INERTIA 
@#endif

//---------------------------------------
// Planner solves Ramsey problem
//---------------------------------------

@#if OMP==1
    @#include "HANK_DynamicsRamsey.mod"
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
logY=log(Y);
logZ=log(Z);
steady_state(m)/(1-tau)= (Upsilon/(1 - alpha)) * Yn^((psi + alpha)/(1 - alpha)) * Z^(-(psi + 1)/(1 - alpha)) * K^(-alpha*(psi + 1)/(1 - alpha)) * (Yn - delta*K)^zeta;

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

@#include "HANK_SSfilecall.mod"

logY        = log(Y);
logZ        = log(Z);
TRCdev       = 0;
Q           = 1;
Welf        = 0;
rho_hh_t    = rho_hh;
q           = 1;
R           = q*r+delta*q;            
G           = 0;
tau_t       = tau;
Yn          = Y;


//----------------------------------------
// SS multipliers of Ramsey problem
//----------------------------------------
@#if OMP == 1
    @#include "HANK_SSmultfilecall.mod"
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

    oo_.exo_simul([UBT],strmatch('MPS',M_.exo_names(:,:)))= rho_hh+pibar+shocksize./10000;
end

// Time preference shock
//----------------------------------------

if shocktype==2  //Time preference shock    
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
    oo_.endo_simul(strmatch('Q',M_.endo_names(:,:)),1)=start(strmatch('Q',M_.endo_names(:,:)))*(1+shocksize/10000);
end



//----------------------------------------
// LINEAR APPROXIMATION or NON-LINEAR
//----------------------------------------

if linearize==1 // Linear approximation
   perfect_foresight_solver(maxit=50,tolx=1e-9,tolf=1e-9,linear_approximation);
elseif linearize==0 // Non-linear solution
   perfect_foresight_solver(maxit=8,tolx=1e-10,tolf=1e-8);
elseif linearize==2 // Use linear approxiation as first guess for non-linear solution
   perfect_foresight_solver(maxit=50,tolx=1e-9,tolf=1e-9,linear_approximation);
   options_.linear_approximation = false;
   perfect_foresight_solver(maxit=8,tolx=1e-8,tolf=1e-8);
end


//----------------------------------------
// FIGURES
//----------------------------------------

postsimulation;



