%------------------------------------------------------------------------
%
% Model Parameters and Options Configuration
% Template file for setting up HANK/RANK model simulations
%
% This file should be customized for your specific needs and saved as
% 'model_parameters.m' in your working directory.
%
% Usage:
%   1. Set model type and options below
%   2. Run Buildfile.m to generate Dynare files
%   3. Run run_mod.m to execute simulation
%
% January 2026
%
% Replication codes: 
% " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
%
% (c) Beatriz Gonzalez, Galo Nuno, Dominik Thaler
% 
% This file sets up the parameters of the model
%
%------------------------------------------------------------------------

global params

%% ========================================================================
% POLICY REGIME and SHOCK CONFIGURATION
% ========================================================================

OMP =1;                     % 0= Taylor rule; 
                            % 1= OMP; 
                            % 2= Exogenous policy Pi=0; 

timeless  = 1;              % 1=timeless; 0=time 0  [only if OMP=1]

shocktype = 2;              % 0=No shock;
                            % 1=MP shock (transitory); 
                            % 2=Time preference shock (transitory)
                            % 3=Cost Push shock (transitory)
                            % 4=TFP shock (transitory)

shocksize = -100;           % Additive shock (positive) in basis points 
UBT = 2;                    % Period when shock hits >=2 (2 means annoucment=implementation period)


%% ========================================================================
% SOLUTION METHOD
% ========================================================================

linearize = 0;             % Solution method:
                           % 0 - Nonlinear perfect foresight
                           % 1 - Linear approximation
                           % 2 - Linear first, then nonlinear (more robust to initial guess)


%% ========================================================================
% SIMULATION HORIZON
% ========================================================================

Deltat = 1/12;             % Time step in years
                           % 1/12 = monthly (typical)
                           % 1/4 = quarterly
                          
N_period  = round(200/Deltat); % Number of periods to have 200 years 


%% ========================================================================
%  MODEL PARAMETERS
% ========================================================================

rho_hh   = .01;               %  Rate of time preference of HH
zeta     = 1;                 %  Inverse of IES (risk aversion coefficient) 1=log 
psi      = 1;                 %  Inverse Frisch elasticity

varpsi   = 0.1;               % Size at entry
eta      = 0.1;               % Exogenous exit shock
rho      = eta*(1-varpsi);    % aux parameter for KF: net dividend rate

gamma    = 1.56;              %  Borrowing constraint parameter

alpha    = 0.35;              %  Capital share in production
delta    = 0.06;              %  Capital depreciation rate

L        = 1;                 %  Normalization of labor supply Ls=1
Upsilon  = 1;                 %  Labor weight - endogenized st L=1  - this value is just an initialization
Kss      = 1;                 % Fixed capital, first guess

thetaP   = 100;               %  Price adjustment cost
epsilon  = 10;                %  Markup
pibar    = 0;                 %  SS inflation

tau      = 1/epsilon;         % SS subsidy  

phiTR    = 1.5;               % Taylor rule coefficient 
rhoTR    = 0.2;               % Inertia in TR: quarterly DT = 0.95, CT (1-0.95)/(1/4)=0.2
eta_q    = 0.2;               % Aggregate shock persistence, yearly DT = 0.8, CT (1-0.8)/1=0.2

zlb      = 0;                 % Zero-lower bound level

vaux= zeta + (psi + alpha) / (1 - alpha); % Aux parameter for FOCS



%% ========================================================================
% TIMING INDICES (for Dynare)
% ========================================================================
% These control timing conventions in the discretized model
% Usually don't need to change these

Is=0; %timing of state vars
If=0; %timing of forward vars
Ic=0; %timing of contemporaneous vars

%% ========================================================================
% HETEROGENEOUS PRODUCTIVITY
% ========================================================================

% Creating productivity grid
%-------------------------------------------------

N    = 200;               % number of grid points of productivity z
n    = 100;               % number of points concentrated around zstar(ss)
zmin = 0.001;             % min grid z
zmax = 48;                % max grid z


differencetype=1;% 1 UNEVEN grid, central
                 % 2 EVEN grid, backward
                 % 3 EVEN grid, forward
                 % 4 EVEN grid, central
                 % 5 EVEN grid, upwind                 

% Creates general grid in z                 
if differencetype>1
    z    = linspace(zmin,zmax,N)'; % Equally spaced grids
else
    z    =  exp(linspace(log(zmin),log(zmax),N-n))'; %log spacing
end

% Include fine grid that covers movement of zstar --> reduce approximation error
% Can be set if we know the range, or run twice and endogenously find it 
% Changes with shock type and shock size

if exist('second_iter')
    % If we set a loop to find [min_zstar, max_zstar] endogenously
    % Need to set min_zstar<min(zstar) and max_zstar>max(zstar) but
    % sufficiently close
        z=sort([z;linspace(min_zstar, max_zstar,n)']); 
else
    % If it is a known shock and we know the range, just set it
    
    if timeless==0          % timeless response to no shock
        min_zstar= 4.11;
        max_zstar= 5;
    elseif  shocktype==1    % MP
        min_zstar= 4.147;
        max_zstar= 4.18;
    elseif  shocktype==2    % TP shock for OMP
        min_zstar= 4.0599;
        max_zstar= 4.1500;
    elseif  shocktype==3    % CP shock for OMP
       min_zstar= 4.12;
       max_zstar= 4.155;
    elseif shocktype==4     % TFP shock for OMP 
        min_zstar= 4.14;
        max_zstar= 4.165;
    else                    % Default range
        min_zstar= 4.09;
        max_zstar= 4.15;
    end

    z=sort([z;linspace(min_zstar, max_zstar,n)']); % add extra points in the relevant range for zstar, where the "relevant range" depends on the shock and can not be knon ex ante

end

dz=z(2:end)-z(1:end-1);% Distance between points
dzf   = [dz;dz(end)];  % Distance between points forward   
dzb   = [dz(1); dz];   % Distance between points backward      
dzc=(dzf+dzb)/2;       % Distance between points average  
dz2=dzc.^2;
d2z=z(3:end)-z(1:end-2);% Distance between points j-1 and j+1 
d2zf=[d2z;2*(z(end)-z(end-1));2*(z(end)-z(end-1))];% Distance between points j and j+2
d2zb=[2*(z(2)-z(1));2*(z(2)-z(1));d2z;];% Distance between points j-2 and jJ
d2z=[(z(2)-z(1))*2;d2z;(z(end)-z(end-1))*2];
dzb2  = dzb.^2;
dzf2  = dzf.^2;

%% ========================================================================
% IDIOSYNCRATIC PRODUCTIVITY PROCESS
% ========================================================================
%%% -- OU PROCESS   --------------------
correl  = 0.8274653;               % persistence (at corr_dt frequency)
corr_dt = 1;                       % Frequency of corr --> quarterly
theta   = -log(correl)/corr_dt;    % Drift in logs: theta-->0: persisitence-->1; theta-->infty: persistence-->0
sig     = 0.7261;                  % Std dev of the shock of the process logs 
sig2    = sig^2;

% Using Ito's Lemma, transform the OU process of log z into levels
% dz = (-theta_z*log(z) + sig2/2) z dt + z^2 sig2 dW
mu    = (-theta*log(z) + sig2/2).*z;      % drift of process in levels 
s2    = z.^2.*sig2;                       % variance of process in levels
s1    = s2.^0.5;                          % std dev of process in levels


% Preliminaries for B matrix --> Constant for all iterations
varrho=zeros(N,1);xi=zeros(N,1);beta_base=zeros(N,1);

if differencetype==1
    %central 1st derivative for UNEVEN grid
    beta_base(:) = -(s2(:)./(dzb.*dzf));
    varrho(:)    =  mu./[d2z(2:end);d2z(end)]   +   s2(:)./(dzf.*d2zf);
    xi(:)        = -mu./[d2z(1);d2z(1:end-1)]   +   s2(:)./(dzb.*d2zb);  
elseif differencetype==2
    % backward 1st derivative, EVEN grid
    varrho(:)    =  mu(:)./dzc+s2(:)./(2*dz2);
    beta_base(:) = -mu(:)./dzc-s2(:)./(dz2);
    xi(:)        =             s2(:)./(2*dz2);
elseif differencetype==3
    % forward 1st derivative, EVEN grid
    varrho(:)    =             s2(:)./(2*dz2);
    xi(:)        = -mu(:)./dzc+s2(:)./(2*dz2);
    beta_base(:) = +mu(:)./dzc-s2(:)./(dz2);
elseif differencetype==4    
    %central  1st derivative, EVEN grid
    varrho(:)    = +mu(:)./(2*dzc)+s2(:)./(2*dz2);
    xi(:)        = -mu(:)./(2*dzc)+s2(:)./(2*dz2);
    beta_base(:) =                -s2(:)./(dz2);
elseif differencetype==5    
    % upwind 1st derivative, EVEN grid
    ind=(mu>0);
    varrho(:)    =  mu(:)./dzc.*ind                            +s2(:)./(2*dz2);
    xi(:)        =                     -mu(:)./dzc.*(1-ind)    +s2(:)./(2*dz2);
    beta_base(:) = -mu(:)./dzc.*ind    +mu(:)./dzc.*(1-ind)    -s2(:)./(dz2);
end

fz=lognpdf(z,(sig^2/(2*theta))^.5);

%% ========================================================================
% OTHER AUXILIARY INDICATORS
% ========================================================================
auxil1    = 1;              % 0=turns off the endogeniety of Z ; 1=baseline
I_zlb     = 0;              % 0=no ZLB; 1=ZLB at zlb level (set later)
I_save    = 1;              % Save mat file at end of simulation

%%

params     =[[rho rho_hh delta gamma alpha psi sig theta Deltat nan Upsilon thetaP epsilon phiTR pibar eta_q timeless rho rhoTR 0 N N_period I_zlb zlb 1 zeta differencetype eta varpsi nan shocktype linearize 1 Kss tau auxil1 vaux 1 1]';z; varrho; xi; beta_base; dzc;fz; dz; nan];