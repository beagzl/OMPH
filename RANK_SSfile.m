function [varargout] = RANK_SSfile(pi,firstcall)

%------------------------------------------------------------------------
%
% RANK steady state file
%
% January 2026
%
% Replication codes: 
% " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
%
% (c) Beatriz Gonzalez, Galo Nuno, Dominik Thaler
%
% This file solves for the SS in RANK for a given pi (input)
%
%     - firstcall = 0 if we are solving for initial SS where we normalize
%       upsilon st L=1, and find Kss such that q=1 
%
%     - firstcall = 1 if Upsilon is treated as a parameter,L is endogenous
%
% Includes function SS_given_L, residual and findzstar
%
%------------------------------------------------------------------------


%initiate params
global params Zhat wedge

pi
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
rhoTR    = params(19);
zeta     = params(26);
tau= params(35);

Zhat     = 1;
wedge    = 1;

%Initialize
K        = 0;
w        = 0;
Y        = 0;
T        = 0;
C        = 0;
I        = 0;

% -- STEADY STATE  --------------------

m        = (epsilon-1)/epsilon + rho_hh*pi*thetaP/epsilon ; % CT 
mbar     = m/(1-tau); % = m if no subsidy; =(m/((epsilon-1)/epsilon)) if subsidy
R        = (rho_hh+delta);
r        = (R-delta);
i        = r + pi;
q        = 1;
rho_hh_t = rho_hh;
Q        = 1;
Welf     = 0;


L        = 1;
SS_given_L;


if firstcall==0
    opt2 = optimoptions('fsolve','FunctionTolerance',1e-14,'OptimalityTolerance',1e-14,'Display','off');
    sol = fsolve(@residual,1,opt2);   
    residual(sol)
end


varargout{1}=L;
varargout{2}=r;
varargout{3}=m;
varargout{4}=K;
varargout{5}=w;
varargout{6}=i;
varargout{7}=Y;
varargout{8}=T;
varargout{9}=C;
varargout{10}=rho_hh_t;
varargout{11}=Q;
varargout{12}=Welf;
varargout{13}=q;
varargout{14}=R;
varargout{15}=mbar;

varargout{end+1}=L^(-psi)*w/(C^zeta);


if firstcall==1
    varargout{1}=cell2mat(varargout);
end

%--------------------------------------------------------------       
function [out] =SS_given_L
    K        = (R/(wedge*alpha*mbar*Zhat*L^(1-alpha)))^(1/(alpha-1));
    w        = (1-alpha)*Zhat*K^alpha*L^(-alpha)*(mbar);
    Y        = Zhat*K^alpha*L^(1-alpha);
    T        = (1-mbar)*Y-thetaP/2*pi^2*Y+(R/wedge-R)*K;
    C        = ((R-delta*q)*K+w*L+T);
    I        = delta*K;  
end
%--------------------------------------------------------------       

%--------------------------------------------------------------       
function out = residual(x_in) % to find L; calls the SS_give_L function
        L        = x_in;
        SS_given_L;
        out      = L - (Upsilon/(w)*C^zeta)^(1/-psi);
end
%--------------------------------------------------------------       



end








