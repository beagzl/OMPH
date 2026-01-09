function [varargout] = HANK_SSfile(pi,firstcall)

%------------------------------------------------------------------------
%
% HANK steady state file
%
% January 2026
%
% Replication codes: 
% " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
%
% (c) Beatriz Gonzalez, Galo Nuno, Dominik Thaler
%
% This file solves for the SS in HANK for a given pi (input)
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
global params

pi
rho       = params(1);
rho_hh    = params(2);
delta     = params(3);
gamma     = params(4);
alpha     = params(5);
psi       = params(6);
sig       = params(7);
theta     = params(8);
Deltat    = params(9);
%dz        = params(10);
Upsilon   = params(11);
thetaP    = params(12);
epsilon   = params(13);
phiTR     = params(14);
pibar     = params(15);
eta_q     = params(16);
timeless  = params(17);
share_d   = params(18);
rhoTR     = params(19);
phiK      = params(20);
J         = params(21);
IP        = params(23);
I_zlb     = params(24);
zlb       = params(25);
zeta      = params(26);
differencetype= params(27);
eta      = params(28);
varpsi   = params(29);
SSsubsidy= params(30);
Kss      = params(34);
tau      = params(35);
nocost   = params(38);
Ibirth   = params(39);


z         = params(end-7*J+1:end-6*J);
varrho    = params(end-6*J+1:end-5*J);
xi        = params(end-5*J+1:end-4*J);
beta_base = params(end-4*J+1:end-3*J);
dzc       = params(end-3*J+1:end-2*J);
fz        = params(end-2*J+1:end-1*J);
dz        = params(end-1*J+1:end-1);

dz2       = dz.^2;
sig2  = sig.^2;       %  variance shock process in logs

% -- PRELIMINARIES  --------------------
maxit = 10000;
tol   = 1e-11; %increased from 1e-11 on 19/03/2024

opt     = optimoptions('fsolve','FunctionTolerance',1e-13,'OptimalityTolerance',1e-13,'Display','off');

% Preallocation
c=zeros(J,1);s=zeros(J,1);shat=zeros(J,1);
omega=zeros(J,1);
f0=zeros(J,1);cdf0=zeros(J,1);f_d=zeros(J,1);
cdf=zeros(J,1) ;onecdf=zeros(J,1);
B=zeros(J,J);
beta_add=zeros(J,1);
y=zeros(J,1); Xaux_vec=zeros(J,1); X=0;Xaux=0;xstar=0;varphi=0;Zaux=0;
A=0;D=0;w=0;C=0;Y=0;T=0;

% Initialize guesses
L       = 1;                       % Normalization of labor supply Ls=1
f0(:)   = lognpdf(z,0,sqrt(sig2)); % Guess for initial distribution
f0      = f0 ./sum(f0.*dzc);
fnew    = f0;

% % Boundary conditions
   beta_base(1) = beta_base(1)+xi(1);
   beta_base(J) = beta_base(J)+varrho(J);

% Construct part of B Matrix that is invariant
B_base = spdiags(beta_base(:),0,J,J)+spdiags(varrho(:),-1,J,J)+spdiags([0;xi(2:J)],1,J,J);


% -- STEADY STATE  --------------------

r       = rho_hh;
i       = rho_hh + pi;
m       = (epsilon-1)/epsilon + rho_hh*pi*thetaP/epsilon;
mbar    = m/(1-tau); %=m if no subsidy, =m/((epsilon-1)/epsilon) if subsidy
R       = (rho_hh+delta);

zstar=.7;

SS_given_L;


if firstcall==0
    sol=fsolve(@residual,1,opt);
    residual(sol);
end


% CHECKS
% This should be 1
sumationerror=sum(f_d.*dzc);
% This should be 0
sum(shat.*f_d.*dzc);


Z      = ((X)/(1-Omega))^alpha;
Omega  = -xstar/gamma+1;
phi    = varphi;
X      = Xaux;
i      = r+pi;


varargout{1}    = Omega;
varargout{2}    = phi;
varargout{3}    = w;
varargout{4}    = r;
varargout{5}    = A;
varargout{6}    = X;
varargout{7}    = Kss;
varargout{8}    = D;
varargout{9}    = Z;
varargout{10}   = C;
varargout{11}   = zstar;
varargout{12}   = L;
varargout{13}   = m;
varargout{14}   = i;
varargout{15}   = Y;
varargout{16}   = T;
varargout{17}   = sumationerror;
varargout{18}   = mbar;
for i=1:J
    varargout{end+1} = f_d(i);
end
varargout{end+1}=L^(-psi)*w/(C^zeta); % Upsilon:parameter multiplying disutility of labor

if firstcall==1 & nargout ==1
    varargout{1} = cell2mat(varargout);
end
if firstcall==1 & nargout ==19
    varargout{19} = f_d;
end
%Bb;
%----------------------------------------------------------------
    function [out] =SS_given_L % to find distribution given L;
        
        for n=1:maxit
            
            f_d  = fnew;
            
            % 1- Compute cdf
            if n==1;zstar=.7;end
            findzstar2   = @(x) findzstar(x,z,f_d,dz,dzc,rho_hh,delta,gamma,rho,alpha,epsilon,1);
            [zstar fval exitflag]        = fsolve(findzstar2,zstar,opt);
            Aux          =  findzstar(zstar,z,f_d,dz,dzc,rho_hh,delta,gamma,rho,alpha,epsilon,2);
            Omega        = Aux(1);
            X            = Aux(2);
            
            % 3- Define aux variable X, xstar and Z
            Xaux      = X;
            xstar     = (1-Omega)*gamma;
            Zaux      = (gamma*Xaux)^alpha;
            
            % 4- Get stock A, K and D
            A       = (((rho_hh+delta))/(alpha*mbar*Zaux*L^(1-alpha)*zstar/(gamma*X)))^(1/(alpha-1));
            Kss       = xstar*A;
            D       = Kss-A;
            
            % 5- Get wages, interest rates and the aux variable varphi
            w      = (1-alpha)*mbar*Zaux*A^alpha*L^(-alpha);
            varphi = alpha*((1-alpha)/w)^((1-alpha)/alpha)*mbar^(1/alpha);
            
            % 5.5- Get consumption
            Y      = Zaux*A^alpha*L^(1-alpha);
            T      = ((1-mbar)*Y) - thetaP/2*pi^2*Y + share_d*A;
            C      = w*L+(R-delta)*D+T;
            
            % 6- Get savings
            s(:)    = (gamma*max(z(:)*varphi-R,0)+R-delta-eta*(1-varpsi*Ibirth));
            shat(:) = s(:); %  shat=(s-Adot/A), in SS Adot=0
            
            % 7- Construct B matrix
            beta_add(:) = shat(:);
            
            
            % B matrix:
            Bb = spdiags(beta_add(:),0,J,J)+B_base;
            
            % To check it is well computed: sum(B,1) needs to be = to shat(:)

            Check=full(Bb); sum_check=sum(Check,1); B_error=max(abs(sum_check(:)-shat(:)));
            if max(abs(sum(dzc'*full(B_base),1)'))>1e-7 | max(abs(sum(dzc'*full(Bb),1)-dzc'.*shat'))>1e-7
                disp('-- B Matrix not well computed ---')
                max(max(abs(sum(dzc'*full(B_base),1)')),                max(abs(sum(dzc'*full(Bb),1)-dzc'.*shat')))
                % error('matrix')
            end
            
            B          = Bb;
            Baux       = (speye(J)-Deltat*B);
            ff    = spdiags(1+Deltat*shat(:),0,J,J)*f_d+fz*varpsi*eta*(1-Ibirth)*Deltat;  %savings explicit
            ff    = (speye(J)-Deltat*B_base)\ff;        %brownian motion implicit
            f_sum = dzc'*(ff);
          %  fnew  = ff./f_sum;
            fnew  = ff./1;
                   
            
            % Convergenge criteria
            fchange = (fnew - f_d)./f_d;
            dist(n) = max(max(abs(fchange)),max(abs(fchange.*f_d)));
            dist(n);
            n;
                               
            if dist(n)<tol
                if pi==0
                    disp('Wealth shares Converged, Iteration , Max Error ')
                    n
                    max(abs(-B_base*ff-(spdiags(shat(:),0,J,J)*ff)))
                end
                f_d=fnew;
                break
            end
            
            if n==maxit
                disp('-- Wealth shares did NOT converge --')
                pi
                dist(n)
            end
            
            
        end
        
    end
%----------------------------------------------------------------

%----------------------------------------------------------------
    function out=residual(x_in) % to find L; calls the SS_give_L function
        L=x_in;
        SS_given_L;
        out=L- (Upsilon/w*C^zeta)^(1/-psi);

    end
%----------------------------------------------------------------

%----------------------------------------------------------------
    function [ resid ] = findzstar(zstar,z,f_d,dz,dzc,rho_hh,delta,gamma,rho,alpha,epsilon,ind)
        
        loc_star  = find(z<zstar,1,'last');
        
        if (zstar<z(1));loc_star=1;end %extrapolate
        
        if (zstar>z(end));loc_star=size(f_d,1)-1;end %extrapolate
        
        Omega=sum(f_d(z<zstar).*dzc(z<zstar))...
            +((f_d(loc_star)+0.5*(zstar-z(loc_star))/dz(loc_star)*(f_d(loc_star+1)-f_d(loc_star)))*(zstar-z(loc_star))/dz(loc_star)-f_d(loc_star)/2)*dz(loc_star);
        
        X=sum(f_d(z>=zstar).*z(z>=zstar).*dzc(z>=zstar))...
            +((f_d(loc_star+1)*z(loc_star+1)+0.5*(z(loc_star+1)-zstar)/dz(loc_star)*(f_d(loc_star)*z(loc_star)-f_d(loc_star+1)*z(loc_star+1)))*(z(loc_star+1)-zstar)/dz(loc_star)-f_d(loc_star+1)*z(loc_star+1)/2)*dz(loc_star);
       
        if ind==1
            resid=((rho_hh+delta)*gamma*(1-Omega)+(rho-rho_hh)).*zstar./(gamma*X)-(rho_hh+delta);
        else
            resid=[Omega,X];
        end
    end
%----------------------------------------------------------------


end


