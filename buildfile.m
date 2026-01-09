%------------------------------------------------------------------------
%
% HANK/RANK Model Builder for Dynare
%
% This script generates Dynare code for HANK (Heterogeneous Agent New Keynesian)
% and RANK (Representative Agent New Keynesian) models with Ramsey optimal policy.
%
% Output files:
%   - [model]_SSmultfile.m: Steady-state multiplier calculation
%   - [model]_DynamicsEquCond.mod: Conditional equilibrium dynamics
%   - [model]_DynamicsRamsey.mod: Ramsey optimal policy dynamics
%   - [model]_DeclareVars.mod: Variable declarations
%   - [model]_DeclareMult.mod: Multiplier declarations
%
% To include a new model, add a new case. User needs to:
%   1. Declare variables according to their type: Backward-looking,
%   forward-looking or static variables
%   2. Include the full set of equlibirum condition:
%       - funcODE: Dynamic equations (ODEs)
%       - funcSTAT: Static equilibrium conditions
%       - Obj: Objective function (welfare)
%   3. Call for the appropriate steady state functions
%   4. Run the buildfile to create the appropriate mod files, and create a
%   mod file that is consistent with this model and calls the created sub
%   mod files.
% 
% January 2026
%
% Replication codes: 
% " Firm Heterogeneity, Capital Misallocation and Optimal Monetary Policy"
%
% (c) Beatriz Gonzalez, Galo Nuno, Dominik Thaler
%
%------------------------------------------------------------------------


%% Initialize
model_parameters % load parameters and keep only those you need
clearvars -except modelfilename N
clc; %close all

% Note: there should be only one time derivative per equation and it should appear as a linear ly additive term. Its Ok to have more than one timederivative in one equation only if the sum of those terms can be written as a single time derivative
    
switch modelfilename
case 'HANK'
    % declare vars, pi cannot be used as varname
    % you cant use variable names starting in v_ or dot_
    omega = sym('omega', [N 1]);
    dot_omega = sym('dot_omega', [N 1]); %dot variable needs to be declared manually only for vectors %DONEST WORK
    syms   Q  K Kss A D C q Pi Y Omega phi w X Z zstar L TRCdev i m T R mbar rho_hh_t tau_t  

    % declare var types
    statevars=[ Q  A D tau_t rho_hh_t transpose(omega) ];
    forwardvars=[C q Pi Y];
    contvars=[K Omega phi w X Z zstar L TRCdev i m T R mbar ];


    %declare params and shock innovations
    syms rho rho_hh delta gamma alpha  psi Upsilon  epsilon  thetaP eta_q phiK Zhat wedge SSsubsidy zeta eta varpsi tau Deltat K 
    z           = sym('z', [N 1]);
    Deltazc     = sym('Deltazc', [N 1]);
    Deltaz      = sym('Deltaz', [N 1]);
    zstarzind   = sym('zstarzind', [N 1]); %zstarzind=(z>=zstar)
    returnind   = sym('returnind', [N 1]); %returnind=((z*Q*phi-R)>0)
    bettaA      = sym('bettaA', [N 1]);
    xi          = sym('xi', [N 1]);
    varrho      = sym('varrho', [N 1]);
    f           = sym('f', [N 1]);
    syms TPS MPS CPS CPS EXRSHOCK

    %generate time drivatives of variables automatically
    dot_statevars=[];
    for varname=statevars
        eval(['syms dot_', char(varname),';'])
        dot_statevars=[ dot_statevars eval(['dot_',char(varname)])];
    end
    dot_forwardvars=[];
    for varname=forwardvars
        eval(['syms dot_', char(varname),';'])
        dot_forwardvars=[ dot_forwardvars eval(['dot_',char(varname)])];
    end
    dot_contvars=[]; %these varaibles are generated for convenience but never appear in the output
    for varname=contvars
        eval(['syms dot_', char(varname),';'])
        dot_contvars=[ dot_contvars eval(['dot_',char(varname)])];
    end

    %%% dynamic equations 
    % declare state equations first, then the forward equations
    % states and forward vars must not be mixed in one equation!!!
    % spaces in expressions are somtimes causing issues, operators need to have a space on both sides or on no side!
    funcODE=[
    -dot_tau_t     + tau_t*(-eta_q)    + eta_q*tau     + CPS/Deltat
    -dot_Q  +  Q*(-eta_q)*log(Q)
    -dot_rho_hh_t  + rho_hh_t*(-eta_q) + eta_q*rho_hh + TPS/Deltat;
    -dot_A  +  1/q*((gamma*(1-Omega))*(mbar*alpha*Z*K^(alpha-1)*L^(1-alpha)-R)+R-delta*q-q*rho)*A
    -dot_D  +  ((R-delta*q)*D+w*L-C+T)/q
    -dot_omega(1)     + (1/q*(gamma*(z(1)*Q*phi-R)*returnind(1)        + R-q*delta-q*eta)-dot_A/A)*omega(1)    + eta*varpsi*omega(1)     +  bettaA(1)*omega(1)          + xi(1)*omega(1) + xi(2)*omega(2) 
    -dot_omega(N)     + (1/q*(gamma*(z(N)*Q*phi-R)*returnind(N)        + R-q*delta-q*eta)-dot_A/A)*omega(N)    + eta*varpsi*omega(N)     +  bettaA(N)*omega(N)          + varrho(N-1)*omega(N-1) +  varrho(N)*omega(N)                                        
    -dot_omega(2:N-1) + (1/q*(gamma*(z(2:N-1)*Q*phi-R).*returnind(2:N-1)+R-q*delta-q*eta)-dot_A/A).*omega(2:N-1)+eta*varpsi*omega(2:N-1) +  bettaA(2:N-1).*omega(2:N-1) + varrho(1:N-2).*omega(1:N-2) + xi(3:N).*omega(3:N)
    -i      +  (R-delta*q+dot_q)/q+Pi
    -dot_q*C^(-zeta)+zeta*q*C^(-zeta-1)*dot_C  +  rho_hh_t*q*C^(-zeta)  -  q*C^(-zeta)*(R-delta*q)/q
    -epsilon/thetaP*(m-(epsilon-1)/epsilon)*C^-zeta*Y  +  rho_hh_t*C^-zeta*Y*Pi  -  dot_Pi*C^-zeta*Y+dot_C*zeta*C^(-zeta-1)*Y*Pi-dot_Y*C^-zeta*Pi
    ];

    mstate=5+N; % Nr of backward looking euqations
    mfor=3; % Nr of forward looking euqations


    %%% static equations 

    funcSTAT=[
    -phi    + alpha*((1-alpha)/w)^((1-alpha)/(alpha))*mbar^(1/alpha)
    -w      + (1-alpha)*Z*K^alpha*L^(-alpha)*mbar
    -R      + alpha*Z*K^(alpha-1)*L^(1-alpha)*zstar/(X/(1-(Omega)))*mbar
    -K      + D+A
     Kss    - K
    -D      + A*(gamma*(1-Omega)-1)
    -Z*(1-Omega)^alpha + (Q*X)^alpha
    -Y      + Z*K^alpha*L^(1-alpha)
    -T      + ((1-mbar)*Y-thetaP/2*Pi^2*Y+q*rho*A)+(delta*q-delta)*K % Profits of retailers + price adj cost + dividends + profits cap good producer
    -TRCdev + Y-C-thetaP/2*Pi^2*Y-delta*K
    -X      + sum(z.*omega.*Deltazc.*zstarzind)  + sum(((omega(2:N).*z(2:N)+0.5*(z(2:N)-zstar)./Deltaz(1:N-1).*(omega(1:N-1).*z(1:N-1)-omega(2:N).*z(2:N))).*(z(2:N)-zstar)./Deltaz(1:N-1)-omega(2:N).*z(2:N)/2).*zstarzind(2:N).*(1-zstarzind(1:N-1)).*Deltaz(1:N-1))  
    -Omega  + sum(omega.*Deltazc.*(1-zstarzind))+sum(((omega(1:N-1)+0.5.*(zstar-z(1:N-1))./Deltaz(1:N-1).*(omega(2:N)-omega(1:N-1))).*(zstar-z(1:N-1))./Deltaz(1:N-1)-omega(1:N-1)./2).*zstarzind(2:N).*(1-zstarzind(1:N-1)).*Deltaz(1:N-1)) 
    -mbar   + m/(1-tau_t)
    -w      + Upsilon*L^psi*C^zeta
    ];

    %Objective
    Obj= C^(1-zeta)/(1-zeta)-Upsilon*(L^(1+psi))/(1+psi);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
case 'RANK'
    % declare vars, pi cannot be used as varname
    % you cant use variable names starting in v_ or dot_
    syms   Q C q Pi Y w  L i m R mbar cons_t rho_hh_t tau_t K
    syms omega % not used in RANK
    
    % declar var types
    statevars=[ K Q tau_t rho_hh_t];
    forwardvars=[C q  Pi Y];
    contvars=[ w L i m R mbar];

    %declare params and shock innovations
    syms rho rho_hh delta gamma alpha  psi Upsilon  epsilon  thetaP eta_q phiK Zhat wedge SSsubsidy zeta eta varpsi  tau Deltat K Kss
    syms TPS MPS CPS EXRSHOCK CPS


    %generate dot vairables automatically
    dot_statevars=[];
    for varname=statevars
        eval(['syms dot_', char(varname),';'])
        dot_statevars=[ dot_statevars eval(['dot_',char(varname)])];
    end
    dot_forwardvars=[];
    for varname=forwardvars
        eval(['syms dot_', char(varname),';'])
        dot_forwardvars=[ dot_forwardvars eval(['dot_',char(varname)])];
    end
    dot_contvars=[]; %these varaibles are generated for convenience but never appear in the output
    for varname=contvars
        eval(['syms dot_', char(varname),';'])
        dot_contvars=[ dot_contvars eval(['dot_',char(varname)])];
    end


    %%% dynamic equations 
    % declare state equations first, then the forward equations
    % dot_variables must appear linearly (could be relaxed)!!! 
    % dot_variables must only be multiplied by other dynamic variables, not static ones (could be relaxed)!!!
    % states and forward vars must not be mixed in one equation!!!
    % spaces in expressions are sometimes causing issues, operators need to have a space on both sides or on no side!
    funcODE=[
    -dot_tau_t     + tau_t*(-eta_q)    + eta_q*tau     + CPS/Deltat
    -dot_Q  +  Q*(-eta_q)*log(Q)
    -dot_rho_hh_t  + rho_hh_t*(-eta_q) + eta_q*rho_hh + TPS/Deltat;
    -i      +  (R-delta*q+dot_q)/q+Pi
    -dot_q*C^(-zeta)+zeta*q*C^(-zeta-1)*dot_C  +  rho_hh_t*q*C^(-zeta)  -  q*C^(-zeta)*(R-delta*q)/q
    -epsilon/thetaP*(m-(epsilon-1)/epsilon)*C^-zeta*Y  +  rho_hh_t*C^-zeta*Y*Pi  -  dot_Pi*C^-zeta*Y+dot_C*zeta*C^(-zeta-1)*Y*Pi-dot_Y*C^-zeta*Pi

    ];

    mstate=3; % Nr of backward looking euqations
    mfor=3; % Nr of forward looking euqations


    %%% static equations 

    funcSTAT=[
    -w      + (1-alpha)*Zhat*Q^alpha*K^alpha*L^(-alpha)*mbar  
    -R      + alpha*Zhat*Q^alpha*K^(alpha-1)*L^(1-alpha)*mbar*wedge
    Kss     - K
    -Y      + Zhat*(Q*K)^alpha*L^(1-alpha)
    -C+Y-(delta*K + thetaP*1/2*Pi^2*Y)
    -mbar   + m/(1-tau_t)
    -w      + Upsilon*L^psi*C^zeta
    ];

    %Objective
    Obj= C^(1-zeta)/(1-zeta)-Upsilon*(L^(1+psi))/(1+psi);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate Ramsey FOCs
nstate=length(statevars);
nfor=length(forwardvars);
ncont=length(contvars);
ndyn=nstate+nfor;

mdyn=length(funcODE);
mstat=length(funcSTAT);

if mstate+mfor~=mdyn, error('Nr ODEs doesnt match sum of Nr of forward and backward ODEs'),end
if -mdyn-mstat+ndyn+ncont<1, error('No degree of freedom for the planner'),end

display(['Private equilibrium has ',num2str(mdyn+mstat),' equations and ',num2str(ndyn+ncont),' variables, so ',num2str(-mdyn-mstat+ndyn+ncont),' degree(s) of freedom for the planner'])

MULT=[];dot_MULT=[];
for ii=1:length(funcSTAT)+length(funcODE)
    eval(['syms MULT',num2str(ii),';']);
    eval(['syms dot_MULT',num2str(ii),';']);
    eval(['MULT=[MULT MULT',num2str(ii),'];']);
    eval(['dot_MULT=[dot_MULT dot_MULT',num2str(ii),'];']);
end

% construct Hamiltonian
H=Obj+MULT*[funcODE ; funcSTAT];

%find FOCs in CT
JHx_dot=jacobian(H,[dot_statevars,dot_forwardvars]); %derivative wrt dot_variables to get the term premultiplying the dot_variable
JHx_dot_dot=jacobian(JHx_dot,[MULT,statevars,forwardvars])*transpose([dot_MULT,dot_statevars,dot_forwardvars]); %time derivative of the above term
FOCdyn=transpose(jacobian(H,[statevars forwardvars]))+rho_hh_t*transpose(JHx_dot)-JHx_dot_dot;
FOCstate=FOCdyn(1:nstate);
FOCforward=FOCdyn(nstate+[1:nfor]);
FOCcont=transpose(jacobian(H,contvars));

funcODE=[funcODE(1:nstate);FOCforward;funcODE(nstate+1:end);FOCstate]
funcSTAT=[funcSTAT;FOCcont]
m2dyn=length(funcODE);
m2stat=length(funcSTAT);

%% find SS multipliers conditional on knowing the Ramsey SS inflation rate
% In case one doesnt know the Ramsey SS inflation rate, one would have to find it numerically.
FOCsm1=subs([FOCforward;FOCstate;FOCcont(1:end)],[dot_forwardvars dot_statevars dot_MULT],0*[dot_forwardvars dot_statevars dot_MULT]);
FOCsmat=jacobian(FOCsm1,[MULT]); % coefficients in front of multipliers
FOCscon=subs(FOCsm1,MULT,MULT*0);% constant

FOCsconc=char(FOCscon);
FOCsmatc=char(FOCsmat);
output_args=char(sort(symvar(MULT)));
input_args=char(symvar([FOCsmat,FOCscon]));
omega_args=char(sort(symvar(omega))); % vector of
vars_args=char(symvar([statevars, forwardvars ,contvars]));
if verLessThan('matlab','9.13') %syntax depends on matlab version
    omega_args=omega_args(10:end-3);
    output_args=output_args(10:end-3);
    input_args=input_args(10:end-3);
    vars_args=vars_args(10:end-3);
    FOCsconc=['FOCscon=',FOCsconc(8:end-1),';'];
    FOCsmatc=[FOCsmatc(8:end-1)];
else
    omega_args=omega_args(2:end-1);
    output_args=output_args(2:end-1);
    input_args=input_args(2:end-1);
    vars_args=vars_args(2:end-1);
    FOCsconc=['FOCscon=',FOCsconc,''';'];
end
for ii=N:-1:1
    FOCsmatc=strrep(FOCsmatc,['returnind',num2str(ii)],['((z',num2str(ii),'*Q*phi-R)>0)']);
    FOCsmatc=strrep(FOCsmatc,['zstarzind',num2str(ii)],['(z',num2str(ii),'>=zstar)']);
    FOCsconc=strrep(FOCsconc,['returnind',num2str(ii)],['((z',num2str(ii),'*Q*phi-R)>0)']);
    FOCsconc=strrep(FOCsconc,['zstarzind',num2str(ii)],['(z',num2str(ii),'>=zstar)']);
    input_args=strrep(input_args,[', zstarzind',num2str(ii)],'');
    input_args=strrep(input_args,[', returnind',num2str(ii)],'');
end
FOCsmatc=strrep(FOCsmatc,'], [','] \n [');
FOCsmatc=['FOCsmat=',FOCsmatc,';'];

header1=['[',output_args,' ] = ',modelfilename,'_SSmultfile( ',input_args,' )'];
header=['function ',header1];
invertMULT='MULT=FOCsmat\\-FOCscon'';';
splitMULT='for i=1:length(MULT) \n eval([''MULT'',num2str(i),''=MULT(i);'']) \n end \n';
fid = fopen([modelfilename,'_SSmultfile.m'],'wt');
fprintf(fid,[header,'\n',FOCsconc,'\n',FOCsmatc,'\n',invertMULT,'\n',splitMULT,'\n','end']); %print SS multiplier file
fclose(fid);% 
fid = fopen([modelfilename,'_SSmultfilecall.mod'],'wt'); 
fprintf(fid,[header1,';']); %print line that calls SS multiplier file
fclose(fid);
switch modelfilename
case  {'HANK'}
    fid = fopen([modelfilename,'_SSfilecall.mod'],'wt');
    text=['[Omega, phi, w, r, A, X, K, D, Z, C, zstar, L, m, i, Y, T, Aux2,mbar,' omega_args,' ]=HANK_SSfile(Pi,0);'];
    fprintf(fid,text); %print line that calls manually geerated SS file
    fclose(fid);
end
fid = fopen([modelfilename,'_DeclareVars.mod'],'wt'); 
text=['var ',vars_args,';'];
fprintf(fid,text); %print line that calls manually geerated SS file
fclose(fid);
fid = fopen([modelfilename,'_DeclareMult.mod'],'wt'); 
text=['var ',output_args,';'];
fprintf(fid,text); %print line that calls manually geerated SS file
fclose(fid);


%% generate discretized dynare code
funcODE2=funcODE;
funcSTAT2=funcSTAT;
clearvars system
for variable=[forwardvars statevars contvars MULT] %rename variabels at t by vvariable
    variablec=char(variable);
    eval(['syms v_',variablec]);
    funcODE2=subs(funcODE2,variable,eval(['v_',variablec]));
    funcSTAT2=subs(funcSTAT2,variable,eval(['v_',variablec]));
end

iii=0;clearvars variablec %generate character vectors of all variables
for variable=[forwardvars statevars contvars MULT dot_forwardvars dot_statevars dot_MULT]
    iii=iii+1;
    variablec{iii}=char(variable);
end
forwardvarsec=variablec(1:nfor);
statevarsec=variablec(nfor+[1:nstate]);
contvarsec=variablec(ndyn+[1:ncont]);
MULTec=variablec(ndyn+ncont+[1:length(MULT)]);
dot_forwardvarsec=variablec(ndyn+ncont+length(MULT)+[1:nfor]);
dot_statevarsec=variablec(ndyn+ncont+length(MULT)+nfor+[1:nstate]);
dot_MULTec=variablec((ndyn+ncont+length(MULT)+nfor+nstate+1):end);

for iii=1:(m2dyn+m2stat)
    if iii<=m2dyn %dyn equs
        RHS=char(funcODE2(iii));
    else
        RHS=char(funcSTAT2(iii-m2dyn));
    end
    %replace indicators
    for ii=N:-1:1
        RHS=strrep(RHS,['returnind',num2str(ii)],['((z',num2str(ii),'*Q*phi-R)>0)']);
        RHS=strrep(RHS,['zstarzind',num2str(ii)],['(z',num2str(ii),'>=zstar)']);
    end

    %timing of cont
    for variable=flip([contvarsec MULTec(mstate+mfor+1:end)])
        RHS=strrep(RHS,['v_',variable{1}],variable{1});
    end
    
    %timing of states
    for variable=flip([statevarsec MULTec(mstate+[1:mfor])])
        RHS=strrep(RHS,['v_',variable{1}],[variable{1},'(@{Is})']);
    end
    
        %timing of forwards
    for variable=flip([forwardvarsec MULTec(1:mstate)])
        RHS=strrep(RHS,['v_',variable{1}],[variable{1},'(@{If})']);
    end
    

    for variable=flip([dot_statevarsec dot_MULTec(mstate+[1:mfor])])
        variablec2=strrep(variable{1},'dot_','');
        RHS=strrep(RHS,variable{1},['(',variablec2,'-',variablec2,'(-1)',')/Deltat']);
    end
    
    for variable=flip([dot_forwardvarsec dot_MULTec(1:mstate)])
        variablec2=strrep(variable{1},'dot_','');
        RHS=strrep(RHS,variable{1},['(',variablec2,'(+1)-',variablec2,')/Deltat']);
    end
    system{iii,1}=strcat('0=',RHS,';');
end



fid = fopen([modelfilename,'_DynamicsEquCond.mod'],'wt');
fid = fopen([modelfilename,'_DynamicsEquCond.mod'],'at');
for iii=[1:mstate  mstate+nfor+[1:mfor]  mstate+nfor+mfor+nstate+[1:mstat] ]
       fprintf(fid,['\n', system{iii}]);
end
fclose(fid);


fid = fopen([modelfilename,'_DynamicsRamsey.mod'],'wt');
fid = fopen([modelfilename,'_DynamicsRamsey.mod'],'at');
for iii=[mstate+[1:nfor] mstate+nfor+mfor+[1:nstate] mstate+nfor+mfor+nstate+mstat+[1:ncont]]
       fprintf(fid,['\n', system{iii}]);
end
fclose(fid);

