

%% Create summary figure with main IRFs
send_endogenous_variables_to_workspace;

iniT=2; %First period for flows  graphs
maxy=10;
maxT=maxy/Deltat;            % Last year we want for the graphs

Y_dev=(Y(1:end)-Y(end))/Y(end)*100;
exTFP=(Q).^alpha ;

try
    TFP=(Q.*X./(1-Omega)).^alpha ;
    endTFP=(X./(1-Omega)).^alpha ;
catch
    endTFP=exTFP*0+1 ;
    TFP=exTFP ;
    Z=0*C+1;
end

TFP_dev=(TFP(1:end)-TFP(end))/TFP(end)*100;

switch modelfilename
case 'HANK'
    slope=(alpha*((1-alpha)./w).^((1-alpha)/alpha).*mbar.^(1/alpha)).*gamma./q.*Q;
case 'RANK'
    slope=w.*0+1;
end

Tn=N_period*Deltat;               % Time in years 
years=[Deltat:Deltat:Tn]-Deltat; % Vector of time in years

figure(16)
subplot(3,3,1),  hold on, plot(years(iniT:maxT),(Pi(iniT:maxT)-Pi(end))*100,'LineWidth',2)         , title('Pi')
subplot(3,3,2),  hold on, plot(years(iniT:maxT),(mbar(iniT:maxT)-mbar(end))/mbar(end)*100,'LineWidth',2)           , title('m')
subplot(3,3,3),  hold on, plot(years(iniT:maxT),(w(iniT:maxT)-w(end))/w(end)*100,'LineWidth',2), title('w')
subplot(3,3,4),  hold on, plot(years(iniT:maxT),(r(iniT:maxT)-r(end))*100,'LineWidth',2)             , title('r')
subplot(3,3,5),  hold on, plot(years(iniT:maxT),((q(iniT:maxT)-q(end))/q(end)*100),'LineWidth',2), title('q')
subplot(3,3,6),  hold on, plot(years(iniT:maxT),Y_dev(iniT:maxT),'LineWidth',2)         , title('Y')
subplot(3,3,7),  hold on, plot(years(iniT:maxT),TFP_dev(iniT:maxT),'LineWidth',2), title('TFP')
subplot(3,3,8),  hold on, plot(years(iniT:maxT),((zstar(iniT:maxT)-zstar(end))/zstar(end)*100),'LineWidth',2), title('zstar')
subplot(3,3,9),  hold on, plot(years(iniT:maxT),((slope(iniT:maxT)-slope(end))/slope(end)*100),'LineWidth',2), title('slope')


if I_save==1
    save([modelfilename,'_shocktype',num2str(shocktype),'_policy',num2str(OMP)])
end

