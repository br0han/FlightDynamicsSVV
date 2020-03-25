%close all;
rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
p0     = 101325;          % pressure sea level [pascal]
gamma  = 1.4;             
S      = 30.00;           %(m^2)
c      = 2.0569;	      % mean aerodynamic cord [m]
D      = 0.686;             % random engine diameter [m]


CmTc = -0.0064;          % thrust control derivative
Ws = 60500/g;      % standard aircraft weight [kg]

Scm_delta = -1.0738;
Fcm_delta = -1.4070;

Scm_alpha = -0.5145;
Fcm_alpha = -0.6405;


%%%%Data
SDat = load("Sample_Elevator.mat"); %Sample Data excluding mass in kg
SMass = load("SampleMasses.mat"); 
SPMass=SMass.PeopleMasses; %Sample Mass of the people in kg

FDat=load("Flight_Set2.mat"); %Measurements from flight
FMass=load("FlightMassData"); %Mass data from flight
FPMass=FMass.FlightData;


EmptyMass=9165; %(lbs)
EmptyMass=EmptyMass*0.453592; %(kg)
PreFuel=4100; %(lbs)
PreFuel=PreFuel*0.453592; %(kg)


STotPeopleMass=sum(SPMass); %(kg)
FTotPeopleMass=sum(FPMass);

STotPreMass=PreFuel+EmptyMass +STotPeopleMass;
FTotPreMass=PreFuel+EmptyMass +FTotPeopleMass;

%SAMPLE DATA
SAlt=SDat.SampleData1(1:7,1);    %Sample Altitude               (ft)
SVelo=SDat.SampleData1(1:7,2) - 2;   %Sample Velocity               (kts)
SAlpha=SDat.SampleData1(1:7,3);  %Sample Angle of Attack        (degrees)
SDelta=SDat.SampleData1(1:7,4);  %Sample Elevator Deflection    (degrees) 
SFe = SDat.SampleData1(1:7,6);   %Sample Elevator Stick Force   (N) 
SFFl=SDat.SampleData1(1:7,7);    %Sample Fuel flow left         (lbs/hr)
SFFr=SDat.SampleData1(1:7,8);    %Sample Fuel flow right        (lbs/hr)
SFUsed=SDat.SampleData1(1:7,9);  %Sample Fuel used              (lbs)
STemp=SDat.SampleData1(1:7,10);  %Sample Temperature            (Celsius)

SAlt= 0.3048*SAlt ;     %Sample Altitude            (m)
SVelo= 0.514444*SVelo;  %Sample Velocity            (m/s)
SAlpha= pi*SAlpha/180;  %Sample Angle of Attack     (radians)
SDelta = pi*SDelta/180  %Sample Elevator deflection (radians)
SFFl= 0.000125998*SFFl; %Sample Fuel flow left      (kg/s)
SFFr= 0.000125998*SFFr; %Sample Fuel flow right     (kg/s)
SFUsed=0.453592*SFUsed; %Sample Fuel used           (kg)
STemp=273.15+STemp ;    %Sample Temperature         (Kelvin)  

%FLIGHT DATA
FAlt=FDat.Data(1:7,1);   %Flight Altitude               (m)
FVelo=FDat.Data(1:7,2) - 2;  %Flight Velocity               (m/s)
FAlpha=FDat.Data(1:7,3); %Flight Angle of Attack        (degrees)
FDelta=FDat.Data(1:7,4); %Flight Elevator Deflection    (degrees) 
FFe = FDat.Data(1:7,6);  %FLight Elevator Stick Force   (N)
FFFl=FDat.Data(1:7,7);   %Flight Fuel Flow left         (lbs/hr)
FFFr=FDat.Data(1:7,8);   %Flight Fuel flow right        (lbs/hr)
FFUsed=FDat.Data(1:7,9); %Flight Fuel Used              (lbs)
FTemp=FDat.Data(1:7,10); %Flight Temperature            (Celsius)

FAlt=0.3048*FAlt ;       %Flight Altitude             (m)
FVelo=0.514444*FVelo ;   %Flight Velocity             (m/s)
FAlpha=pi*FAlpha/180 ;   %Flight AoA                  (radians)
FDelta = pi*FDelta/180 ; %%Flight Elevator Deflection (radians)
FFFl=0.000125998*FFFl ;  %Flight Fuel Flow Left       (kg/s)
FFFr=0.000125998*FFFr ;  %Flight Fuel Flow right      (kg/s)
FFUsed=0.453592*FFUsed ; %Flight Fuel used            (kg)
FTemp=FTemp+273.15 ;     %Flight Temp                 (Kelvin)


STempISA=[276.15; 275.58; 275.18; 274.53; 275.96; 276.65; 277.64];
FTempISA=[260.66; 260.42; 260.40; 260.07; 261.81; 262.60; 263.55]; 

sz=size(SAlt); 

% SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data
% FCL=zeros(sz(1),sz(2)); %Matrix containing CL for Flight data
% for i= 1:sz(1)
%    SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVelo(i,1)^2*S);
%    FCL(i,1)=2* (FTotPreMass-FFUsed(i,1))*g/(Frho(i,1)*FVelo(i,1)^2*S);
% end


SDeltaTemp = STemp-STempISA;
FDeltaTemp = FTemp-FTempISA;


% equivelent airspeed calculations 

Fp = ones(sz(1),sz(2));
Sp = zeros(sz(1),sz(2));

% pressure at different altitudes 
for i= 1:sz(1)
    Sp(i,1) = p0*(1+(lambda*SAlt(i,1)/Temp0))^(-g/(lambda*R));
    Fp(i,1) = p0*(1+(lambda*FAlt(i,1)/Temp0))^(-g/(lambda*R));
end

%Mach number at different altitudes
FM = zeros(sz(1), sz(2)); 
SM = zeros(sz(1), sz(2));
    
for i= 1:sz(1)
    SM(i,1) = sqrt(2/(gamma-1) *((1+ p0/Sp(i,1) *((1+ (gamma-1)*rho0*(SVelo(i,1)^2)/(2*gamma *p0))^(gamma/(gamma-1))-1))^((gamma -1)/gamma) -1));
    FM(i,1) = sqrt(2/(gamma-1) *((1+ p0/Fp(i,1) *((1+ (gamma-1)*rho0*(FVelo(i,1)^2)/(2*gamma *p0))^(gamma/(gamma-1))-1))^((gamma -1)/gamma) -1));
end 


% correcting the measured temperature. 
FT = zeros(sz(1), sz(2));
ST = zeros(sz(1),sz(2));

for i= 1:sz(1)
    FT(i,1) = FTemp(i,1)/(1+(gamma - 1)/2 *FM(i,1)^2);
    ST(i,1) = STemp(i,1)/(1+(gamma - 1)/2 *SM(i,1)^2);
end

% calculating the density using ideal gas law
Srho= zeros(sz(1),sz(2));
Frho= zeros(sz(1),sz(2));   

for i= 1:sz(1)%used to be Srho1 and Frho2
    Srho(i,1) = Sp(i,1)/(R*ST(i,1));
    Frho(i,1) = Fp(i,1)/(R*FT(i,1)); 
end

% calculating the true airspeed and equivelant airspeed 
FVeq = zeros(sz(1),sz(2));
SVeq = zeros(sz(1),sz(2));

FVet = zeros(sz(1),sz(2));
SVet = zeros(sz(1),sz(2));

for i = 1:sz(1)
    FVet(i,1) = sqrt(gamma*R*FT(i,1))*FM(i,1);
    SVet(i,1) = sqrt(gamma*R*ST(i,1))*SM(i,1);
    
    FVeq(i,1) = FVet(i,1)*sqrt(Frho(i,1)/rho0); % true airspeed = M*a , Veq = Vt*sqrt(rho/rho0)
    SVeq(i,1) = SVet(i,1)*sqrt(Srho(i,1)/rho0);
end


%Actual Thrust
SampleThrustParams = [SAlt SM SDeltaTemp SFFl SFFr]; 
FlightThrustParams = [FAlt FM FDeltaTemp FFFl FFFr];

save TrimCurve_Thrust\STP SampleThrustParams;
save TrimCurve_Thrust\FTP FlightThrustParams;

run('TrimCurve_Thrust\Thruster');

S_ActualThrustCoef = (S_Thrust(:, 1) + S_Thrust(:, 2))./(0.5*(D^2)*Srho.*(SVet.^2));
F_ActualThrustCoef = (F_Thrust(:, 1) + F_Thrust(:, 2))./(0.5*(D^2)*Frho.*(FVet.^2));


%Reduced thrust
%Standard Fuel flow
FFs = 0.048*ones(7,1);

SampleThrustParams = [SAlt SM SDeltaTemp FFs FFs]; 
FlightThrustParams = [FAlt FM FDeltaTemp FFs FFs];

save TrimCurve_Thrust\STP SampleThrustParams;
save TrimCurve_Thrust\FTP FlightThrustParams;

run('TrimCurve_Thrust\Thruster.m');

SVeqr = SVeq.*sqrt(Ws./(STotPreMass - SFUsed));
FVeqr = FVeq.*sqrt(Ws./(FTotPreMass - FFUsed));


S_ReducedThrustCoef = (S_Thrust(:, 1) + S_Thrust(:, 2))./(0.5*(D^2)*Srho.*(SVet.^2));
F_ReducedThrustCoef = (F_Thrust(:, 1) + F_Thrust(:, 2))./(0.5*(D^2)*Frho.*(FVet.^2));


%trim curve
SDelta_eq = SDelta - (CmTc/Scm_delta)*(S_ReducedThrustCoef - S_ActualThrustCoef);
FDelta_eq = FDelta - (CmTc/Fcm_delta)*(F_ReducedThrustCoef - F_ActualThrustCoef);




%force curve
SFe_eq = SFe.*(Ws./(STotPreMass - SFUsed));
FFe_eq = FFe.*(Ws./(FTotPreMass - FFUsed));


FVeqrplot = sort(FVeqr);
FDelta_eqplot = sort(FDelta_eq);%(sorting);
FFe_eqplot = sort(FFe_eq);
SVeqrplot = sort(SVeqr);
SDelta_eqplot = sort(SDelta_eq);
SFe_eqplot = sort(SFe_eq);

%plots
figure();

plot(FVeqrplot, FDelta_eqplot*180/pi, 'o-');
hold on;
plot(SVeqrplot, SDelta_eqplot*180/pi, 'x-');
set(gca, 'YDir','reverse')
title('Trim Curves')
grid on
xlabel("$$\tilde{V}_{eq}$$ [m/s]", "Interpreter", "latex", 'FontSize', 15)
ylabel("$$\delta_{eq}^*$$ [deg]", "Interpreter", "latex",  'FontSize', 15)
legend(['Flight data' newline 'x_{cg} = 0.363 [mac] '], ['Reference data' newline 'x_{cg} = 0.368 [mac]'], 'FontSize', 15)

figure();

plot(FVeqrplot, FFe_eqplot, 'o-');
hold on;
plot(SVeqrplot, SFe_eqplot, 'x-');
set(gca, 'YDir','reverse')
title('Force Curves')
grid on
xlabel("$$\tilde{V}_{eq}$$ [m/s]", "Interpreter", "latex", 'FontSize', 15)
ylabel("$$F_{e_{aer}}^{*}$$ [N]", "Interpreter", "latex", 'FontSize', 15)
legend(['Flight data' newline '\delta_{t_e} = 2.4 [deg]' newline 'x_{cg} = 0.363 [mac] '], ['Reference data' newline '\delta_{t_e} = 2.8 [deg]' newline 'x_{cg} = 0.368 [mac]'], 'FontSize', 15)


% figure();
% subplot(1,2,1);
% plot(SVeqrplot, SDelta_eqplot*180/pi, 'o-');
% set(gca, 'YDir','reverse')
% title('Trim Curve (sample)')
% grid on
% xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
% ylabel("$$\delta_{eq}^*$$ [deg]", "Interpreter", "latex")
% 
% subplot(1,2,2);
% plot(SVeqrplot, SFe_eqplot, 'x-');
% set(gca, 'YDir','reverse')
% title('Force Curve (sample)')
% grid on
% xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
% ylabel("$$F_{e_{aer}}^{*}$$ [deg]", "Interpreter", "latex")






