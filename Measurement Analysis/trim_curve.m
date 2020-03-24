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
Ws = 60500/9.80665;      % standard aircraft weight [N]

Scm_delta = -1.046329977000884;
Fcm_delta = -1.366101036028307;

Scm_alpha = -0.100463799265000;
Fcm_alpha =   -0.138995602327464;


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
SVelo=SDat.SampleData1(1:7,2);   %Sample Velocity               (kts)
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
FVelo=FDat.Data(1:7,2);  %Flight Velocity               (m/s)
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


STempISA=[278.23; 278.21; 278.21; 278.21; 278.21; 278.04; 278.04];
FTempISA=[264.4; 264.42; 264.4; 264.4; 264.4; 264.4; 264.4]; %CHANGE 

sz=size(SAlt); 

SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data
FCL=zeros(sz(1),sz(2)); %Matrix containing CL for Flight data
for i= 1:sz(1)
   SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVelo(i,1)^2*S);
   FCL(i,1)=2* (FTotPreMass-FFUsed(i,1))*g/(Frho(i,1)*FVelo(i,1)^2*S);
end


SDeltaTemp= STemp-STempISA;
FDeltaTemp= FTemp-FTempISA;


% equivelent airspeed calculations 

Fp = ones(sz(1),sz(2));
Sp = zeros(sz(1),sz(2));

% pressure at different altitudes 
for i= 1:sz(1)
    Sp(i,1) = p0*(1+lambda*SAlt(i,1)/Temp0)^(-g/(lambda*R));
    Fp(i,1) = p0*(1+lambda*FAlt(i,1)/Temp0)^(-g/(lambda*R));
end

%Mach number at different altitudes
FM = zeros(sz(1), sz(2)); 
SM = zeros(sz(1), sz(2));
    
for i= 1:sz(1)
    SM(i,1) = sqrt(2/(gamma-1) *((1+ p0/Sp(i,1) *((1+ (gamma-1)*rho0*SVelo(i,1)^2/(2*gamma *p0))^(gamma/(gamma-1))-1))^((gamma -1)/gamma) -1));
    FM(i,1) = sqrt(2/(gamma-1) *((1+ p0/Fp(i,1) *((1+ (gamma-1)*rho0*FVelo(i,1)^2/(2*gamma *p0))^(gamma/(gamma-1))-1))^((gamma -1)/gamma) -1));
end 


% correcting the measured temperature. 
FT = zeros(sz(1), sz(2));
ST = zeros(sz(1),sz(2));

for i= 1:sz(1)
    FT(i,1) = FTemp(i,1)/(1+(gamma -1)/2 *FM(i,1)^2);
    ST(i,1) = STemp(i,1)/(1+(gamma -1)/2 *SM(i,1)^2);
end

% calculating the density using ideal gas low
Srho= zeros(sz(1),sz(2));
Frho= zeros(sz(1),sz(2));   

for i= 1:sz(1)%used to be Srho1 and Frho2
    Srho(i,1) = Sp(i,1)/(R*ST(i,1));
    Frho(i,1) = Fp(i,1)/(R*FT(i,1)); 
end

% calculating the equivelant airspeed 
FVeq = zeros(sz(1),sz(2));
SVeq = zeros(sz(1),sz(2));

for i = 1:sz(1)
    FVeq(i,1) = sqrt(gamma*R*FT(i,1))*FM(i,1)*sqrt(Frho(i,1)/rho0); % true airspeed = M*a , Veq = Vt*sqrt(rho/rho0)
    SVeq(i,1) = sqrt(gamma*R*ST(i,1))*SM(i,1)*sqrt(Srho(i,1)/rho0);
end


%Actual Thrust
SampleThrustParams = [SAlt SM SDeltaTemp SFFl SFFr]; 
FlightThrustParams = [FAlt FM FDeltaTemp FFFl FFFr];

save TrimCurve_Thrust\STP SampleThrustParams;
save TrimCurve_Thrust\FTP FlightThrustParams;

run('TrimCurve_Thrust\Thruster');

S_ActualThrustCoef = (S_Thrust(:, 1) + S_Thrust(:, 2))./(0.5*(D^2)*Srho.*SVelo.^2);
F_ActualThrustCoef = (F_Thrust(:, 1) + F_Thrust(:, 2))./(0.5*(D^2)*Frho.*FVelo.^2);


%Reduced thrust
%Standard Fuel flow
FFs = 0.048*ones(7,1);

SampleThrustParams = [SAlt SM SDeltaTemp FFs FFs]; 
FlightThrustParams = [FAlt FM FDeltaTemp FFs FFs];

save TrimCurve_Thrust\STP SampleThrustParams;
save TrimCurve_Thrust\FTP FlightThrustParams;

run('TrimCurve_Thrust\Thruster.m');

SVeqt = SVeq.*sqrt(Ws./(STotPreMass - SFUsed));
FVeqt = FVeq.*sqrt(Ws./(FTotPreMass - FFUsed));


S_ReducedThrustCoef = (S_Thrust(:, 1) + S_Thrust(:, 2))./(0.5*(D^2)*Srho.*SVelo.^2);
F_ReducedThrustCoef = (F_Thrust(:, 1) + F_Thrust(:, 2))./(0.5*(D^2)*Frho.*FVelo.^2);


%trim curve
SDelta_eq = SDelta - (CmTc/Scm_delta)*(S_ReducedThrustCoef - S_ActualThrustCoef);
FDelta_eq = FDelta - (CmTc/Fcm_delta)*(F_ReducedThrustCoef - F_ActualThrustCoef);

SVeqt = SVeq.*sqrt(Ws./(STotPreMass - SFUsed));
FVeqt = FVeq.*sqrt(Ws./(FTotPreMass - FFUsed));


%force curve
SFe_eq = SFe.*(Ws./(STotPreMass - SFUsed));
FFe_eq = FFe.*(Ws./(FTotPreMass - FFUsed));


%plots
figure();
subplot(1,2,1);
plot(FVeqt, FDelta_eq*180/pi, 'o');
set(gca, 'YDir','reverse')
title('Trim Curve (flight)')
grid on
xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
ylabel("$$\delta_{eq}^*$$ [deg]", "Interpreter", "latex")

subplot(1,2,2);
plot(FVeqt, FFe_eq, 'x');
set(gca, 'YDir','reverse')
title('Force Curve (flight)')
grid on
xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
ylabel("$$F_{e_{aer}}^{*}$$ [deg]", "Interpreter", "latex")


% figure();
% subplot(1,2,1);
% plot(SVeqt, SDelta_eq*180/pi, 'o');
% set(gca, 'YDir','reverse')
% title('Trim Curve (sample)')
% grid on
% xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
% ylabel("$$\delta_{eq}^*$$ [deg]", "Interpreter", "latex")
% 
% subplot(1,2,2);
% plot(SVeqt, SFe_eq, 'x');
% set(gca, 'YDir','reverse')
% title('Force Curve (sample)')
% grid on
% xlabel("$$\tilde{V_{eq}}$$ [m/s]", "Interpreter", "latex")
% ylabel("$$F_{e_{aer}}^{*}$$ [deg]", "Interpreter", "latex")






