rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
S = 30.00; %(m^2)

%%%%Data
SDat=load("SampleData.mat"); %Sample Data excluding mass in kg
SMass=load("SampleMasses.mat"); 
SPMass=SMass.PeopleMasses; %Sample Mass of the people in kg


FDat=load("LiftFlightData.mat"); %Measurements from flight
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
SAlt=SDat.SampleData(:,1);    %Sample Altitude        (ft)
SVelo=SDat.SampleData(:,2);   %Sample Velocity        (kts)
SAlpha=SDat.SampleData(:,3);  %Sample Angle of Attack (degrees)
SFFl=SDat.SampleData(:,4);    %Sample Fuel flow left  (lbs/hr)
SFFr=SDat.SampleData(:,5);    %Sample Fuel flow right (lbs/hr)
SFUsed=SDat.SampleData(:,6);  %Sample Fuel used       (lbs)
STemp=SDat.SampleData(:,7);   %Sample Temperature     (Celsius)

SAlt= 0.3048*SAlt ; %Sample Altitude        (m)
SVelo= 0.514444*SVelo; %Sample Velocity        (m/s)
SAlpha= pi*SAlpha/180; %Sample Angle of Attack (radians)
SFFl= 0.000125998*SFFl; %Sample Fuel flow left  (kg/s)
SFFr= 0.000125998*SFFr; %Sample Fuel flow right (kg/s)
SFUsed=0.453592*SFUsed; %Sample Fuel used       (kg)
STemp=273.15+STemp ;  %Sample Temperature     (Kelvin)  

%FLIGHT DATA
FAlt=FDat.LiftFlightData(:,1);  %Flight Altitude (m)
FVelo=FDat.LiftFlightData(:,2); %Flight Velocity (m/s)
FAlpha=FDat.LiftFlightData(:,3); %Flight Angle of Attack (degrees)
FFFl=FDat.LiftFlightData(:,4);   %Flight Fuel Flow left (lbs/hr)
FFFr=FDat.LiftFlightData(:,5);   %Flight Fuel flow right (lbs/hr)
FFUsed=FDat.LiftFlightData(:,6);  %Flight Fuel Used (lbs)
FTemp=FDat.LiftFlightData(:,7);   %Flight Temperature (Celsius)

FAlt=0.3048*FAlt ; %Flight Altitude (m)
FVelo=0.514444*FVelo ; %Flight Velocity (m/s)
FAlpha=pi*FAlpha/180 ; %Flight AoA (radians)
FFFl=0.000125998*FFFl ; %Flight Fuel Flow Left (kg/s)
FFFr=0.000125998*FFFr ; %Flight Fuel Flow right (kg/s)
FFUsed=0.453592*FFUsed ; %Flight Fuel used (kg)
FTemp=FTemp+273.15 ; %(Kelvin)


STempISA=[278.23; 278.21; 278.21; 278.21; 278.21; 278.04];
FTempISA=[264.4; 264.42; 264.4; 264.4; 264.4; 264.4;];

sz=size(SAlt); 

Srho=zeros(sz(1),sz(2)); %Density for Sample data points
Frho=zeros(sz(1),sz(2)); %Density for Flight Data points
for i= 1:sz(1)
    Srho(i,1)=rho0*((1+(lambda*SAlt(i)/Temp0)))^(-((g/(lambda*R))+1));
    Frho(i,1)=rho0*((1+(lambda*FAlt(i)/Temp0)))^(-((g/(lambda*R))+1));
end

SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data
FCL=zeros(sz(1),sz(2)); %Matrix containing CL for Flight data
for i= 1:sz(1)
   SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVelo(i,1)^2*S);
   FCL(i,1)=2* (FTotPreMass-FFUsed(i,1))*g/(Frho(i,1)*FVelo(i,1)^2*S);
end


SDeltaTemp= STemp-STempISA;
FDeltaTemp= FTemp-FTempISA;



