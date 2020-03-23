run("Parameters.m")
%%%%Sample Data
SDat=load("Sample_Elevator.mat"); %Sample Data excluding mass in kg
SMass=load("SampleMasses.mat"); 
SPMass=SMass.PeopleMasses; %Sample Mass of the people in kg

%Flight data 
FDat=load("Flight_Elevator.mat"); %Measurements from flight
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


%Sample Data
SAlt=SDat.SampleData1(:,1);    %Sample Altitude        (ft)
SVelo=SDat.SampleData1(:,2);   %Sample Velocity        (kts)
SAlpha=SDat.SampleData1(:,3);  %Sample Angle of Attack (degrees)
Sde  =SDat.SampleData1(:,4);   % Sample Deflection     (degress)

SFe = SDat.SampleData1(:,7); 
SFFl=SDat.SampleData1(:,7);    %Sample Fuel flow left  (lbs/hr)
SFFr=SDat.SampleData1(:,8);    %Sample Fuel flow right (lbs/hr)
SFUsed=SDat.SampleData1(:,9);  %Sample Fuel used       (lbs)
STemp=SDat.SampleData1(:,10);   %Measured Sample Temperature     (Celsius)

SAlt= 0.3048*SAlt ; %Sample Altitude        (m)
SVelo= 0.514444*SVelo; %Sample Velocity        (m/s)
SAlpha= pi*SAlpha/180; %Sample Angle of Attack (radians)
SFFl= 0.000125998*SFFl; %Sample Fuel flow left  (kg/s)
SFFr= 0.000125998*SFFr; %Sample Fuel flow right (kg/s)
SFUsed= 0.453592*SFUsed; %Sample Fuel used       (kg)
STemp= 273.15+STemp ;  %Measured Sample Temperature     (Kelvin)  


%Flight data 

FAlt  = FDat.Flight_Elevator(:,1);    %Sample Altitude        (ft)
FVelo =FDat.Flight_Elevator(:,2);   %Sample Velocity        (kts)
FAlpha=FDat.Flight_Elevator(:,3);  %Sample Angle of Attack (degrees)
Fde  =FDat.Flight_Elevator(:,4);   % Sample Deflection     (degress)

FFe = FDat.Flight_Elevator(:,7); 
FFFl=FDat.Flight_Elevator(:,7);    %Sample Fuel flow left  (lbs/hr)
FFFr=FDat.Flight_Elevator(:,8);    %Sample Fuel flow right (lbs/hr)
FFUsed=FDat.Flight_Elevator(:,9);  %Sample Fuel used       (lbs)
FTemp=FDat.Flight_Elevator(:,10);   %Measured Sample Temperature     (Celsius)

FAlt= 0.3048*FAlt ; %Sample Altitude        (m)
FVelo= 0.514444*FVelo; %Sample Velocity        (m/s)
FAlpha= pi*FAlpha/180; %Sample Angle of Attack (radians)
FFFl= 0.000125998*FFFl; %Sample Fuel flow left  (kg/s)
FFFr= 0.000125998*FFFr; %Sample Fuel flow right (kg/s)
FFUsed= 0.453592*FFUsed; %Sample Fuel used       (kg)
FTemp= 273.15+FTemp ;  %Measured Sample Temperature     (Kelvin



sz=size(SAlt); 

Srho=zeros(sz(1),sz(2)); %Density for Sample data points
Frho=zeros(sz(1),sz(2)); %Density for Flight Data points

SStaticP=zeros(sz(1),sz(2));  %Static pressure Sample data [Pa]
FStaticP=zeros(sz(1),sz(2));  %Static Pressure flight data [Pa]

SM=zeros(sz(1),sz(2));      %Mach number Sample data
FM=zeros(sz(1),sz(2));      %Mach number Flight data

SSound=zeros(sz(1),sz(2)); %Speed of sound Sample Data
FSound=zeros(sz(1),sz(2)); %Speed of sound flight data

STempstatic=zeros(sz(1),sz(2)); %Static Air temperature Sample Data
FTempstatic=zeros(sz(1),sz(2)); %Static Air temperature Flight Data

SVeloT=zeros(sz(1),sz(2));
FVeloT=zeros(sz(1),sz(2));

SVeloEq=zeros(sz(1),sz(2));     % Equivelent speed sample data
FVeloEq=zeros(sz(1),sz(2));     % Equivelent speed sample data



for i=1:sz(1)
    SStaticP(i,1)=p0*(1+(lambda*SAlt(i,1)/Temp0))^(-g/(lambda*R));
    FStaticP(i,1)=p0*(1+(lambda*FAlt(i,1)/Temp0))^(-g/(lambda*R));

    SM(i,1) = sqrt(2/(gamma-1)*((1+ p0/SStaticP(i,1) *((1+(gamma-1)*rho0*SVelo(i,1)^2/(2*gamma*p0))^(gamma/(gamma-1)) -1))^((gamma-1)/gamma) -1));
    FM(i,1) = sqrt(2/(gamma-1)*((1+ p0/FStaticP(i,1) *((1+(gamma-1)*rho0*FVelo(i,1)^2/(2*gamma*p0))^(gamma/(gamma-1)) -1))^((gamma-1)/gamma) -1));
    
    STempstatic(i,1)=STemp(i,1)/(1+(gamma-1)*SM(i,1)^2/2);
    FTempstatic(i,1)=FTemp(i,1)/(1+(gamma-1)*FM(i,1)^2/2);
    
    SSound(i,1)=sqrt(gamma*R*STempstatic(i,1));
    FSound(i,1)=sqrt(gamma*R*FTempstatic(i,1));
    
    SVeloT(i,1)=SM(i,1)*SSound(i,1);
    FVeloT(i,1)=FM(i,1)*FSound(i,1);
    
    Srho(i,1)=SStaticP(i,1)/(R* STempstatic(i,1));
    Frho(i,1)=FStaticP(i,1)/(R* FTempstatic(i,1));
    
    SVeloEq(i,1) = SVeloT(i,1) * sqrt(Srho(i,1)/rho0);
    FVeloEq(i,1) = FVeloT(i,1) * sqrt(Frho(i,1)/rho0);

    
end

SCL=zeros(sz(1),sz(2)); %Matrix containing CL for Sample data

for i= 1:sz(1)
   SCL(i,1)=2* (STotPreMass-SFUsed(i,1))*g/(Srho(i,1)*SVeloT(i,1)^2*S);
   FCL(i,1)=2* (FTotPreMass-FFUsed(i,1))*g/(Frho(i,1)*FVeloT(i,1)^2*S);
end




Sdelta_cg = -0.03818;        %[m] sample data  change in cg (first measumernt)
Fdelta_cg =  -0.05474  ;     %[m] flight data  change in cg (first measumernt)  

Sdelta_elev = (Sde(9,1)   -  Sde(8,1))*pi/180;             %[rad] sample change in elevator deflection 
Fdelta_elev = (Fde(9,1)   -  Fde(8,1))*pi/180;             %[rad] flight change in elevator deflection 

Sdelta_alpha= SAlpha(9,1) - SAlpha(8,1)*pi/180;            %[rad] sample change in alpha   
Fdelta_alpha= FAlpha(9,1) - FAlpha(8,1)*pi/180;            %[rad] flight change in alpha       


Scm_delta = -1/Sdelta_elev * SCL(9,1) *Sdelta_cg/c;        % sample cm_delta 
Fcm_delta = -1/Fdelta_elev * FCL(9,1) *Fdelta_cg/c;         % Flight  cm_delta 


Scm_alpha = -Sdelta_elev/Sdelta_alpha *Scm_delta;          % Sample cm_alpha 
Fcm_alpha = -Fdelta_elev/Fdelta_alpha *Fcm_delta;           % Flight cm_alpha 




