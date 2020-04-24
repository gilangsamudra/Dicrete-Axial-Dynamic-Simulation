clear ;
% close all;
clc;

load('Input Accel.mat');
load('Viscous Force 1KM 4Segments.mat');
Hookawal = ones(1,31)*1;
hookmin = ones(1,21)*-0;
nol = zeros(1,10);
pengurang = [nol hookmin];
Hook = Trajectory(2:end);

%% DRAG FORCES
Loaded_Data = readtable('Drag 1KM 4 Segments.xlsx'); % in table format
drag        = Loaded_Data{:,5};                      % drag forces (lbf)
drag        = drag(~isnan(drag));

% Panjang in meter
Total = 1000;
ns    = 4;
LengthS = Total/ns;
mod_el = 206842718795.3;          % Elastic modulus, N/m^2

% String Table, inches
ODm = 5 * 0.0254; %in meter now
ODj = 6.5 * 0.0254;
ODc = 6.5 * 0.0254;
ODb = 6.5 * 0.0254;
IDm = 4.276 * 0.0254;
IDj = 3.75 * 0.0254;
IDc = 2.813 * 0.0254;
IDb = 3.89 * 0.0254;
holeD = (7+(0/8)) * 0.0254;

% in SI unit
pipe_den = 7800.001722079;        %Pipe density, in kg/m^3
rho  = 1840;                      %Mud density, in kg/m^3

alp = ODm/holeD;
Kc = (alp^2 - sqrt((alp^4 +alp)/(1+alp)))/(1-alp^2);
n = 0.55; 
k = 0.13;
h1 = (holeD-ODm)/2; %The width of the annular space
h2 = (holeD-ODj)/2;
hC = (holeD-ODc)/2;
hB = (holeD-ODb)/2;
aa = 0.138; %Rheology parameter for annulus
bb = 0.312;
Au1 = 2*pi*0.95*LengthS*(ODm/2 + holeD/2);
Au2 = 2*pi*0.05*LengthS*(ODj/2 + holeD/2);
AuC = 2*pi*0.911*LengthS*(ODc/2 + holeD/2);
AuB = 2*pi*0.088*LengthS*(ODb/2 + holeD/2);

% Drillstring Parameters
Stiff_bha = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthS);
Stiff_coll= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthS);
Stiff_main= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthS);
Stiff_sec = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthS);

stiffness1 = Stiff_main*Stiff_sec/(Stiff_main+Stiff_sec);     %stiffness of drill pipe
stiffness2 = stiffness1;
stiffness3 = stiffness1;
stiffness4 = Stiff_bha*Stiff_coll/(Stiff_bha+Stiff_coll);
stiffness = [stiffness1 stiffness2 stiffness3 stiffness4];

mass_top  = 20000;
mass_bha  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthS)*pipe_den;        
mass_coll = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthS)*pipe_den;
mass_main = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthS)*pipe_den;
mass_sec  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthS)*pipe_den;

mass1 = mass_main+mass_sec; %mass of drill pipe
mass2 = mass1;
mass3 = mass1;
mass4 = mass_bha+mass_coll;
mass  = [mass_top mass1 mass2 mass3 mass4];


r1 = stiffness1/mass_top;
r12= stiffness1/mass1;
r2 = stiffness2/mass1;
r23= stiffness2/mass2;
r3 = stiffness3/mass2;
r34= stiffness3/mass3;
r4 = stiffness4/mass3;
r45= stiffness4/mass4;
r = [r1 r12 r2 r23 r3 r34 r4 r45];

A = [0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0
    r12 0 (-r12-r2) -pp1(1)/mass1 r2 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0
    0 0 r23 0 (-r23-r3) -pp2(1)/mass2 r3 0 0 0
    0 0 0 0 0 0 0 1 0 0
    0 0 0 0 r34 0 (-r34-r4) -pp3(1)/mass3 r4 0
    0 0 0 0 0 0 0 0 0 1
    0 0 0 0 0 0 r45 0 -r45 -pp4(1)/mass4];

% Matrix B
B = [0 0 0 0 0 0 0 0 0 0
    1 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 -1/mass1 0 0 0 -1/mass1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 -1/mass2 0 0 0 -1/mass2 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 -1/mass3 0 0 0 -1/mass3 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 0 -1/mass4 0 0 0 -1/mass4 0];
% 
C = eye(10);
D = zeros(10,10);
fre = 50;

t = 0:(1/fre):29;
drag1 = drag(1)*ones(numel(t),1);
drag2 = drag(2)*ones(numel(t),1);
drag3 = drag(3)*ones(numel(t),1);
drag4 = drag(4)*ones(numel(t),1);
vis1  = pp1(2)*ones(numel(t),1);
vis2  = pp2(2)*ones(numel(t),1);
vis3  = pp3(2)*ones(numel(t),1);
vis4  = pp4(2)*ones(numel(t),1);

hook2 = [Hook(1)];
for i=1:numel(Hook)-1
    interpol = linspace(Hook(i),Hook(i+1),fre);
    hook2 = [hook2 interpol];
end
input = [hook2' drag1 drag2 drag3 drag4 vis1 vis2 vis3 vis4 zeros(numel(t),1)];

syscon = ss(A,B,C,D);
xo = zeros(10,1);
modeld2 = c2d(syscon,(1/fre));
z=lsim(modeld2,input,t,xo);
