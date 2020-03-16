function dx=EoM(t,x)
%% Documentation
%Classification: Unclassified
%Program writen by: Brian Wade
%Date 9 Sept 2015

%% Description and setup
%This program runs as a subprogram of the TBM_Flight.m program.  This
%function solves the 12 simultatious ordinary differential equations that
%describe the ballistic motion of a projectile. These equations are based 
%on the equations found in McCoy, 1998 (see citation below).  The input
%aerodynamic coefficeints are from McCoy's book as well.

%Input data from: 
%1. McCoy, RL, Modern Exerioor Balistics: The Launch and Flight Dynamics 
%of Symmetric Projectiles, Schiffer Military History, Atglen, PA, 1998.

%% Program
%Get globals from the main program.
global cdo_wpn
global clo_wpn
global cmo_wpn
global cmqao_wpn
global cd2_wpn
global cl2_wpn
global cm2_wpn
global cmqa2_wpn
global std_atm
global d
global It
global m
global R
global gravity

%Find atmospheric variables
rho=interp1(std_atm(:,1),std_atm(:,7),x(11)/1000); %atmpospheric density in
%kg/m^3

a=interp1(std_atm(:,1),std_atm(:,8),x(11)/1000); %speed of sound in m/s

%Find refrence area of projectile.
S=(pi/4)*d^2; %refrence area in m^2

%% Body Forces
%find total velocity (does not account for winds)
V=sqrt(x(1)^2+x(2)^2+x(3)^2);

%Find total angle of attack
alpha_t=acos((x(1)*x(7)+x(2)*x(8)+x(3)*x(9))/V);

%find mach number
mach=V/a;

%find dynamic coefficient
q=rho*V*S;

%Aeroforces - all in nonrotating frame
%lookup cd,cl,cm,cmqa
if mach>cdo_wpn(size(cdo_wpn,1),1)
    M=cdo_wpn(size(cdo_wpn,1),1);
else
    M=mach;
end

%Look up the aerodynamic coefficients
cdo=interp1(cdo_wpn(:,1),cdo_wpn(:,2),M);
cd2=interp1(cd2_wpn(:,1),cd2_wpn(:,2),M);
clo=interp1(clo_wpn(:,1),clo_wpn(:,2),M);
cl2=interp1(cl2_wpn(:,1),cl2_wpn(:,2),M);
cmo=interp1(cmo_wpn(:,1),cmo_wpn(:,2),M);
cm2=interp1(cm2_wpn(:,1),cm2_wpn(:,2),M);
cmqao=interp1(cmqao_wpn(:,1),cmqao_wpn(:,2),M);
cmqa2=interp1(cmqa2_wpn(:,1),cmqa2_wpn(:,2),M);

cd=cdo+cd2*(sin(alpha_t))^2; %Drag
cl=clo+cl2*(sin(alpha_t))^2; %lift
cm=cmo+cm2*(sin(alpha_t))^2; %moment
cmqa=cmqao+cmqa2*(sin(alpha_t))^2; %overturning moment

Cd=(q*cd)/(2*m); %Body drag
Cl=(q*cl)/(2*m); %Body lift
Cm=(q*d*cm)/(2*It); %body moment
Cmq=(q*d^2*cmqa)/(2*It); %Body overturning moment

%gravity compenents in each axis

g1=-gravity*x(10)/R;
g2=-gravity*(1-(2*x(11)/R));
g3=0;
% g1=0;
% g2=gravity;
% g3=0;

%Intermediate term for below equations
IpP_It=x(4)*x(7)+x(5)*x(8)+x(6)*x(9);

%Equations of motion:
    %x(1) = x-velocity with respect (wrt) intertial frame (m/s).
    %x(2) = y-velocity wrt intertial frame (m/s).
    %x(3) = z-velocity wrt intertial frame (m/s).
    %x(4) = roll rate wrt intertial frame (rad/s).
    %x(5) = pitch rate wrt intertial frame (m/s).
    %x(6) = yaw rate wrt intertial frame (m/s).
    %x(7) = x component of projectile unit vector wrt intertial frame
    %x(8) = y component of projectile unit vector wrt intertial frame
    %x(9) = z component of projectile unit vector wrt intertial frame
    %x(10) = position of munition center of gravity (CG) wrt intertial
    %frame x-axis (m).  This is the range.
    %x(11) = position of munition CG wrt intertial frame y-axis (m). 
    %This is the altitude.
    %x(12) = position of munition CG wrt intertial frame z-axis (m). 
    %This is the cross-range.
    

dx(1,1)=-Cd*x(1)+Cl*(V^2*x(7)-V*x(1)*cos(alpha_t))+g1;
dx(2,1)=-Cd*x(2)+Cl*(V^2*x(8)-V*x(2)*cos(alpha_t))+g2;
dx(3,1)=-Cd*x(3)+Cl*(V^2*x(9)-V*x(3)*cos(alpha_t))+g3;

dx(4,1)=Cm*(x(2)*x(9)-x(3)*x(8))+Cmq*(x(4)-IpP_It*x(7));
dx(5,1)=Cm*(x(3)*x(7)-x(1)*x(9))+Cmq*(x(5)-IpP_It*x(8));
dx(6,1)=Cm*(x(1)*x(8)-x(2)*x(7))+Cmq*(x(6)-IpP_It*x(9));

dx(7,1)=x(5)*x(9)-x(6)*x(8);
dx(8,1)=x(6)*x(7)-x(4)*x(9);
dx(9,1)=x(4)*x(8)-x(5)*x(7);

dx(10,1)=x(1); 
dx(11,1)=x(2)+(x(1)^2)/(2*R); 
dx(12,1)=x(3); 

end
