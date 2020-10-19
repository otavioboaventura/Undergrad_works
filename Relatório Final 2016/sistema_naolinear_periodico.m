function yprime = sistema_naolinear_periodico(t,z)
% sistema linear com oscilador não ideal
yprime=zeros(4,1);
% =================== Parâmetros  ====================
c1= 0.002;
m1= 1;
knl1= 0.5;

c2= 0.002;
m2= 0.05;
knl2= 0.0025;

u=0.5;
E=0.5;
alpha=2.9; %2.5 - 2.9
B=1.5;

% Igualar de  acordo com os deslocamentos
F = 4.5;
w = 1; 

% ============================= Sistema em forma de estado ========================
yprime(1)=z(2);
yprime(2)= (-c1*z(2) - c2*(z(2)-z(4)) - knl1*z(1)^3 - knl2*(z(1)-z(3))^3)+ F*sin(w*t);
yprime(3)=z(4);
yprime(4)= - c2/m2*(z(4)-z(2)) - knl2/m2*(z(3)-z(1))^3;
% ================================================================