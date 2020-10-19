function yprime = sistema_naolinear_naoideal(t,z)

yprime=zeros(6,1);
% =================== Parâmetros  ====================
c1= 0.002;
m1= 1;
knl1= 0.45;

c2= 0.002;
m2= 0.05;
knl2= 0.003;

u=0.5;
E=0.5;
alpha=2.9; %2.5 - 2.9
B=1.5;

% ============================= Sistema ========================
yprime(1)=z(2);
yprime(2)= (-c1*z(2) - c2*(z(2)-z(4)) - knl1*z(1)^3 - knl2*(z(1)-z(3))^3)+u*z(6)^2*cos(z(5))+u*sin(z(5))*(alpha-B*z(6))/(m1-u*E*(sin(z(5))^2));
yprime(3)=z(4);
yprime(4)= - c2/m2*(z(4)-z(2)) - knl2/m2*(z(3)-z(1))^3;
yprime(5)=z(6);
yprime(6)=(E*sin(z(5))*(-knl1*z(1)-knl2*(z(1)-z(3)))+u*(z(6)^2)*E*cos(z(5))*sin(z(5))+u*E*sin(z(5)^2)+(m1-u*E*sin(z(5))^2)*(alpha-B*z(6)))/(m1-u*E*(sin(z(5))^2));

% ================================================================

% yprime(2)= (-c1*z(2) - c2*(z(2)-z(4)) - knl1*z(1)^3 - knl2*(z(1)-z(3))^3)+u*z(6)^2*cos((pi/180)*z(5))+u*sin((pi/180)*z(5))*(alpha-B*z(6))/(m1-u*E*(sin((pi/180)*z(5))^2)); % graus
% yprime(6)=(E*sin((pi/180)*z(5))*(-knl1*z(1)-knl2*(z(1)-z(3)))+u*z(6)^2*E*cos((pi/180)*z(5))*sin((pi/180)*z(5))+u*E*sin((pi/180)*z(5)^2)+(m1-u*E*sin((pi/180)*z(5))^2)*(alpha-B*z(6)))/(m1-u*E*(sin((pi/180)*z(5))^2)) % graus

         