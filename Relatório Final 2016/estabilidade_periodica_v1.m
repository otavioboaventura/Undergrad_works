clc;
clear all;
close all;
format long;

% Definição do sistema simbolicamente
syms x1 x2 x3 x4 knl1 knl2 c1 m1 c2 m2 F w t;
 
x11=x2;
x22= (-c1*x2 - c2*(x2-x4) - knl1*x1^3 - knl2*(x1-x3)^3)+ F*cos(w*t);
x33=x4;
x44= - c2/m2*(x4-x2) - knl2/m2*(x3-x1)^3;
 
%Cálculo da matriz Jacobiana
J=jacobian([x11; x22; x33 ; x44],[x1 x2 x3 x4]);
 
%========== PARAMETROS ===============
c1= 0.002;
m1= 1;

c2= 0.002;
m2= 0.05;

u=0.5;
E=0.5;
alpha=2.9; %2.5 - 2.9
B=1.5;

F = 1; 
w = 4;  
t = 1;

%Condições iniciais
knl1= 0.4;   % k11=0.5;
knl2= 0.001; % knl2= 0.0025;
knl2in = 0.001;
x1=0.1;
x2=0;
x3=0;
x4=0;

% Definição do passo de knl1 e knl2
passo1 = 0.01;
passo2 = 0.00015;
 
%Análise da Estabilidade da Estrutura
while knl1<=(0.61)
 while knl2<=(0.004)
j11= 0; 
j12= 1;
j13= 0;
j14= 0;
j21= -3*knl1*x1^2-3*knl2*(x1-x3)^2; 
j22= -c1-c2;
j23= 3*knl2*(x1-x3)^2;
j24= c2;
j31= 0; 
j32= 0;
j33= 0;
j34= 1;
j41= 3*knl2/m2*(x3-x1)^2; 
j42= c2/m2;
j43= -3*knl2/m2*(x3-x1)^2;
j44= -c2/m2;
 
  A=[j11 j12 j13 j14;j21 j22 j23 j24;j31 j32 j33 j34;j41 j42 j43 j44];
     fprintf('%f,%f',knl1,knl2)
     autovalor=eig(A)

%========== PARA A ESTRUTURA ===============   
 if (autovalor(1) ~= autovalor(2)) & abs(autovalor(1)*autovalor(2))>1e-7 & (abs(imag(autovalor(1)))<1e-7 & abs(imag(autovalor(2)))<1e-7) % diferentes, nenhum nulo, reais
    if (real(autovalor(1))*real(autovalor(2)))>1e-7 & real(autovalor(1))<-1e-7 & real(autovalor(2))<-1e-7 % parte real negativa, no hiperbolico assintoticamente estavel
      figure(1)
      hold on
      plot(knl1,knl2,'go') %estavel
    end
    if (real(autovalor(1))*real(autovalor(2)))>1e-7 & real(autovalor(1))>1e-7 & real(autovalor(2))>1e-7 % no hiperbolico instavel
      figure(1)
      hold on
      plot(knl1,knl2,'g*') %instavel
    end
    if (real(autovalor(1))*real(autovalor(2)))<-1e-7 % sela hiperbolica instavel
      figure(1)
      hold on
      plot(knl1,knl2,'m*') %instavel
    end
 end

 if  imag(autovalor(1)) == -imag(autovalor(2)) & abs(imag(autovalor(1)))>1e-7 & abs(imag(autovalor(2)))>1e-7 % complexo conjugado, complexos
  
  if abs(real(autovalor(1)))>1e-7 & abs(real(autovalor(2)))>1e-7 % re(v1,v2!=0)
      
      if real(autovalor(1))>1e-7 & real(autovalor(2))>1e-7
          figure(1)
          hold on
          plot(knl1,knl2,'b*') %instavel - foco hiperbolico
      end
      if real(autovalor(1))<-1e-7 & real(autovalor(2))<-1e-7
          figure(1)
          hold on
          plot(knl1,knl2,'bo') %estavel - foco hiperbolico
     end
  end
  
  if abs(real(autovalor(1)))<1e-7 & abs(real(autovalor(2)))<1e-7
      figure(1)  
      hold on
      plot(knl1,knl2,'ko') %estavel - centro eliptico 
  end
 end

  if abs(autovalor(1)*autovalor(2))<1e-7 %v1*v2=0 Nulo - caso degenerados
      figure(1)
      hold on
      plot(knl1,knl2,'rs') % Degenerado (nem estável ou instável)
  end

  if (autovalor(1) == autovalor(2) & abs(real(autovalor(1)))>1e-7 & abs(real(autovalor(2)))>1e-7 & abs(imag(autovalor(1)))<1e-7 & abs(imag(autovalor(2)))<1e-7) % v1=v2, re(v1,v2)>0, eR
      if real(autovalor(1))>1e-7 % Nó impróprio - instavel
        figure(1)
        hold on
        plot(knl1,k21,'y*') %instavel 
     end 
     if real(autovalor(1))<-1e-7  % Nó impróprio - assintoticamente estavel
        figure(1)
        hold on
        plot(knl1,knl2,'yo') %estavel
     end
  end
  
%========== PARA O TMD ===============
  if (autovalor(3) ~= autovalor(4)) & abs(autovalor(3)*autovalor(4))>1e-7 & (abs(imag(autovalor(3)))<1e-7 & abs(imag(autovalor(4)))<1e-7) % diferentes, nenhum nulo, reais
    if (real(autovalor(3))*real(autovalor(4)))>1e-7 & real(autovalor(3))<-1e-7 & real(autovalor(4))<-1e-7 % parte real negativa, no hiperbolico assintoticamente estavel
      figure(2)
      hold on
      plot(knl1,knl2,'go') %estavel
    end
    if (real(autovalor(3))*real(autovalor(4)))>1e-7 & real(autovalor(3))>1e-7 & real(autovalor(4))>1e-7 % no hiperbolico instavel
      figure(2)
      hold on
      plot(knl1,knl2,'g*') %instavel
    end
    if (real(autovalor(3))*real(autovalor(4)))<-1e-7 % sela hiperbolica instavel
      figure(2)
      hold on
      plot(knl1,knl2,'m*') %instavel
    end
 end

 if  imag(autovalor(3)) == -imag(autovalor(4)) & abs(imag(autovalor(3)))>1e-7 & abs(imag(autovalor(4)))>1e-7 % complexo conjugado, complexos
  
  if abs(real(autovalor(3)))>1e-7 & abs(real(autovalor(4)))>1e-7 % re(v1,v2!=0)
      
      if real(autovalor(3))>1e-7 & real(autovalor(4))>1e-7
          figure(2)
          hold on
          plot(knl1,knl2,'b*') %instavel - foco hiperbolico
      end
      if real(autovalor(3))<-1e-7 & real(autovalor(4))<-1e-7
          figure(2)
          hold on
          plot(knl1,knl2,'bo') %estavel - foco hiperbolico
     end
  end
  
  if abs(real(autovalor(3)))<1e-7 & abs(real(autovalor(4)))<1e-7
      figure(2)  
      hold on
      plot(knl1,knl2,'ko') %estavel - centro eliptico 
  end
 end

  if abs(autovalor(3)*autovalor(4))<1e-7 %v1*v2=0 Nulo - caso degenerados
      figure(2)
      hold on
      plot(knl1,knl2,'rs') % Degenerado (nem estável ou instável)
  end

  if (autovalor(3) == autovalor(4) & abs(real(autovalor(3)))>1e-7 & abs(real(autovalor(4)))>1e-7 & abs(imag(autovalor(3)))<1e-7 & abs(imag(autovalor(4)))<1e-7) % v1=v2, re(v1,v2)>0, eR
      if real(autovalor(3))>1e-7 % Nó impróprio - instavel
        figure(2)
        hold on
        plot(knl1,k21,'y*') %instavel 
     end 
     if real(autovalor(3))<-1e-7  % Nó impróprio - assintoticamente estavel
        figure(2)
        hold on
        plot(knl1,knl2,'yo') %estavel
     end
  end
        knl2=knl2+passo2;
  end
    knl2=knl2in;
    knl1=knl1+passo1;
end

figure(1)
title('Diagrama da Estabilidade do Oscilador Não Ideal (Massa 1)');
xlabel('Rigidez k_{nl1}')
ylabel('Rigidez k_{nl2}')

figure(2)
title('Diagrama da Estabilidade do TMD (Massa 2)');
xlabel('Rigidez k_{nl1}')
ylabel('Rigidez k_{nl2}')

  
