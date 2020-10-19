clc;
clear all;
close all;
format long;

% Parâmetros do sistema
c1= 0.002;
m1= 1;

c2= 0.002;
m2= 0.05;

u=0.5;
E=0.5;
alpha=2.9; %2.5 - 2.9
B=1.5;


% % Do artigo
% c1= 0.002;
% m1= 1;
% 
% c2= 0.002;
% m2= 0.5;
% 
% u=0.5;  % u = raio
% E=0.06; % E= u*d (raio x inercia)
% alpha=5; %2.5 - 2.9 
% B=1.5;


% Definição do sistema simbolicamente
syms x1 x2 x3 x4 x5 x6 knl1 knl2; % variáveis simbólicas
x11=x2;
x22= (-c1*x2 - c2*(x2-x4) - knl1*x1^3 - knl2*(x1-x3)^3)+u*x6^2*cos(x5)+u*sin(x5)*(alpha-B*x6)/(m1-u*E*(sin(x5)^2));
x33=x4;
x44= - c2/m2*(x4-x2) - knl2/m2*(x3-x1)^3;
x55=x6;
x66=(E*sin(x5)*(-knl1*x1-knl2*(x1-x3))+u*x6^2*E*cos(x5)*sin(x5)+u*E*sin(x5^2)+(m1-u*E*sin(x5)^2)*(alpha-B*x6))/(m1-u*E*(sin(x5)^2));

% Calculo da matriz Jacobina
J=jacobian([x11; x22; x33 ; x44; x55; x66],[x1 x2 x3 x4 x5 x6]);

% Calculo dos autovalores analíticos
autovalor1=eig(J)

% Definição do passo de knl1 e knl2
passo1 = 0.01;
passo2 = 0.00015;

% Substituição das variáveis simbólicas por variáveis float
knl1= 0.4;   % k11=0.5;
knl2= 0.001; % knl2= 0.0025;
knl2in = 0.001;
x1=0.1;
x2=0;
x3=0;
x4=0;
x5=0.7;
x6=0;

  disp('Plotando Grafico, aguarde...');
while knl1<=(0.61)
 while knl2<=(0.004)
    % Substituição das variáveis simbólicas por variáveis float
j11= 0; 
j12= 1;
j13= 0;
j14= 0;
j15= 0;
j16= 0;
j21=  - 3*knl1*x1^2 - 3*knl2*(x1 - x3)^2; 
j22= -1/250;
j23=  3*knl2*(x1 - x3)^2;
j24= 1/500;
j25= -1/2*x6^2*sin(x5)+1/2*cos(x5)*(29/10-3/2*x6)/(1-1/4*sin(x5)^2)+1/4*sin(x5)^2*(29/10-3/2*x6)/(1-1/4*sin(x5)^2)^2*cos(x5);
j26= x6*cos(x5)-3/4*sin(x5)/(1-1/4*sin(x5)^2);
j31= 0; 
j32= 0;
j33= 0;
j34= 1;
j35= 0; 
j36= 0;
j41=  60*knl2*(x1 - x3)^2; 
j42= 1/25;
j43= -60*knl2*(x1 - x3)^2;
j44= -1/25;
j45= 0;
j46= 0;
j51= 0;
j52= 0;
j53= 0;
j54= 0;
j55= 0;
j56= 1;
j61= 1/2*sin(x5)*(-knl1-knl2)/(1-1/4*sin(x5)^2);
j62= 0;
j63= 1/2*sin(x5)*knl2/(1-1/4*sin(x5)^2);
j64= 0;
j65= (1/2*cos(x5)*(-knl1*x1-knl2*(x1-x3))-1/4*x6^2*sin(x5)^2+1/4*x6^2*cos(x5)^2+1/2*cos(x5^2)*x5-1/2*sin(x5)*cos(x5)*(29/10-3/2*x6))/(1-1/4*sin(x5)^2)+1/2*(1/2*sin(x5)*(-knl1*x1-knl2*(x1-x3))+1/4*x6^2*cos(x5)*sin(x5)+1/4*sin(x5^2)+(1-1/4*sin(x5)^2)*(29/10-3/2*x6))/(1-1/4*sin(x5)^2)^2*sin(x5)*cos(x5);
j66= (1/2*x6*cos(x5)*sin(x5)-3/2+3/8*sin(x5)^2)/(1-1/4*sin(x5)^2);
    A = [j11 j12 j13 j14 j15 j16;
         j21 j22 j23 j24 j25 j26;
         j31 j32 j33 j34 j35 j36;
         j41 j42 j43 j44 j45 j46;
         j51 j52 j53 j54 j55 j56;
         j61 j62 j63 j64 j65 j66];
    
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
      plot(knl1,knl2,'c*') %instavel
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
      plot(knl1,knl2,'c*') %instavel
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
      figure(1)
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
  
%==========  PARA A FONTE    =================
 if (autovalor(5) ~= autovalor(6)) & abs(autovalor(5)*autovalor(6))>1e-7 & (abs(imag(autovalor(5)))<1e-7 & abs(imag(autovalor(6)))<1e-7) % diferentes, nenhum nulo, reais
    if (real(autovalor(5))*real(autovalor(6)))>1e-7 & real(autovalor(5))<-1e-7 & real(autovalor(6))<-1e-7 % parte real negativa, no hiperbolico assintoticamente estavel
      figure(3)
      hold on
      plot(knl1,knl2,'go') %estavel
    end
    if (real(autovalor(5))*real(autovalor(6)))>1e-7 & real(autovalor(5))>1e-7 & real(autovalor(6))>1e-7 % no hiperbolico instavel
      figure(3)
      hold on
      plot(knl1,knl2,'g*') %instavel
    end
    if (real(autovalor(5))*real(autovalor(6)))<-1e-7 % sela hiperbolica instavel
      figure(3)
      hold on
      plot(knl1,knl2,'c*') %instavel
    end
 end

 if  imag(autovalor(5)) == -imag(autovalor(6)) & abs(imag(autovalor(5)))>1e-7 & abs(imag(autovalor(6)))>1e-7 % complexo conjugado, complexos
  
  if abs(real(autovalor(5)))>1e-7 & abs(real(autovalor(6)))>1e-7 % re(v1,v2!=0)
      
      if real(autovalor(5))>1e-7 & real(autovalor(6))>1e-7
          figure(3)
          hold on
          plot(knl1,knl2,'b*') %instavel - foco hiperbolico
      end
      if real(autovalor(5))<-1e-7 & real(autovalor(6))<-1e-7
          figure(3)
          hold on
          plot(knl1,knl2,'bo') %estavel - foco hiperbolico
     end
  end
  
  if abs(real(autovalor(5)))<1e-7 & abs(real(autovalor(6)))<1e-7
      figure(3)  
      hold on
      plot(knl1,knl2,'ko') %estavel - centro eliptico 
  end
 end

  if abs(autovalor(5)*autovalor(6))<1e-7 %v1*v2=0 Nulo - caso degenerados
      figure(3)
      hold on
      plot(knl1,knl2,'rs') %estavel
  end

  if (autovalor(5) == autovalor(6) & abs(real(autovalor(5)))>1e-7 & abs(real(autovalor(6)))>1e-7 & abs(imag(autovalor(5)))<1e-7 & abs(imag(autovalor(6)))<1e-7) % v1=v2, re(v1,v2)>0, eR
      if real(autovalor(5))>1e-7 % Nó impróprio - instavel
        figure(3)
        hold on
        plot(knl1,k21,'y*') %instavel 
     end 
     if real(autovalor(5))<-1e-7  % Nó impróprio - assintoticamente estavel
        figure(3)
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

figure(3)
title('Diagrama da Estabilidade da Fonte Não Ideal');
xlabel('Rigidez k_{nl1}')
ylabel('Rigidez k_{nl2}')


