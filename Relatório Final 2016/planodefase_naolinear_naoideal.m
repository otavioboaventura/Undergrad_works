%PLANO DE FASE 
clear all
close all
clc
format short

ti=0;                   % tempo inicial de integração
tf=200;                 % tempo final de integração
ptos=4000;              % numero de pontos desejados
passo=(tf-ti)/ptos;     % passo de integração
lixo= 3000;             % pontos de regime transiente

[t,y]=ode45(@sistema_naolinear_naoideal,[ti:passo:tf],[0.1 0 0 0 2 0]);
   kk1=t;      % Tempo
   kk2=y(:,1); % Posição de m1
   kk3=y(:,2); % Velocidade de m1
   kk4=y(:,3); % Posição de m2
   kk5=y(:,4); % Velocidade de m2
   kk6=y(:,5); % Posição do pendulo
   kk7=y(:,6); % Velocidade do pendulo
   
%---- Históricos no tempo ----%
    figure(1)
      subplot(3,2,1);
      plot(kk1,kk2,'r');
      title('(a)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Posição [taxa]','fontsize',10); %x1
      set (gcf, 'Color', 'white');
      box on;
     subplot(3,2,2); 
      plot(kk1,kk3,'r');
      title('(b)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
     subplot(3,2,3);
      plot(kk1,kk4,'b');
      title('(c)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Posição [taxa]','fontsize',10); %x1
      set (gcf, 'Color', 'white');
      box on;
    subplot(3,2,4);
      plot(kk1,kk5,'b');
      title('(d)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
    subplot(3,2,5);     
      plot(kk1,kk6,'k');
      title('(e)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Posição [taxa]','fontsize',10); %x1
      set (gcf, 'Color', 'white');
      box on;
    subplot(3,2,6);  
      plot(kk1,kk7,'k');
      title('(f)','fontsize',10);
      xlabel('Tempo [taxa]','fontsize',10);
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(1),'-deps','history1');
      
 %---- Retratos de Fase -----%
     figure(2) 
      subplot(3,1,1);
      plot(kk2,kk3,'r');
      title('(a)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      
     subplot(3,1,2);
      plot(kk4,kk5,'b');
      title('(b)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(2),'-deps','plan1');
      
     subplot(3,1,3);
      plot(kk6,kk7,'k');
      title('(c)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(2),'-deps','plan1');
     
      
%--- Removendo o regime transiente -----%
     figure (3)
      subplot(3,1,1);
      plot(kk2(lixo:4000),kk3(lixo:4000),'r')
      title('(a)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      
     subplot(3,1,2);
      plot(kk4(lixo:4000),kk5(lixo:4000),'b')
      title('(b)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(3),'-deps','plan1');
      
     subplot(3,1,3);
      plot(kk6(lixo:4000),kk7(lixo:4000),'k')
      title('(c)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(3),'-deps','plan1');
     
      
%---- Espectros de frequencia -----%

% N=256;
% X=abs(fft(k3,N));
% X=fftshift(X)/100;
% F=[-N/2:N/2-1]/N*10;
% figure(11);
% plot(F,X,'k');
% %xlim([0. 3]);
% xlabel('Frequencia','fontsize',14);
% ylabel('Deslocamento (x_2)','fontsize',14);
% title('Espectro de Frequencia','fontsize',14)
% 
% N=256;
% X=abs(fft(k5,N));
% X=fftshift(X)/1000;
% F=[-N/2:N/2-1]/N*10;
% figure(12);
% plot(F,X,'k');
% %xlim([0. 3]);
% xlabel('Frequencia','fontsize',14);
% ylabel('Velocidade (x_4)','fontsize',14);
% title('Espectro de Frequencia','fontsize',14)
% 
% 
% N=256;
% X=abs(fft(k7,N));
% X=fftshift(X)/1000;
% F=[-N/2:N/2-1]/N*10;
% figure(13);
% plot(F,X,'k');
% %xlim([0. 3]);
% xlabel('Frequencia','fontsize',14);
% ylabel('Velocidade (x_6)','fontsize',14);
% title('Espectro de Frequencia','fontsize',14)
% 
% %print -deps figure(1)
% 
% % Cláudio Basqueroto
% % x0=[0.1 0 0 0 0 0]'; % condição inicial
% % w=0.5;
% % option=odeset ('RelTol', 1e-8, 'AbsTol',1e-8); % argumento de entrada
% % % Poincaré%
% % [t y]=ode45('sistema_naolinear_naoideal',-(20000/w)*pi:(2/w)*pi:(20000/w)*pi,x0,option);
% % plot(x(:,1),x(:,2),'.','MarkerSize',6)
% break
% 
% clear all; close all; clc; format long;
% %%%%%%%% initial conditions %%%%%%%%
% t0=0;
% tf=1000;
% x0=[0.1 0 0 0 0.5 0]' ;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% w=0.5; % frequencia de entrada
% Te=2*pi/w; % cycle time  -- perioodo do ciclo limite.
% RES=1000; % resolution - number of sample points per cycle
% Tspan=[t0:(Te/RES):tf];
% option=odeset ('RelTol', 1e-8, 'AbsTol',1e-8); % option argument can be omitted
% [t,X]=ode45 (@sistema_naolinear_naoideal, Tspan, x0, option);
% m=length(t); % tamanho total do tempo


% %----- POINCARE SAMPLING ------%
% per=0.5; % inicio do estado estacionario ---beginning of steady state
% poincare_x=X(round(m*per):RES:m , 1) ;
% poincare_y=X(round(m*per):RES:m , 2) ;
% figure (1);
% subplot (3,1,1) , plot (t,X(:,1));
% xlabel ('t');ylabel ('x1');grid;
% subplot (3,1,2) , plot (X(round(m*per):m,1),X(round(m*per):m,2));
% xlabel ('x1');ylabel ('x2');grid;
% subplot (3,1,3) , plot (poincare_x, poincare_y,'+');
% xlabel ('x1');ylabel ('x2');grid;



