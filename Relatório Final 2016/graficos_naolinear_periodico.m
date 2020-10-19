%====================== GRÁFICOS ==========================%
clear all
close all
clc
format short

ti=0;                   % tempo inicial de integração
tf=200;                 % tempo final de integração
ptos=8000;              % numero de pontos desejados
passo=(tf-ti)/ptos;     % passo de integração
lixo= 1000;             % pontos de regime transiente

[t,y]=ode45(@sistema_naolinear_periodico,[ti:passo:tf],[0.1 0 0 0]);
   kk1=t;
   kk2=y(:,1); % Posição de m1
   kk3=y(:,2); % Velocidade de m1
   kk4=y(:,3); % Posição de m2
   kk5=y(:,4); % Velocidade de m2
   
%---- Históricos no tempo ----%
   figure(1)   
     subplot(2,2,1); % Histórico no tempo da posição de m1
      plot(kk1,kk2,'r');
      title('a)(t,x1)','fontsize',10);
      xlabel('Tempo [s]','fontsize',10);
      ylabel('Posição [m]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
     subplot(2,2,2); 
      plot(kk1,kk3,'r'); % Histórico no tempo da velocidade de m1
      title('b)(t,x1´) ','fontsize',10);
      xlabel('Tempo [s]','fontsize',10);
      ylabel('Velocidade [m/s]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
     subplot(2,2,3); % Histórico no tempo da posição de m2
      plot(kk1,kk4,'b');
      title('c)(t,x2)','fontsize',10);
      xlabel('Tempo [s]','fontsize',10);
      ylabel('Posição [m]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
    subplot(2,2,4);
      plot(kk1,kk5,'b'); % Histórico no tempo da velocidade de m2
      title('d)(t,x2´)','fontsize',10);
      xlabel('Tempo [s]','fontsize',10);
      ylabel('Velocidade [m/s]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
    
 %---- Retratos de Fase -----%
      figure(2) 
     subplot(2,1,1);
      plot(kk2,kk3,'r'); % massa 1 
      title('(a)Estrutura','fontsize',10);
      xlabel('Posição [rad]','fontsize',10);
      ylabel('Velocidade [m/s]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
     subplot(2,1,2); % massa 2
      plot(kk4,kk5,'b');
      title('(b)TMD','fontsize',10);
      xlabel('Posição [rad]','fontsize',10);
      ylabel('Velocidade [m/s]','fontsize',10);
      set (gcf, 'Color', 'white');
      box on;
      print (figure(2),'-deps','plan1');
      
      
%--- Removendo o regime transiente -----%
     figure (3)
      subplot(2,1,1);
      plot(kk2(lixo:4000),kk3(lixo:4000),'r')
      title('(a)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      
     subplot(2,1,2);
      plot(kk4(lixo:4000),kk5(lixo:4000),'b')
      title('(b)','fontsize',10);
      xlabel('Posição [taxa]','fontsize',10); %x1
      ylabel('Velocidade [taxa]','fontsize',10); %x2
      set (gcf, 'Color', 'white');
      box on;
      print (figure(3),'-deps','plan1');
      

