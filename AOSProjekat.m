%% Procena nepoznatog sistema
close all;
clear all;
clc;
fs=1000;
t=1/fs:1/fs:1;
%sig=sin(150*2*pi*t) + 0.1*randn(size(t));
N_usrednji=200;
sigmav_2=1;
N_eksperimenta=10000;
mi1=0.001;%const
mi2=0.001;%Lorentzian
mi2p=0.001;
dMi2=zeros(N_eksperimenta,N_usrednji);
mi3=0.001;%novi Lorentzian
mi3p=0.001;
dMi3=zeros(N_eksperimenta,N_usrednji);
alfa=0.0015;
delta=0.05;
M=100; % duzina adaptivnog filtra

%% IIR nepoznat sistem
[b,a]=ellip(5,1,40,[0.2 0.4]);
%[z,p,k] = ellip(6,3,50,10/500,'low');
%sos = zp2sos(z,p,k);
[H,f]=freqz(b,a,1000);
%[H,f]=freqz(sos,1000,"whole");
J1=zeros(N_eksperimenta,N_usrednji);
J2=zeros(N_eksperimenta,N_usrednji);
J3=zeros(N_eksperimenta,N_usrednji);
w_pamti_11=zeros(N_eksperimenta,M);
w_pamti_12=zeros(N_eksperimenta,M);
w_pamti_13=zeros(N_eksperimenta,M);
w_pamti1=zeros(N_eksperimenta,N_usrednji,M);
w_pamti2=zeros(N_eksperimenta,N_usrednji,M);
w_pamti3=zeros(N_eksperimenta,N_usrednji,M);
%% Racunanje procena, odstupanja i usrednjavanje
for br_usrednji=1:N_usrednji
    v=sqrt(sigmav_2)*randn(N_eksperimenta,1);
    d=filter(b,a,v);
    x=v+sqrt(sigmav_2/100)*randn(N_eksperimenta,1);
    w1=zeros(1,M);
    w2=zeros(1,M);
    w3=zeros(1,M);
    xn=zeros(M,1);
    e3p=d(1)-w3*[x(1);xn(1:end-1)];
    for br=1:N_eksperimenta
        xn=[x(br);xn(1:end-1)];
        d_procena1=w1*xn;
        d_procena2=w2*xn;
        d_procena3=w3*xn;
        e1=d(br)-d_procena1;
        e2=d(br)-d_procena2;
        e3=d(br)-d_procena3;

        mi2=alfa* log(1+ ((e2/delta).^2)/2);
%         if (br==1000) && (br_usrednji==1) %glitch
%             mi2=1;%do 1 je stabilno
%         end
         dMi2(br,br_usrednji)=abs(mi2-mi2p);
         mi3=alfa* log(1+ ((abs(e3*e3p)/((delta).^2)*2)));
%         if (br==4000) && (br_usrednji==1) %glitch
%             mi3=1;%do 1 je stabilno
%         end
        dMi3(br,br_usrednji)=abs(mi3-mi3p);
        J1(br,br_usrednji)=e1.^2;
        J2(br,br_usrednji)=e2.^2;
        J3(br,br_usrednji)=e3.^2;
        w1=w1+mi1*xn'*e1;
        w2=w2+mi2*xn'*e2;
        w3=w3+mi3*xn'*e3;
        w_pamti_11(br,:)=sign(real(w1)).*abs(w1); % za jedan slučajni signal
        w_pamti1(br,br_usrednji,:)=sign(real(w1)).*abs(w1); % usrednjeno
        w_pamti_12(br,:)=sign(real(w2)).*abs(w2);
        w_pamti2(br,br_usrednji,:)=sign(real(w2)).*abs(w2);
        w_pamti_13(br,:)=sign(real(w3)).*abs(w3);
        w_pamti3(br,br_usrednji,:)=sign(real(w3)).*abs(w3);
        e3p=e3;
        mi2p=mi2;
        mi3p=mi3;
    end
end
%% Iscrtavanje
figure,subplot(2,1,1);semilogy(1:N_eksperimenta,dMi2(:,1));
xlabel('n');
ylabel('\mu_2[n]');
ylim([10^(-6) 10^(-2)]);
legend('\mu_2 var');
title('Promena \mu_2 i \mu_3');
subplot(2,1,2);semilogy(1:N_eksperimenta,dMi3);
xlabel('n');
ylabel('\mu_3[n]');
ylim([10^(-6) 10^(-2)]);
legend('\mu_3 var');

%% Cost funkcija i koeficijenti filtra
figure,subplot(3,1,1);semilogy(1:N_eksperimenta,J1);
xlabel('n');
ylabel('J[n]');
ylim([10^(-3) 10]);
legend(['\mu_1=' num2str(mi1)]);
title('Cost funkcija za svaki eksperiment');
subplot(3,1,2);semilogy(1:N_eksperimenta,J2);
xlabel('n');
ylabel('J[n]');
ylim([10^(-3) 10]);
legend('\mu_2 var');
subplot(3,1,3);semilogy(1:N_eksperimenta,J3);
xlabel('n');
ylabel('J[n]');
ylim([10^(-3) 10]);
legend('\mu_3 var');

figure,subplot(3,1,1);plot(w_pamti_11);
xlabel('n');
ylabel('w[n]');
legend(['\mu_1=' num2str(mi1)]);
title('Koeficijenti filtra za jedan eksperiment');
subplot(3,1,2);plot(w_pamti_12);
xlabel('n');
ylabel('w[n]');
legend('\mu_2 var');
subplot(3,1,3);plot(w_pamti_13);
xlabel('n');
ylabel('w[n]');
legend('\mu_3 var');

w_usrednjeno1(1:N_eksperimenta,1:M)=squeeze(sum(w_pamti1,2)/N_usrednji);

figure,subplot(3,1,1);plot(1:N_eksperimenta,w_usrednjeno1);
xlabel('n');
ylabel('w_k[n]');
legend(['\mu_1=' num2str(mi1)]);
title('Koeficijenti filtra w_k[n] usrednjeni');
w_usrednjeno2(1:N_eksperimenta,1:M)=squeeze(sum(w_pamti2,2)/N_usrednji);
subplot(3,1,2);plot(1:N_eksperimenta,w_usrednjeno2);
xlabel('n');
ylabel('w_k[n]');
legend('\mu_2 var');
w_usrednjeno3(1:N_eksperimenta,1:M)=squeeze(sum(w_pamti3,2)/N_usrednji);
subplot(3,1,3);plot(1:N_eksperimenta,w_usrednjeno3);
xlabel('n');
ylabel('w_k[n]');
legend('\mu_3 var');

figure,subplot(3,1,1);plot(1:N_eksperimenta,w_pamti1(:,:,1)); hold on;
plot(1:N_eksperimenta,w_usrednjeno1(:,1),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_0[n]');
legend(['\mu_1=' num2str(mi1)]);
title('Koeficijent filtra w_0[n] za svaki eksperiment');
subplot(3,1,2);plot(1:N_eksperimenta,w_pamti2(:,:,1)); hold on;
plot(1:N_eksperimenta,w_usrednjeno2(:,1),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_0[n]');
legend('\mu_2 var');
subplot(3,1,3);plot(1:N_eksperimenta,w_pamti3(:,:,1)); hold on;
plot(1:N_eksperimenta,w_usrednjeno3(:,1),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_0[n]');
legend('\mu_3 var');

k0=10;

figure,subplot(3,1,1);plot(1:N_eksperimenta,w_pamti1(:,:,k0)); hold on
plot(1:N_eksperimenta,w_usrednjeno1(:,k0),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_{k0}[n]');
legend(['\mu_1=' num2str(mi1)]);
title(['Koeficijent filtra w_', num2str(k0-1), '[n] za svaki eksperiment']);
subplot(3,1,2);plot(1:N_eksperimenta,w_pamti2(:,:,k0)); hold on
plot(1:N_eksperimenta,w_usrednjeno2(:,k0),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_{k0}[n]');
legend('\mu_2 var');
subplot(3,1,3);plot(1:N_eksperimenta,w_pamti3(:,:,k0)); hold on
plot(1:N_eksperimenta,w_usrednjeno3(:,k0),'k--','linewidt',1.5);
xlabel('n');
ylabel('w_{k0}[n]');
legend('\mu_3 var');

figure,semilogy(1:N_eksperimenta,mean(J1,2),1:N_eksperimenta,mean(J2,2),1:N_eksperimenta,mean(J3,2));
xlabel('n');
ylabel('J[n]');
title('Usrednjena Cost funkcija');
legend(['\mu_1=' num2str(mi1)],'\mu_2 var','\mu_3 var');
ylim([10^(-3) 1]);

w_kraj1(1:N_usrednji,1:M)=w_pamti1(end,:,:);
w_kraj_sr1=mean(w_kraj1);
[W1,f]=freqz(w_kraj_sr1,1,1000);
w_kraj2(1:N_usrednji,1:M)=w_pamti2(end,:,:);
w_kraj_sr2=mean(w_kraj2);
[W2,f]=freqz(w_kraj_sr2,1,1000);
w_kraj3(1:N_usrednji,1:M)=w_pamti3(end,:,:);
w_kraj_sr3=mean(w_kraj3);
[W3,f]=freqz(w_kraj_sr3,1,1000);

%% Procena sistema
figure,plot(f/pi,20*log10(abs(H)),f/pi,20*log10(abs(W1)),f/pi,20*log10(abs(W2)),f/pi,20*log10(abs(W3)));
xlabel('\omega/\pi');
ylabel('M(\omega)');
legend('nepoznat',['Procenjen \mu_1=' num2str(mi1)],'Procenjen \mu_2 var','Procenjen \mu_3 var')

figure,subplot(2,1,1),stem(impz(b,a));
xlim([0 100]);
xlabel('n');
ylabel('h [n]');
legend('nepoznati sistem');
subplot(2,1,2),stem(0:M-1,[w_kraj_sr1' w_kraj_sr2' w_kraj_sr3']);
xlim([0 100]);
xlabel('n');
ylabel('w [end]');
legend(['\mu_1=' num2str(mi1)], '\mu_2 var', '\mu_3 var');