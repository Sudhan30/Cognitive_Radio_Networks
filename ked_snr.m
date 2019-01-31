% ----- IED KMM -------
clear all;
Pf = 0.001;
SNR_dB= -15:1:5; % SNR in dB 
gam=power(10,SNR_dB/10); % SNR ratio
N=100; % No of samples
nu = 0.1;
ro = power(10,nu/10);

avg_noise_pow = 1;
%----- uniformly distributed random numbers-----
a = avg_noise_pow/ro;
b = avg_noise_pow*ro;
% r = a + (b-a).*rand(1,1);
% act_noise_pow = r;
% Pd calculation using formula
% for i = 1:length(Pf) 
%     lamda(i) = avg_noise_pow*(N +((sqrt(2*N))*qfuncinv(Pf(i))));
%     v1 = sqrt((2*N*(gam+1)^2*avg_noise_pow^2));
%     v2 = lamda(i) -  (N*(1+gam)*avg_noise_pow);
%     Pd(i) = qfunc(v2/v1);
% end
run = 1000;
for i = 1:length(SNR_dB) 
    gam=power(10,SNR_dB(i)/10);
    lamda = avg_noise_pow*(N +((sqrt(2*N))*qfuncinv(Pf)));
    for j = 1:run
        r = a + (b-a).*rand(1,1);
        act_noise_pow = r;
        v1 = sqrt((2*N*(gam+1)^2*act_noise_pow^2));
        v2 = lamda -  (N*(1+gam)*act_noise_pow);
        Pdc1(j) = qfunc(v2/v1);
        v3 = sqrt((2*N*act_noise_pow^2));
        v4 = lamda -  (N*act_noise_pow);
        Pfc1(j) = qfunc(v4/v3);
    end
    Pdc(i) = mean(Pdc1);
    Pfc(i) = mean(Pfc1);
end

Pe = Pfc+1-Pdc;
plot(SNR_dB,Pe); hold on;


M = 10;
L = 10;
mu_avg = ((M/L)*N*(1+gam)*avg_noise_pow)+ (((L-M)/L)*N*avg_noise_pow);
sig_avg = ((M/(L^2))* (2*N)*(((1+gam)^2)*(avg_noise_pow^2))+(((L-M)/(L^2))*(2*N)*(avg_noise_pow^2)));

for i = 1:length(SNR_dB) 
    gam=power(10,SNR_dB(i)/10);
    lamda1 = ((((qfuncinv(Pf))*(sqrt(2*N)))+N)*(avg_noise_pow/ro));
    lamda2 = ((((qfuncinv(Pf))*(sqrt(2*N)))+N)*(avg_noise_pow*ro));
    for j = 1:run
        r = a + (b-a).*rand(1,1);
        act_noise_pow = r;
        v1 = sqrt((2*N*(gam+1)^2*act_noise_pow^2));
        v2 = lamda2 - (N*(1+gam)*act_noise_pow);
        Pd1(j) = qfunc(v2/v1);
        v3 = sqrt((2*N*act_noise_pow^2));
        v4 = lamda2 -  (N*act_noise_pow);
        Pf1(j) = qfunc(v4/v3);
    end
    Pdt(i) = mean(Pd1);
    Pfa(i) = mean(Pf1);
    td = (lamda1-(N*(1+gam)*act_noise_pow))/v1;
    tf = (lamda1-(N*act_noise_pow))/sqrt(2*N*act_noise_pow^2);
    ta1 = (lamda2-mu_avg)/sqrt(sig_avg);
    ta2 = (lamda2-mu_avg)/sqrt( sig_avg);
    Pdked1(i) = Pdt(i)+((1-Pdt(i))*(qfunc(td)-Pdt(i))*qfunc(ta1));
    Pfked1(i) = Pfa(i)+((1-Pfa(i))*(qfunc(tf)-Pfa(i))*qfunc(ta2));
end

Pek = Pfked1+1-Pdked1;

plot(SNR_dB,Pek,'m');


% plot(Pfc,Pdc,'b','linewidth',2); hold on
% plot(Pf,Pd,'m','linewidth',2);hold on
% plot(Pfa,Pdt,'g','linewidth',2); hold on
% 
