% ----- IED KMM -------
clear all;
Pf = 0:0.01:1;
% Pf = 0.5;
SNR_dB= -13; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=1000; % No of samples

f=6e6; % carrier freq for modn
bw = 6e6; % bandwidth of interest
noise_figure = 11; % noise figure in dB (as specified in IEEE 802.22)
noise_power = -174+10*log(bw)+noise_figure; %noise power in dBm conversion
avg_noise_pow=double((1e-3)*power(10,noise_power/10));% received power in watts

%----- uniformly distributed random numbers-----

p=2;

M = 10;
L = 10;

mu0 = N*avg_noise_pow;
var0 = 2*N*(avg_noise_pow^2);
mu1 = N*(1+gam)*avg_noise_pow;
var1 = 2*N*((1+gam)^2)*(avg_noise_pow^2);

 mu_avg = ((M/L)*mu1)+ (((L-  M)/L)*mu0);
 sig_avg =((M/L^2)*var1)+ (((L-M)/L^2)*var0) ;
 
for i = 1:length(Pf) 
    lamda(i) = ((qfuncinv(Pf(i)))*(sqrt(var0)))+ mu0;
    Pd(i) = qfunc((lamda(i)-mu1)/sqrt(var1));
    Pdied(i) = Pd(i)+(Pd(i)*((1-Pd(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
    Pfied(i) = Pf(i)+(Pf(i)*((1-Pf(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
end
% plot(Pfc,Pdc,'b','linewidth',2); hold on
% plot(Pfa,Pdt,'b','linewidth',2);hold on
plot(Pfied,Pdied,'-*b','linewidth',2); hold on
plot(Pf,Pd,'b');hold on