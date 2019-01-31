% ----- IED KMM -------
clear all;
Pf = 0:0.1:1;
% Pf = 0.5;
SNR_dB= -5; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=100; % No of samples

% avg_noise_pow = 1;
f=6e6; % carrier freq for modn
bw = 6e6; % bandwidth of interest
noise_figure = 11; % noise figure in dB (as specified in IEEE 802.22)
noise_power = -174+10*log(bw)+noise_figure; %noise power in dBm conversion
avg_noise_pow=double((1e-3)*power(10,noise_power/10));% received power in watts

%----- uniformly distributed random numbers-----

m = [0.5000 2.1300  1.8700  1.7000  1.5700  1.4600  1.3600    1.2500    1.1400    0.9900    0.5000];

% mu0 = (2^(p/2))*gamma((p+1)/2)*(noise_pow^(p/2))/sqrt(pi);
% var0 = (2^p)*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))*avg_noise_pow^p/sqrt(pi);
% mu1 = (2^(p/2))*((1+gam)^(p/2))*gamma((p+1)/2)*(noise_pow^(p/2))/sqrt(pi);
% var1 = (2^p)*((1+gam)^(p))*(gamma((2*p+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))*avg_noise_pow^p/sqrt(pi);

M = 10;
L = 10;
 
for i = 1:length(Pf) 
    p = m(i);
    mu0 = N*((2^(p/2))*gamma((p+1)/2)*(avg_noise_pow^(p/2))/sqrt(pi));
    var0 = N*(((2^p)*(gamma(((2*p)+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))*(avg_noise_pow^p))/sqrt(pi));
    mu1 = N*((2^(p/2))*((1+gam)^(p/2))*gamma((p+1)/2)*(avg_noise_pow^(p/2))/sqrt(pi));
    var1 = N*(((2^p)*((1+gam)^p)*(gamma(((2*p)+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))*(avg_noise_pow^p))/sqrt(pi));
    mu_avg = ((M/L)*mu1)+ (((L-M)/L)*mu0);
    sig_avg =((M/L^2)*var1)+ (((L-M)/L^2)*var0) ;
    lamda(i) = ((qfuncinv(Pf(i)))*(sqrt(var0)))+ mu0;
    Pd(i) = qfunc((lamda(i)-mu1)/sqrt(var1));
    Pdied(i) = Pd(i)+(Pd(i)*((1-Pd(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
    Pfied(i) = Pf(i)+(Pf(i)*((1-Pf(i))*qfunc((lamda(i)-mu_avg)/sqrt(sig_avg))));
       
end
% plot(Pfc,Pdc,'b','linewidth',2); hold on
% plot(Pfa,Pdt,'b','linewidth',2);hold on
plot(Pfied,Pdied,'b','linewidth',1); hold on
% plot(Pf,Pd,'r');hold on;