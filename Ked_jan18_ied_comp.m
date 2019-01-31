% ----- IED KMM -------
clear all;
Pf = 0:0.01:1;
% Pf = 0.5;
SNR_dB= -10; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=100; % No of samples
nu = 0;
ro = power(10,nu/10);

avg_noise_pow = 1e-6;
%----- uniformly distributed random numbers-----
a = avg_noise_pow/ro;
b = avg_noise_pow*ro;
% r = a + (b-a).*rand(1,1);
% act_noise_pow = r;
% Pd calculation using formula
run = 1000;
M = 10;
L = 10;
for k = 1:L
    rp(k) = a + (b-a).*rand(1,1);
end

past_avg_np = mean(rp); 
mu_avg = ((M/L)*N*(1+gam)*past_avg_np)+ (((L-M)/L)*N*past_avg_np);
sig_avg = ((M/(L^2))* (2*N)*(((1+gam)^2)*(past_avg_np^2))+(((L-M)/(L^2))*(2*N)*(past_avg_np^2)));

for i = 1:length(Pf) 
    lamda(i) = ((((qfuncinv(Pf(i)))*(sqrt(2*N)))+N)*(avg_noise_pow));
    for j = 1:run
        r = a + (b-a).*rand(1,1);
        act_noise_pow = r;
        v1 = sqrt((2*N*(gam+1)^2*act_noise_pow^2));
        v2 = lamda(i) - (N*(1+gam)*act_noise_pow);
        Pd1(j) = qfunc(v2/v1);
        v3 = sqrt((2*N*act_noise_pow^2));
        v4 = lamda(i) -  (N*act_noise_pow);
        Pf1(j) = qfunc(v4/v3);    
        ta = (lamda(i)-mu_avg)/sqrt(sig_avg);
        Pdked1(j) = Pd1(j)+((1-Pd1(j))*Pd1(j)*qfunc(ta));
        Pfked1(j) = Pf1(j)+((1-Pf1(j))*Pf1(j)*qfunc(ta));
    end
    Pdt(i) = mean(Pd1);
    Pfa(i) = mean(Pf1);
    Pdked(i) = mean(Pdked1);
    Pfked(i) = mean(Pfked1);
end
% plot(Pfc,Pdc,'b','linewidth',2); hold on
% plot(Pf,Pd,'m','linewidth',2);hold on
plot(Pfked,Pdked,'m','linewidth',2); hold on
