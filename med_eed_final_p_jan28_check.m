% ----- IED KMM -------
clear all;
Pf = 0.1;
% Pf = 0.5;
SNR_dB= -5; % SNR in dB
gam=power(10,SNR_dB/10); % SNR ratio
N=100; % No of samples

m =0.5:0.01:5;
 
for i = 1:length(m) 
    p = m(i);
    mu0 = ((2^(p/2))*gamma((p+1)/2))/sqrt(pi);
    var0 = (((2^p)*(gamma(((2*p)+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi))))/sqrt(pi));
    mu1 = ((2^(p/2))*((1+gam)^(p/2))*gamma((p+1)/2))/sqrt(pi);
    var1 = (((2^p)*((1+gam)^p)*(gamma(((2*p)+1)/2)-((gamma((p+1)/2)^2)/sqrt(pi)))))/sqrt(pi);
    lamda = ((qfuncinv(Pf)*sqrt(var0))/sqrt(N))+mu0;
    Pd(i) = qfunc(sqrt(N)*(lamda-mu1)/sqrt(var1));     
end

% [v1 p1]=max(Pd);
% opt_m1 = m(p1)

plot(m,Pd,'b','linewidth',2); hold on
axis([0.5 5 0.4 1])
% plot(opt_m1,v1,'*k','markersize',10);hold on;
grid on