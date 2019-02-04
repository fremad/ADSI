clear
clc


%% Variables

N=100;

r= 0.9;
w0 = pi/3;
K = 8;
w = linspace(0,1,1000)*pi;

%where to get b ??
b=0.0274^4;

mu = 0;
sigma = 2;
n = linspace(-5,5,N);

w_1 = 0.34*pi;
w_2 = 0.6*pi;


%% Signals



s1 = normpdf(n,mu,sigma);
n= 0:4*N-1;

xn = zeros(1,length(n));
xn(1:N) = s1.*cos(w_1*(0:N-1));
xn(N+1:2*N) = s1.*cos(w_2*(0:N-1));




bcas = b^(1/K);
acas = [1 -2*r*cos(w0) r^2];




Hcas = freqz(bcas,acas,w);

gd_cas = K*grpdelay(bcas,acas,w);

Hcas_mag_db = 10*log10(abs(Hcas).^2)*K;

yn_cas = xn;

for ii=1:K
yn_cas = filter(bcas,acas,yn_cas);
end

X = (fft(xn));


%% Plots

subplot(411)
plot(1:400,xn)

subplot(412)

plot(1:400,yn_cas)
axis([1,400,-0.2,0.2])

subplot(413)
plotyy(w/pi,Hcas_mag_db,w/pi,gd_cas)
grid on
subplot(414)

mags = abs(X);

plot((1:200)*0.005,mags(1:200))

figure 
zplane(bcas,acas)