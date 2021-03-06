clear
clc

%% Variables

syms n z k
N = 20;

%x and y signals
x = heaviside(sym(1));
y = 2*(1/3)^n*heaviside(sym(1));

%z-transforms of x and y
Y = ztrans(y);
X = ztrans(x);

%The impulse response
H = Y/X;
h = iztrans(H);

a = [2,-2];
b = [1, -1/3];

%% New assignment

x = (1/2)^n*heaviside(sym(1));
X = ztrans(x);
Y = H*X;
y(n) = iztrans(Y);


figure('Name','Filter function','Numbertitle','off')
subplot(211)
stem(0:N,y(0:N))
grid on
title('Output Response')
xlabel('samples (n)') 
ylabel('Amplitude') 

subplot(212)
impz(a,b,N)
grid on

figure
freqz(a,b,linspace(0,1,1000)*pi)


%% Plots
%impz(a,b,50)