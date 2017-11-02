%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ImportSamp.m
%
%
% Function to compute the integral 
%
%  \theta = \int_{-\infty}^\infty e^{-|x-a|/D} dx
%
% where a is some real number (an input argument of the function. 
%
% We compute by 
% 
% 1. Standard Monte Carlo quadrature
% 2. Importance sampling using N(0,1) as the proposal distribution
%
% Produce plots of sample locations as well as values of integral
% NB hard coded for N=1000 samples. 
% 
% function [MC, VMC, ImpS, VImpS]=ImportSamp(a, D, L)
% 
% INPUTS
% a     - real valued constant (parameter of integrand)
% D      - real valued constant (width parameter for integrand)
% L     - domain size for finite domain uniform dist. 
% 
% S. L. Dance January 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImportSamp(a, D, L)

N=100;

% compute standard MC estimate over a finite domain [-L, L] 
XMC=rand(N,1)*2.0e0*L-L; % sample
FXMC=2.0*L*exp(-abs(XMC-a*ones(N,1))/D); % integrand
MC=mean(FXMC);
VMC=(std(FXMC)).^2;

% compute importance sampling estimate 
XIS=randn(N,1); % normal(0,1) samples
FXIS=exp(-abs(XIS-a*ones(N,1))/D); % integrand
PXIS=(1.0/sqrt(2.0e0*pi))*exp(-0.5e0*(XIS.^2)); % proposal pdf

ImpS=mean(FXIS./PXIS);
VImpS=(std(FXIS./PXIS)).^2;

% plots

%set up grid to plot integrand
dx=L/50; 
x=-L:dx:L;
F=exp(-abs(x-a*ones(size(x)))/D);

figure;
subplot(3,1,1)
plot(x, F, 'k');
xlabel('x')
ylabel('exp(-|x-a|/D)')

subplot(3,1,2)
plot(XMC, 1.0/(2.0*L)*ones(size(XMC)), 'x')
xlabel('x')
ylabel('p-value')
tstr=['MC estimate ' num2str(MC) ' Variance ' num2str(VMC)]
title(tstr)
axis tight

subplot(3,1,3)
plot(XIS,PXIS, 'x')
xlabel('x')
ylabel('p-value')
tstr=['IS estimate ' num2str(ImpS) ' Variance ' num2str(VImpS)]
title(tstr)
axis tight

figure;
plot(XMC, FXMC*2.0*L, 'bx')
hold on
xlabel('X')
ylabel('F(X)/\pi(X)')
plot(XIS, FXIS./PXIS, 'ro')
legend('MC integrand/proposal', 'IS integrand/proposal')





