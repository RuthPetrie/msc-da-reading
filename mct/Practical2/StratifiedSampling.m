%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StratifiedSampling.m
%
%
% Script to compute the integral 
%
%  \theta = \int_0^1 e^{x^2} dx
%
% by stratified sampling.
%
% Hard-coded to divide [0,1] into 10 intervals with 100 samples each
%
% 
% S. L. Dance January 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% loop over 10 intervals and generate samples
for k=1:10
    X(:,k) = rand(100,1)*0.1 +(k-1.0e0)/10.0e0;
end 

% evaluate integrand at sample points
Y = exp(-X.^2);

% compute Monte Carlo estimate
theta= 0.0010e0*sum(sum(Y))

% compute variance
vartheta = sum(std(X).^2)

