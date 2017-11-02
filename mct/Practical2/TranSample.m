%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TranSample.m
% 
% A simple script to sample using the transformation method
%
% function X = TranSample(N,Finv, seed)  
%
% INPUT ARGUMENTS
% N     - integer value giving number of samples required
% Finv  - inline function describing the inverse cdf for the distribution 
%         we wish to sample from 
% seed  - optional argument - scalar integer seed for the 
%         random number generator to enable repetition of experiments
%
% OUTPUT ARGUMENT
% X     - N-vector containing samples from the required distribution
%
% Example of use with seed (just omit argument if not needed):
%
% 1. create a function to compute Finv - e.g. type  
% 
% function x=g(u)
% x = -log(1-u);
% 
% in a file and save as g.m in your working directory. 
%
% 2. Then call the transformation method at the command prompt
%
% >> X=TranSample(100,@g,1);
%
% Note that the function call uses the @symbol to refer to the function 
% handle.
%
% S L Dance January 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=TranSample(N, Finv, seed)

% set up random generator seed if applicable
if(nargin==3)
    s=seed;
    rand('twister', s);
end 

% generate N unif[0,1] random samples
U=rand(N, 1);

%transform variables
X=Finv(U);


