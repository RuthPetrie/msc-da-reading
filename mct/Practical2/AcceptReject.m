%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcceptReject.m
% 
% A simple script to sample using the Accept-Reject method
% with q(x) ~ unif[0,1]
%
% function [Y, M] = AcceptReject(N,p, c, seed)  
%
% INPUT ARGUMENTS
% N     - integer value giving number of samples required
% p     - inline function describing the pdf for the distribution 
%         we wish to sample from 
% c     - constant c such that p(x) <= c q(x) for all x
%         (in this simple version q(x) = 1 (the unif[0,1] pdf ) 
% seed  - optional argument - scalar integer seed for the 
%         random number generator to enable repetition of experiments
%
% OUTPUT ARGUMENTS
% X     - N-vector containing samples from the required distribution
% M     - number of iterations (total number of samples generated, 
%          including both acceptances and rejections)
%
% Example of use with seed (just omit this argument if not needed):
%
% 1. create a function to compute the pdf  - e.g. type  
% 
% function s=pfun(t)
% s = exp(-t);
% 
% in a file and save as pfun.m in your working directory.  
%
% 2. Then call the accept-reject method at the command prompt
% >> [Y, M]=AcceptReject(100,@pfun,1.0e0, 1);
%
% Note that the function call uses the @symbol to refer to the function 
% handle.
%
% S L Dance January 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, M]=AcceptReject(N, p,c, seed)

% set up random generator seed if applicable
if(nargin==4)
    s=seed;
    rand('twister', s);
end 

% initialize variables
i=1;
M=1;
cinv=1.0e0/c;

% accept-reject loop
while i <= N
    U=rand; X=rand;
    if(U < cinv*p(X))
        Y(i) = X;
        i=i+1;
    end;
    M=M+1;
end;

return;


