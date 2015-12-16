function R = cor_bar(P)
%COR_BAR  Calculates approximate sample correlation matrix.
%
%   S=COR_BAR(P)
%
%   Produces an n-by-n approx correlation matrix based on
%   data of size m-by-n. n columns of different
%   random variables observed at m different times.
%   P has missing data represented by NaNs.
%
%   INPUT:  P  data matrix
%
%   OUTPUT: R  approx sample correlation matrix

[m,n]=size(P);
S=cov_bar(P);
D=diag(1./sqrt(diag(S)));

R=D*S*D;
