% Set up
echo off
format compact
format short
num = 128;
c1 = linspace(1,0,num);
c2 = linspace(1,0.2,num);
c3 = linspace(1,0.6,num);
colormap([c1' c2' c3']);
mymax = 0.2;
mylim = roundsd(mymax,1);

echo on


% Get data and generate an approximate correlation matrix

[P,n,k] = get_data

G = cor_bar(P)        

sort( eig(G) )'

pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Compute the Nearest Correlation Matrix
%

[Gout, X, iter, feval, nrmgrd, ifail] = g02aa(G);

X
iter

sort( eig(X) )'

imagesc(abs(X-G)), colorbar, axis on, set(gca, 'XTick', []), 
set(gca, 'YTick', []), caxis([0 mylim]); 
res = norm( X-G, 'fro' )
title(['g02aa: iterations: ' int2str(iter),', ||G-X||_F = ',num2str(res)] ); hold on



pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Compute the Nearest Correlation Matrix with row and column weighting
%
sort( eig( G(1:k,1:k) ) )'

W = [10,10,10,1,1,1,1,1], opt = 'B'; alpha = 0.001;

[Gout, W, X, iter, feval, nrmgrd, ifail] = g02ab(G, opt, alpha, W);

X
iter

sort( eig(X) )'

imagesc(abs(X-G)) 

res = norm( X-G, 'fro' )
title(['g02ab, alpha = 0.001: iterations: ' int2str(iter),', ||W(G-X)W||_F = ',num2str(res)] );
norm( X(1:k,1:k)-G(1:k,1:k), 'fro' )

pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Compute the Nearest Correlation Matrix with element-wise weighting
% Same alpha as before
%
H=ones(n); maxit = nag_int(500); H(1:k,1:k) = 100,

[Gout, H, X, iter, norm_p, ifail] = g02aj(G, alpha, H, 'maxit', maxit);

X
iter

imagesc(abs(X-G)) 

res = norm( X-G, 'fro' )
title(['g02aj, alpha = 0.001: iterations: ' int2str(iter),', ||H.*(G-X)||_F = ',num2str(res)] );
norm( X(1:3,1:3)-G(1:3,1:3), 'fro' )


pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Compute the Nearest Correlation Matrix with element-wise weighting
% Zero alpha
%
alpha = 0.0;

[Gout, H, X, iter, norm_p, ifail] = g02aj(G, alpha, H, 'maxit', maxit);

X
iter

imagesc(abs(X-G)) 

res = norm( X-G, 'fro' )
title(['g02aj, alpha = 0.0: iterations: ' int2str(iter),', ||H.*(G-X)||_F = ',num2str(res)] );
norm( X(1:3,1:3)-G(1:3,1:3), 'fro' )

pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% Compute a true Correlation Matrix with fixed leading block

[Gout, X, alpha, iter, eigmin, norm_p, ifail] = g02an(G, k);

X
iter
alpha

imagesc(abs(X-G)) 

res = norm( X-G, 'fro' )
title(['g02an: iterations: ' int2str(iter),', ||G-X||_F = ',num2str(res)] );
norm( X(1:3,1:3)-G(1:3,1:3), 'fro' )

echo off