function y = func_PCE_Chaos(x, c, d, p, varargin);
%       Returns value of normalized polynomial chaos expansion from Wiener-Askey scheme
%
%       Syntax: c = func_PCE_Chaos(x, c, d, p, 'Hermite');
%               c = func_PCE_Chaos(x, c, d, p, 'Legendre');
%               c = func_PCE_Chaos(x, c, d, p, 'Laguerre', a);
%               c = func_PCE_Chaos(x, c, d, p, 'Chebyshev');
%               c = func_PCE_Chaos(x, c, d, p, 'Jacobi', a, b);
%
%       Input variables
%		x	argument (d*1 real vector)
%		c	coefficients (n*1 real vector)
%		d	dimension (integer)
%               p       order (integer)
%               a       Laguerre parameter (real scalar > -1)
%               a, b    Jacobi parameters (real scalars > -1)
%
%       Output variables
%               y       (real scalar)
%
%	Description
%		y = c(1)*H_{0}(x)+...+c(n)*H_{n-1}(x)
%		n = sum_{i=0}^{p}(i+d-1)!/i!/(d-1)!
%               Hermite orthogonal w.r.t.       (2\pi)^{-1/2}\exp(-x^{2}/2)     on \real
%               Legendre orthogonal w.r.t.      1                               on ]-1,1[
%               Laguerre orthogonal w.r.t.      x^{a}\exp(-x)                   on ]0,+\infty[
%               Chebyshev orthogonal w.r.t.     (1-x^{2})^{-1/2}                on ]-1,1[
%               Jacobi orthogonal w.r.t.        (1-x)^{a}(1+x)^{b}              on ]-1,1[ 
%
%       See also
%               func_PCE_1DPoly	func_PCE_MultiIndex 

%       Reference
%               C. Soize and R. Ghanem. Physical systems with random uncertainties: 
%               chaos representations with arbitrary probability measure. SIAM Journal 
%               on Scientific Computing, 26:395–410, 2004.
%
%               R. Ghanem and P. Spanos. Stochastic Finite Elements: A Spectral Approach. Springer, 1991.

%       Maarten Arnst, 01/21/2009
%       arnst@usc.edu

% Compute number of polynomials
num_pols = func_PCE_NumPols(d,p);

% Evaluate chaos decomposition at x
y = 0;

for m1=1:num_pols
    MultiIndex=func_PCE_MultiIndex(m1-1, d);
    PolyValue = 1;
    for m2=1:d
        PolyValue=PolyValue*polyval(func_PCE_1DPoly(MultiIndex(m2), varargin),x(m2));
    end    
    y=y+c(m1)*PolyValue;
end
