function c=func_PCE_cijk(n,d,varargin)
%       Returns expectation value of product of multi-dimensional polynomials from Wiener-Askey scheme
%
%       Syntax: c = func_PCE_cijk(n, d, 'Hermite');
%               c = func_PCE_cijk(n, d, 'Legendre');
%               c = func_PCE_cijk(n, d, 'Laguerre', a);
%               c = func_PCE_cijk(n, d, 'Chebyshev');
%               c = func_PCE_cijk(n, d, 'Jacobi', a, b);
%
%       Input variables
%               n       indices (integer vector)
%               d       dimension (integer)
%		a       Laguerre parameter (real scalar > -1)
%               a, b    Jacobi parameters (real scalars > -1)
%
%       Output variables
%               c       (real scalar)
%
%       Description
%               c=E{H_{n(1)}*...*H_{n(end)}}
%               Hermite orthogonal w.r.t.       (2\pi)^{-1/2}\exp(-x^{2}/2)     on \real
%               Legendre orthogonal w.r.t.      1                               on ]-1,1[
%               Laguerre orthogonal w.r.t.      x^{a}\exp(-x)                   on ]0,+\infty[
%               Chebyshev orthogonal w.r.t.     (1-x^{2})^{-1/2}                on ]-1,1[
%               Jacobi orthogonal w.r.t.        (1-x)^{a}(1+x)^{b}              on ]-1,1[
%
%	See also
%		func_PCE_1DPoly	func_PCE_MultiIndex

%       Reference
%               C. Soize and R. Ghanem. Physical systems with random uncertainties: 
%               chaos representations with arbitrary probability measure. SIAM Journal 
%               on Scientific Computing, 26:395$-1ту410, 2004.
%
%               R. Ghanem and P. Spanos. Stochastic Finite Elements: A Spectral Approach. Springer, 1991.

%       Maarten Arnst, 01/21/2009
%       maarten.arnst@ulg.ac.be

% Obtain multiindices 
MultiIndices=zeros(length(n),d);
for m1=1:length(n)
    MultiIndices(m1,:)=func_PCE_MultiIndex(n(m1), d);
end

if strcmp(varargin{1},'Hermite')
    % Compute moments of standard Gaussian by recursive algorithm
    MaxOrder=max(sum(MultiIndices,1));
    moments=zeros(max(2,MaxOrder)+1,1);
    moments(1)=1;moments(2)=0;moments(3)=1;
    for m1=4:length(moments)
        moments(m1)=(m1-2)*moments(m1-2);
    end

elseif strcmp(varargin{1},'Legendre')
    % Compute moments of standard Gaussian by recursive algorithm
    MaxOrder=max(sum(MultiIndices,1));
    moments=zeros(max(2,MaxOrder)+1,1);
    moments(1)=2;
    for m1=2:length(moments)
        if mod(m1-1,2)==0 % even
            moments(m1)=2/(m1);
        else % odd
            moments(m1)=0;
        end
    end
end

% Compute cijk value
c=1;
for m1=1:d
    p=1;
    for m2=1:length(n)
        p=conv(p,func_PCE_1DPoly(MultiIndices(m2,m1),varargin{1}));
    end
    c=c*([moments(length(p):-1:1)]'*p');
end

