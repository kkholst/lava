function [y,x] = bspline_basis(j,n,t,x)

% Copyright (c) 2010-2012, Levente Hunyadi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% B-spline basis function value B(j,n) at x.
%
% Input arguments:
% j:
%    interval index, 0 =< j < numel(t)-n
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x (optional):
%    value where the basis function is to be evaluated
%
% Output arguments:
% y:
%    B-spline basis function value, nonzero for a knot span of n

% Copyright 2010 Levente Hunyadi

validateattributes(j, {'numeric'}, {'nonnegative','integer','scalar'});
validateattributes(n, {'numeric'}, {'positive','integer','scalar'});
validateattributes(t, {'numeric'}, {'real','vector'});
assert(all( t(2:end)-t(1:end-1) >= 0 ), ...
    'Knot vector values should be nondecreasing.');
if nargin < 4
    x = linspace(t(n), t(end-n+1), 100);  % allocate points uniformly
else
    validateattributes(x, {'numeric'}, {'real','vector'});
end
assert(0 <= j && j < numel(t)-n, ...
    'Invalid interval index j = %d, expected 0 =< j < %d (0 =< j < numel(t)-n).', j, numel(t)-n);

y = bspline_basis_recurrence(j,n,t,x);

function y = bspline_basis_recurrence(j,n,t,x)

y = zeros(size(x));
if n > 1
    b = bspline_basis(j,n-1,t,x);
    dn = x - t(j+1);
    dd = t(j+n) - t(j+1);
    if dd ~= 0  % indeterminate forms 0/0 are deemed to be zero
        y = y + b.*(dn./dd);
    end
    b = bspline_basis(j+1,n-1,t,x);
    dn = t(j+n+1) - x;
    dd = t(j+n+1) - t(j+1+1);
    if dd ~= 0
        y = y + b.*(dn./dd);
    end
elseif t(j+2) < t(end)  % treat last element of knot vector as a special case
	y(t(j+1) <= x & x < t(j+2)) = 1;
else
	y(t(j+1) <= x) = 1;
end