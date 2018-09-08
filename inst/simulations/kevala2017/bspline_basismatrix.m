function [B,x] = bspline_basismatrix(n,t,x)

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

% B-spline basis function value matrix B(n) for x.
%
% Input arguments:
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x (optional):
%    an m-dimensional vector of values where the basis function is to be
%    evaluated
%
% Output arguments:
% B:
%    a matrix of m rows and numel(t)-n columns

% Copyright 2010 Levente Hunyadi

if nargin > 2
    B = zeros(numel(x),numel(t)-n);
    for j = 0 : numel(t)-n-1
        B(:,j+1) = bspline_basis(j,n,t,x);
    end
else
    [b,x] = bspline_basis(0,n,t);
    B = zeros(numel(x),numel(t)-n);
    B(:,1) = b;
    for j = 1 : numel(t)-n-1
        B(:,j+1) = bspline_basis(j,n,t,x);
    end
end