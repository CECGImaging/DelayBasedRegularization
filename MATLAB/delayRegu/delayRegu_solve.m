%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function X = delayRegu_solve(B, A, L, pairs, actTimes, lambdaL, lambdaD, X0)
    
    N = size(A,2);
    T = size(B,2);
    
    if nargin < 8
        X0 = zeros(N,T);
    end
    
    AA = A'*A;
    LL = L'*L;
    Ab = reshape(A'*B, [], 1);
    
    D = delayRegu_assembleDifferenceMatrix(pairs, actTimes, T);
    DD = D'*D;
    
    AAinv = inv(AA + lambdaL*LL + lambdaD*speye(N));
    
    pcg_count = 0;
    msg = '';
    
    x0 = reshape(X0, T*N, 1);
    [x,flag] = pcg_quiet(@afun, Ab, 1e-6, 2000, @mfun, [], x0);
    
    fprintf(repmat('\b', 1, numel(msg)));
    fprintf('%i pcg iterations', pcg_count);
    if flag
        fprintf('. PCG FAILED WITH FLAG %i', flag);
    end
    
    X = reshape(x,N,T);
    
    function y = afun(x)
        y = reshape((AA+lambdaL*LL) * reshape(x,N,T), T*N, 1);
        y = y + lambdaD*DD*x;
        if ~mod(pcg_count,10)
            fprintf(repmat('\b', 1, numel(msg)));
            msg = sprintf('%i ', pcg_count);
            fprintf('%s', msg);
        end
        pcg_count = pcg_count+1;
    end
    
    function y = mfun(x)
        y = reshape(AAinv * reshape(x,N,T), T*N, 1);
    end
    
end

function [x,flag,relres,iter,resvec] = pcg_quiet(varargin)
    [x,flag,relres,iter,resvec] = pcg(varargin{:});
end