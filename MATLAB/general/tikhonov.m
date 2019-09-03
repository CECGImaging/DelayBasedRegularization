%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [X,resNorm,solNorm] = tikhonov(B, A, R, lambda, zeroMeanNullspaceRemoval)

if nargin < 5
    zeroMeanNullspaceRemoval = false;
end

AA = A'*A;
RR = R'*R;
AB = A'*B;

if zeroMeanNullspaceRemoval
    AA = AA + repmat(mean(abs(AB(:)))/size(A,2), size(A,2));
end

if numel(lambda) == 1 % time-constant lambda
    
    X = (AA + lambda*RR) \ AB;

    if nargout > 1
        resNorm = norm(A*X-B, 'fro');
    end
    if nargout > 2
        solNorm = norm(R*X, 'fro');
    end
    
elseif numel(lambda) == size(B,2) % time-varying lambda
    
    X = NaN(size(A,2), size(B,2));
    for i = 1:size(B,2)
        X(:,i) = (AA + lambda(i)*RR) \ A'*B(:,i);
    end

    if nargout > 1
        resNorm = sqrt(sum((A*X-B).^2,1));
    end
    if nargout > 2
        resNorm = sqrt(sum((R*X).^2,1));
    end
    
else
    error('numel(lambda) is neither 1 nor size(B,2).');
end

end
