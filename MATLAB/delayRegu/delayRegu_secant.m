%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function [X,s] = delayRegu_secant(s, residualNorm0, B, A, L, pairs, actTimes, lambdaL, X)

% s: history for extrapolation of logLambdaD from previous iterations
% s(1), s(2), s(3): logLambdaD from previous, pre-previous and 
%                   pre-pre-previous iterations, respectively
% s(4): 1/slope of residualNormDiff(logLambdaD) from previous iteration
% for the first iteration, s has to be initialized with NaNs

relTol = 1e-4;  % relative tolerance of residual norm difference
maxSteps = 20;  % maximum number of secant steps

l = NaN(1,2);   % history of logLambdaD for secant method
r = NaN(1,2);   % history of residualNormDiff for secant method
l(1) = -8;      % initial value of logLambdaD

converged = false;

for step = 1:maxSteps
    
    fprintf('secant step %i\t| ', step);
    
    if step == 1
        % no point (l,r) from previous steps available yet
        % -> check how many values of logLambdaD from previous iterations
        %    are available and use them for extrapolation
        switch 3-sum(isnan(s(1:3)))
            case 0
                % use initial value
                logLambdaD = l(1);
            case 1
                % use previous value
                logLambdaD = s(1);
            case 2
                % linear extrapolation
                logLambdaD = s(1) + (s(1)-s(2));
            otherwise
                % quadratic extrapolation
                logLambdaD = s(1) + (s(1)-s(2)) + (s(1)-2*s(2)+s(3))/2;
        end
    else
        % at least one point (l,r) available
        if step == 2
            % one point available
            if ~isnan(s(4))
                % use secant slope from previous iteration
                logLambdaD = l(1) - s(4) * r(1);
            else
                % no slope available, but we can move in the right direction
                logLambdaD = l(1) - sign(r(1)) * 1e-2;
            end
        else
            % two points available -> use actual secant update
            logLambdaD = l(1) - (l(1)-l(2))/(r(1)-r(2)) * r(1);
        end
        
        % limit maximum change of lambdaD to one decade
        if abs(logLambdaD-l(1)) > 1
            logLambdaD = l(1) + sign(logLambdaD-l(1));
        end
        
        % fallback
        if ~isfinite(logLambdaD)
            logLambdaD = l(1) - sign(r(1)) * 1e-4;
        end
    end
    
    fprintf('logLambdaD = %.3f\t| ', logLambdaD);
    lambdaD = 10^logLambdaD;
    
    % compute solution for new lambdaD
    X = delayRegu_solve(B, A, L, pairs, actTimes, lambdaL, lambdaD, X);
    residualNorm = norm(A*X-B,'fro');
    residualNormDiff = residualNorm-residualNorm0;
    
    % update history for secant method
    l(2) = l(1);
    l(1) = logLambdaD;
    r(2) = r(1);
    r(1) = residualNormDiff;
    
    relErr = abs(residualNormDiff)/residualNorm0;
    fprintf('\t| relErr = %.2e\n', relErr);
    
    % check for convergence
    if relErr < relTol
        % update history for extrapolation
        s(3) = s(2);
        s(2) = s(1);
        s(1) = l(1);
        if step > 1
            s(4) = (l(1)-l(2))/(r(1)-r(2));
        else
            s(4) = nan;
        end
        
        converged = true;
        break;
    end
    
end

if ~converged
    fprintf('SECANT METHOD DID NOT CONVERGE WITHIN %i STEPS.\n', maxSteps);
end

end