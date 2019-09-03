%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

%% Version 1:
%  Differences are calculated for all (overlapping and non-overlapping) 
%  timesteps of time-shifted signals. For non-overlapping timesteps, the
%  difference is calculated with the closest timestep of the other signal,
%  i.e. its very first or very last timestep.

function D = delayRegu_assembleDifferenceMatrix(pairs, actTimes, T)

% T: number of timesteps
P = size(pairs,1);      % number of pairs
N = length(actTimes);   % number of nodes

delays = round(actTimes(pairs(:,1))-actTimes(pairs(:,2)));
pairs(delays < 0,:) = fliplr(pairs(delays < 0,:));
delays = abs(delays(:));
cumsumDelays = [0; cumsum(delays)];

i = NaN(1, 2*(P*T+cumsumDelays(end)));
j = i;
v = i;
k = 0;

for p = 1:size(pairs,1)
    
    d = delays(p);
    r = (p-1)*T+cumsumDelays(p)+1;  % start index of rows
    c1 = (pairs(p,1)-1)*T+d+1;      % start index of minuend columns
    c2 = (pairs(p,2)-1)*T+1;        % start index of subtrahend columns
    
    i(k+1:k+T+d) = r:r+T+d-1;
    j(k+1:k+T) = c1-d:c1+T-d-1;
    j(k+T+1:k+T+d) = repmat(c1+T-d-1, 1, d);
    v(k+1:k+T+d) = ones(1,T+d);
    k = k+T+d;
    
    i(k+1:k+T+d) = r:r+T+d-1;
    j(k+1:k+d) = repmat(c2, 1, d);
    j(k+d+1:k+T+d) = c2:c2+T-1;
    v(k+1:k+T+d) = -ones(1,T+d);
    k = k+T+d;
    
end

% Permute columns of D (timestep-major -> node-major):
% (time1 time2 ... timeT)*N -> (node1 node2 ... nodeN)*T
permut = reshape(repmat((1:N:T*N)',1,N) + repmat(0:N-1,T,1), T*N, 1);
j = permut(j);

D = sparse(i, j, v, P*T+cumsumDelays(end), T*N);

end

%% Version 2:
%  Differences are calculated only for overlapping timesteps.
%  No differences are calculated for non-overlapping timesteps.
%  This means some timesteps at the beginning and the end of the signal may
%  not experience any delay-based regularization.

% function D = delayRegu_assembleDifferenceMatrix(pairs, actTimes, T)
% 
% % T: number of timesteps
% P = size(pairs,1);      % number of pairs
% N = length(actTimes);   % number of nodes
% 
% delays = round(actTimes(pairs(:,1))-actTimes(pairs(:,2)));
% pairs(delays < 0,:) = fliplr(pairs(delays < 0,:));
% delays = abs(delays);
% 
% i = NaN(1, 2*(P*T-sum(delays)));
% j = i;
% v = i;
% k = 0;
% 
% for p = 1:P
% 
%     d = delays(p);
%     r = (p-1)*T+1;              % start index of rows
%     c1 = (pairs(p,1)-1)*T+d+1;  % start index of minuend columns
%     c2 = (pairs(p,2)-1)*T+1;    % start index of subtrahend columns
%     l = T-d;                    % number of overlapping timesteps
%     
%     i(k+1:k+l) = r:r+l-1;
%     j(k+1:k+l) = c1:c1+l-1;
%     v(k+1:k+l) = ones(1,l);
%     k = k+l;
%     
%     i(k+1:k+l) = r:r+l-1;
%     j(k+1:k+l) = c2:c2+l-1;
%     v(k+1:k+l) = -ones(1,l);
%     k = k+l;
% 
% end
% 
% % Permute columns of D (timestep-major -> node-major):
% % (time1 time2 ... timeT)*N -> (node1 node2 ... nodeN)*T
% permut = reshape(repmat((1:N:T*N)',1,N) + repmat(0:N-1,T,1), T*N, 1);
% j = permut(j);
% 
% D = sparse(i, j, v, P*T, T*N);
% 
% end
