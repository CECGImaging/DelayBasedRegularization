%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function signal = addwhitenoise(signal, SNR, refTimeInterval)

if nargin < 3
    refTimeInterval = [1 size(signal,2)];
end

sigPower = mean( mean(abs(signal(:,refTimeInterval(1):refTimeInterval(2))).^2,2) );   
sigPowerDecibel = 10*log10(sigPower);

sigma = sqrt(10^((sigPowerDecibel-SNR)/10));
noise = sigma * randn(size(signal));
signal = signal + noise;

end