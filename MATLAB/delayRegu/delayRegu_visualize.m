%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function delayRegu_visualize(actTimes, X, mesh, sensorNodes)

colMat = [0.8   0     0.0  ; ...
       1     0.6   0    ; ...
       0.65  0.6   0    ; ...
       0.533 0.533 0.533; ...
       0.8   0     0.8  ; ...
       0     0     1    ; ...
       0     0.6   1    ; ...
       0.45  0.45  0.8  ; ...
       0     0     0    ; ...
       0     0.733 0    ];
colMat = repmat(colMat, ceil(numel(sensorNodes)/size(colMat,1)), 1);
colMat = colMat(1:numel(sensorNodes),:);
colCell = num2cell(colMat, 2);

subplot(1,5,1:2);
h = plot(X(sensorNodes,:)', 'LineWidth',1.5);
set(h, {'color'}, colCell);
center = (prctile(X,10,'all')+prctile(X,90,'all'))/2;
ylim([center-70 center+70]);
xlim([1 size(X,2)])
xlabel('Time (ms)');
ylabel('TMV (mV)');

subplot(1,5,3:5);
visualizeDataOnMesh(mesh, actTimes, [0 100], 20, [0 0]);
c = colorbar;
c.Label.String = 'AT (ms)';
c.Label.FontSize = 11;
c.FontSize = 10;
hold on
h = scatter3(mesh.p(sensorNodes,1), mesh.p(sensorNodes,2), mesh.p(sensorNodes,3), 72, 'w', 'filled');
h.CData = colMat;
scatter3(mesh.p(sensorNodes,1), mesh.p(sensorNodes,2), mesh.p(sensorNodes,3), 72, 'w');

drawnow;

end