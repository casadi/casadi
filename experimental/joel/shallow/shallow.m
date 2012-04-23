%
%     This file is part of CasADi.
% 
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
% 
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
% 
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
% 
% 



sprad = 0.03;
spheight = 0.01;
g = 9.81;
numboxes = 100;
numboxes_per_plot = 1;
poolwidth = 0.2;
drag = 4.0; % => u(0)
dt = 0.001;
dt = 0.0001;
depth = 0.01; % => u(1)
ntimesteps = 10000;
numsplashes = 6;
numtomiss = 100;
dx = poolwidth/numboxes;
dy = poolwidth/numboxes;

x = linspace(0,poolwidth,numboxes);
y = linspace(0,poolwidth,numboxes);
[X,Y] = meshgrid(x,y);

% Allocate memory
u0 = zeros(numboxes+1, numboxes);
u1 = zeros(numboxes+1, numboxes);
v0 = zeros(numboxes  , numboxes+1);
v1 = zeros(numboxes  , numboxes+1);
h0 = zeros(numboxes  , numboxes);
h1 = zeros(numboxes  , numboxes);


% O = imread('~/Desktop/O.png');
% O = O(:,:,1)==0;
% O = O(end:-1:1,:);
% O = O*spheight;
% 
% % smoothen out
% for k=1:20
%     O(2:end-1,2:end-1) = (4*O(2:end-1,2:end-1) + O(1:end-2,2:end-1) + O(3:end,2:end-1) + O(2:end-1,1:end-2) + O(2:end-1,3:end))/8;
% end
% O = [O(41:end,:);O(1:40,:)];
% O = [O(:,41:end),O(:,1:40)];
% 
% P = imread('~/Desktop/P.png');
% P = P(:,:,1)==0;
% P = P(end:-1:1,:);
% P = P*spheight;

% % smoothen out
% for k=1:20
%     P(2:end-1,2:end-1) = (4*P(2:end-1,2:end-1) + P(1:end-2,2:end-1) + P(3:end,2:end-1) + P(2:end-1,1:end-2) + P(2:end-1,3:end))/8;
% end
% P = [P(41:end,:);P(1:40,:)];
% P = [P(:,41:end),P(:,1:40)];

% Initial conditions
spdist = sqrt( (X-0.04).^2 + (Y-0.04).^2);
I = spdist<sprad/3.0;
h0(I) = spheight * cos(3.0*pi*spdist(I)/(2.0*sprad));
%h0 = O;

% Clear figure
figure(1)
close
clf

% Integrate forward in time
for tt=1:ntimesteps    
%    h0(1:numboxes/10,1:numboxes/10) = cos(2*pi*tt/100)*spheight;

%    if tt==ntimesteps/20
%        h0 = h0 + P;
%    end
%    if mod(tt,ntimesteps/40)==1
%	h0(I) = h0(I) + spheight * cos(3.0*pi*spdist(I)/(2.0*sprad));
%    end


    if mod(tt,5)==1
        sh = surf(X(1:numboxes_per_plot:end,1:numboxes_per_plot:end),Y(1:numboxes_per_plot:end,1:numboxes_per_plot:end),h0(1:numboxes_per_plot:end,1:numboxes_per_plot:end));
        axis([0,poolwidth,0,poolwidth,-spheight,spheight])
        view(18,78)
        colormap('bone')
        set(gca,'CLim',[-spheight/10,spheight/10])
        drawnow
    end
    
    time = tt*dt;
    for ii=1:numboxes-1
        for jj=1:numboxes
            u1(1 + ii,jj) = u0(1 + ii,jj) + ...
                dt*(-g/dx * (h0(1+ii,jj)-h0(ii,jj))- u0(1+ii,jj)*drag);
        end
    end
    
    for ii=1:numboxes
        for jj=1:numboxes-1
            v1(ii,jj+1) = v0(ii,jj+1) + ...
                dt*(-g/dy*(h0(ii,jj+1)-h0(ii,jj))- v0(ii,jj+1)*drag);
        end
    end
    
    for ii=1:numboxes
        for jj=1:numboxes
            h1(ii,jj) = h0(ii,jj) + ...
                (-depth*dt)*(1.0/dx*(u1(1+ii,jj)-u1(ii,jj)) ...
                + 1.0/dy * (v1(ii,jj+1)-v1(ii,jj)));
        end
    end
    
    % Shift the timestep so that h0 contains h1, u0 contains u1 and v0 contains v1
    u0 = u1;
    v0 = v1;
    h0 = h1;
    
    
end



