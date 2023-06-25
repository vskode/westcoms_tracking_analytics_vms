%% ---** create polygons **---

% Clyde Sea  
% -- counts everything inside the Clyde Sea Area 
Clyde = polyshape([-5.971 -4.23 -4.23 -4.938 -5.48 -5.388], [55.06 55.06 56.4 56.4 56.06 55.816]);

% Firth of Lorn 
Lorn = polyshape([-5.48 -6.1 -4.958 -4.938],[56.04 56.3 57.2 56.6]);

%Little Minch 
LMinch = polyshape([-7.6 -7.2 -4.958 -5.530], [56.75 57.5 57.2 56.75]);

%Upper Minch
UMinch = polyshape([-7.2 -6.2 -4.5 -4.958], [57.5 59 59.1 57.2]);

%Outer Hebs 
Hebs = polyshape([-6.2 -9 -9 -7.6 -7.2], [58.6 58.5 56.75 56.75 57.5]);

%North Channel 
Channel = polyshape([-8.5 -5.53 -6.1 -5.48 -5.388 -5.83 -5.553], [56.75 56.75 56.3 56.04 55.82 55.3 55.02]);

%Solway Coast
SolCol = polyshape([-4.9 -5.553 -4.948 -4.365], [55.04 55.02 54.55 54.69]);
