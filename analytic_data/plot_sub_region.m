westcoast = load('../data/ukireland_coastline.NaN.dat');
l = find(isnan(westcoast(:,1)));

lonlim = [-9 -4];
latlim = [54 59]; 


Clyde = polyshape([-5.7 -4.9 -4.202 -4.938 -5.48 -5.388 -5.83], [55.024 55.024 55.9 56.6 56.044 55.816 55.3]);

% Firth of Lorn 
Lorn = polyshape([-5.48 -6.1 -4.958 -4.938],[56.04 56.3 57.2 56.6]);

%Little Minch 
LMinch = polyshape([-7.6 -7.2 -4.958 -5.530], [56.75 57.5 57.2 56.75]);

%Upper Minch
UMinch = polyshape([-7.2 -6.2 -4.20 -4.958], [57.5 58.64 58.64 57.2]);

%Outer Hebs 
Hebs = polyshape([-6.2 -7.89 -7.89 -7.6 -7.2], [58.64 58.64 56.75 56.75 57.5]);

%North Channel 
Channel = polyshape([-7.89 -5.53 -6.1 -5.48 -5.388 -5.83 -5.7], [56.75 56.75 56.3 56.04 55.82 55.3 55.02]);

%Solway Coast
SolCol = polyshape([-4.9 -5.6 -4.948 -4.365], [55.04 55.02 54.55 54.69]);

% region sub-area
sub_area3 = polyshape([-7.9 -4.20 -4.20 -7.9], [55.024 55.024 58.64 58.64]); % north and extends to north channel 

sub_areax1 = Clyde.Vertices(:,1); 
sub_areay1 = Clyde.Vertices(:,2); 

sub_areax2 = Lorn.Vertices(:,1); 
sub_areay2 = Lorn.Vertices(:,2); 

sub_areax3 = LMinch.Vertices(:,1); 
sub_areay3 = LMinch.Vertices(:,2); 

sub_areax4 = UMinch.Vertices(:,1); 
sub_areay4 = UMinch.Vertices(:,2); 

sub_areax5 = Hebs.Vertices(:,1); 
sub_areay5 = Hebs.Vertices(:,2); 

sub_areax6 = Channel.Vertices(:,1); 
sub_areay6 = Channel.Vertices(:,2); 

sub_areax7 = SolCol.Vertices(:,1); 
sub_areay7 = SolCol.Vertices(:,2); 

x = westdat(:,3);
y = westdat(:,2);
value = westdat(:,16);
size = westdat(:,6); 

x = table2array(x);
y = table2array(y);
value = table2array(value);
size = table2array(size);
figure(1)

clf;
hold on

% first draw the coast 
for i = 1:size(l)-1
    plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k', 'LineWidth', 1)
end 

plot(sub_area3, 'FaceAlpha', 0, 'LineWidth', 1)
hold on 
plot(Clyde, 'FaceColor', 'b', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on
plot(Lorn, 'FaceColor', 'c', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on 
plot(LMinch, 'FaceColor', 'g', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on 
plot(UMinch, 'FaceColor', 'r', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on
plot(Hebs, 'FaceColor', 'b', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on 
plot(Channel, 'FaceColor', 'y', 'FaceAlpha', 0.01, 'LineWidth', 1)
hold on 
% legend(plot, {'Firth of Clyde','Firth of Lorn', 'Little Minch', 'Upper Minch', 'Outer Hebrides', 'North Channel'},'Location','northwest','Orientation','vertical')
% FontSize = 12;
set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

figure(2)
plot = scatter(x, y, [], value, 'filled');
colormap(jet(7));