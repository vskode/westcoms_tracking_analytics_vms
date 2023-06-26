%% load files %%

clear
close all;

folder_path = '/home/vincent/Nicnak/third_paper/westcoms_trackingOutput/VMS/v09_3percent_vms/';
str_file = dir(folder_path);
polygons_seaRegion
polys = [Lorn, Channel, Clyde, Hebs, LMinch, UMinch];
polys_names = ["Lorn", "Channel", "Clyde", "Hebs", "LMinch", "UMinch"];

% north and extends to north channel 
sub_area = polyshape([-7.8 -4.23 -4.23 -7.8], [55.06 55.06 58.6 58.6]); 

%% load all files %%
N_particles.names = polys_names;
N_particles.counts = zeros(366, 6);
for i = 1:365
    file = load(join([folder_path, str_file(i+2).name]));
    
    for p = 1:6
        p_area = polys(p);

        [in_subarea] = inpolygon(file.x, ...
                                  file.y, ...
                                  sub_area.Vertices(:, 1), ...
                                  sub_area.Vertices(:, 2));
        
        beached_in_subarea_x = file.x(bitand(file.ashore, in_subarea));
        beached_in_subarea_y = file.y(bitand(file.ashore, in_subarea));

        boundary_x = p_area.Vertices(:,1); 
        boundary_y = p_area.Vertices(:,2); 
        
        [in_poly] = inpolygon(beached_in_subarea_x, ...
                              beached_in_subarea_y, ...
                              boundary_x, ...
                              boundary_y);

        particles_in_poly_x = beached_in_subarea_x(in_poly);
        particles_in_poly_y = beached_in_subarea_y(in_poly);
        
        N_particles.counts(i, p) = numel(particles_in_poly_x);
    end    
end

for i = 0:floor(365/7)-1
    for p = 1:6
        N_particles.weekly(i+1, p) = sum(N_particles.counts(i*7+1:i*7+8, p));
    end
end

save("Particles_per_poly.mat", "N_particles")

%% plot

figure()
title('accumulation daily')
for i = 1:6
    % subplot(6, 1, i)
    hold on
    a(i) = plot(N_particles.counts(:, i));
    % title(N_particles.names(i))
end
names = N_particles.names;
legend(a, names)

%%
figure()
title('new beached daily')
for i = 1:6
    % subplot(6, 1, i)
    hold on
    a(i) = plot(diff(N_particles.counts(305:end, i)));
    % title(N_particles.names(i))
end
xline([1, 8, 15, 22, 29, 36, 43, 50, 57], '-', color='black')
names = N_particles.names;
legend(a, names)

%%
figure()
title('new beached daily')
for i = 1:6
    subplot(6, 1, i)
    heatmap(diff(N_particles.counts(305:end, i))');
    title(N_particles.names(i))
end

%%

figure()
title('new beached daily percentage')
percent_new_beach_daily = N_particles.counts ./ sum(N_particles.counts, 2) * 100;
for i = 1:6
    subplot(6, 1, i)
    heatmap(percent_new_beach_daily(305:end, i)');
    title(N_particles.names(i))
end

%%

figure()
title('new beached daily - average')
for i = 1:6
    hold on
    plot(N_particles.counts(305:end, i)' - mean(N_particles.counts(305:end, i)));
    % title(N_particles.names(i))
end
legend()
%% weekly


figure()
title('new beached weekly')
for i = 1:6
    subplot(6, 1, i)
    heatmap(diff(N_particles.weekly(44:end, i))');
    title(N_particles.names(i))
end

%%

figure()
title('new beached weekly percentage')
percent_new_beach_weekly = N_particles.weekly ./ sum(N_particles.weekly, 2) * 100;
for i = 1:6
    hold on
    plot(percent_new_beach_weekly(44:end, i)');
    % title(N_particles.names(i))
end

%%

figure()
title('new beached weekly percentage')
for i = 1:6
    hold on
    plot(N_particles.weekly(44:end, i)' - mean(N_particles.weekly(44:end, i)));
    % title(N_particles.names(i))
end

%%

figure()
title('pies by polygons')
% for i = 1:8
    % subplot(4, 2, i)
pie(sum(N_particles.weekly(44:52, :)), N_particles.names)
% end


