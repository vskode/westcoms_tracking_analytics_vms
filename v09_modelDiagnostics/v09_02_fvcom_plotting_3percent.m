close all
clear

% v07 with update 16/08/2021
%  - code updated to only include plotting 
%  - allows plotting of output from cluster 

%% set path to WeStCOMS2 data 
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking/');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_setup_particles/');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_setup_model/');
addpath('E:/2020_particle_track_models/2021_WestComs_na/data/');
addpath('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper');
addpath('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent');
% addpath('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/v07_1percent01');

%% read user-defined parameters and setup
% set model parameters
setup_model04_alltest_3percent_resuspend_coast
% set particle postions and properties
setup_particles07_randlocations

lonlim = [-7.9 -4.23];
latlim = [55.06 58.64]; 

% coastline
westcoast = load('../data/ukireland_coastline.NaN.dat');
gap = find(isnan(westcoast(:,1))); % finds the island separators in the dataset

%% various setup
dtday = dtsec/(24*3600);  % tracking time step in days
dthour = dtsec/3600;  % tracking time step in hours
if(particlewind)
    wind = false;  % 'particlewind' not compatible with 'wind' 
end
gifname = ['v09_02_movie_' runname '.gif'];
if(isfile(gifname))
    delete(gifname)
end

%% time-stepping loop
it = 1; % sequential numbering of tracking time steps
timetrack(1) = startloop;
dayprev = 0;
time_fullrun = tic;

while(timetrack(it)<endloop)
    time_onestep = tic; % set timer for one timestep
    day = floor(timetrack(it));
    if(day~=dayprev) % check if need to read a new 25-hour model data file
       timevec = datevec(timetrack(it));
       datafile = ['v09_02_3percent_' num2str(timevec(1)) num2str(timevec(2:3),'%02.f')];
%        load(['./matfiles_25levels/' datafile]);
%        load(['E:/WeStCOMS2/mat_files_25levels2020/' datafile])
       load(['../cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent/' datafile]);
    end

    % advance time and store day of previous time
    dayprev = day;
    it = it+1;
    timetrack(it) = timetrack(it-1)+dtday;

    if(drawplot) 
    % add to the gif movie
    if(mod(it,nmovie)==2)  % add a new movie frame every nmovie steps
      
        % plot map of particles
        clf
        plot(x(active),y(active),'b.')
        hold on
        plot(x(ashore),y(ashore),'r.')
        for k = 1:size(gap)-1
          plot(westcoast(gap(k)+1:gap(k+1)-1,1),westcoast(gap(k)+1:gap(k+1)-1,2),'k')
        end 
        
        date_label = datestr(timetrack(it),'dd mmm yyyy HH:MM:ss');
        title(date_label(1:20))
        set(gca,'xlim',lonlim,'ylim',latlim)
        set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
        pause(.01)
            
        I = getframe(gcf);
        I = frame2im(I);   
          [X, map] = rgb2ind(I, 128);
          if(isfile(gifname))
             imwrite(X, map, gifname, 'GIF', 'WriteMode', 'append', 'DelayTime', 0);
          else
             imwrite(X, map,gifname, 'GIF', 'WriteMode', 'overwrite', 'DelayTime', 0, 'LoopCount', Inf);
          end
          
    else 
        continue; % dont do anything (so no plotting)
    end 
    end
    
    if(mod(it,nmovie)==0)
    fname = ['step' num2str(it)]; 
    save(fname, 'x', 'y'); 
    end
end
    toc(time_onestep)

toc(time_fullrun)
