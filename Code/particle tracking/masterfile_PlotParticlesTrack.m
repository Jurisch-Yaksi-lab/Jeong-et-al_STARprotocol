%% Master Analysis file
% All analysis performed on epifluorescence recordings of motile cilia
% mediated flow in the zebrafish brain


%% step0 is done in imageJ/fiji--------------------------------------------------

% 1) convert image to pixel
% 2) TrackMate, with LoG detector with 10pixel blobs, a threshold of 50, 
%   - initial search radius: 15,0
%   - max search radius: 15,0
%   - max frame gap: 2
% 3) simple laptracker
% 6)export track to XML file

%% step1: import the track to Matlab and set some values ===================================================================================

% Define the path to the data
[data.filename,data.pdir] = uigetfile;
cd(data.pdir); 

%import the data
Tracks = importXMLfile(data.filename);

% input the spatial resolution, 1 pixel = 1.13 micron, this value needs to
% be measured or recovered from the microscope. In our system with a 10x
% objective we have a resolution of 1pixel=1.13micron
data.spatres=1.13;

%% step2: reorganize the tracks for each particles in a cell array ====================================================

% in the Tracks file, all data for each track are separated by some summary
% data that appears like Nan. Finding Nan will tell us where each track
% starts and where each track stops
temp=find(isnan(double(Tracks{:,4}))); % find the Nan that surrounds each track

% identify the track onset and offset
trk.onset=temp(1:2:end)+1; % this is the row in the Track matrix where each track starts
trk.offset=temp(2:2:end)-1; % this is the row in the Track matrix where each track stops

% make a cell array "track" where each cell correspond to the positions of
% one particle over time (first column is x, second column is y, each row
% is one time interval)
 for i=1: size(trk.offset,1)
 track{i,1}=cat(2,Tracks{trk.onset(i):trk.offset(i),4},Tracks{trk.onset(i):trk.offset(i),6}); 
 end
 
 % organize the data in a structure array
 data.track=track;
 
 % clear temp and Tracks from the workspace
 clear temp Tracks trk track
%% step3: plot some tracks to make sure that everything is fine, you can plot few tracks ===================================
% this step is to check that everything is fine, can be omitted 
figure, 
 for i=1:1000
     plot(data.track{i}(:,1),-data.track{i}(:,2)) %here -y is plotted to keep the right orientation
     hold on
 end
 title(' first 1000 tracks plotted')

%% step 4: identify the good tracks based on some threshold values inputed by the experimenter ==================================================================================================
% for instance identify track longer than "tracklength" = e.g. 40
% datapoints and than "trackdistance" = e.g., 20 pixel between the start
% and end of the track

% input the threshold values
data.tracklength=40;
data.trackdistance=20;

% data.trackid will be a matrix with the id of tracks passing the decided
% criteria

count=1;
for i=1: size(data.track,1)
  if size(data.track{i},1)>data.tracklength & pdist(data.track{i}([1,end],:))>data.trackdistance;
    data.trackid(count)=i;
    count=count+1;    
  end
end

% clear i and count from workspace
clear i count
 %% step 5: plot the selected tracks, each track has a different color ==============================================
 figure, 
 for i=data.trackid
    plot(data.track{i}(:,1),-data.track{i}(:,2), 'LineWidth', 2) %here -y is plotted to keep the right orientation
    hold on
 end
title(['total of ' num2str(size(data.trackid,2)) ' tracks with min ' num2str(data.tracklength) ' datapoints and ' num2str(data.trackdistance) ' pixels'])
axis off; 
axis image
line([0 100/data.spatres],[-100 -100],'Color', 'k', 'LineWidth' , 2 );
text (1, -130, '100µm')

%% step 6: plot selected tracks with color based on the angle of the point in relation with 1 dot later====================
% to use this approach, we use the plot function and draw each segment of
% the track with teh angle color. Tuse this approach, we calculate the
% angle between time point t and time point t+1, which is the value "angle"
% and use a color in the hsv colormap. Since angle is in radian (-pi to
% +pi), we first need to convert it into a value between 1 and 256. To do
% so, we take the (angle in radius +pi)/(2*pi)*256

% be aware that plotting takes time

figure;
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8])
ax1=subplot('Position', [0.05 0.1 0.7 0.8]);

color=hsv;
tic % allows to report the time needed to plot

 for i=data.trackid
    for k=1: size(data.track{i},1)-1
    angle= ceil(((atan2(data.track{i}(k+1,2)-data.track{i}(k,2),data.track{i}(k+1,1)-data.track{i}(k,1)))+pi)/(2*pi)*256);
    plot(data.track{i}(k:k+1,1),-data.track{i}(k:k+1,2), 'Color', color(angle,:), 'LineWidth', 3)
    hold on
 end
 end
title(['direction of a total of ' num2str(size(data.trackid,2)) ' tracks with min ' num2str(data.tracklength) ' datapoints and ' num2str(data.trackdistance) ' pixels'])
axis off; 
axis image
line([0 100/data.spatres],[-100 -100],'Color', 'k', 'LineWidth' , 2 );
text (10, -120, '100µm')
toc

% this will plot the direction of the tracks in a colorwheel
ax2=subplot('Position', [0.85 0.7 0.1 0.1]);
create_color_wheel()
set(gca,'xdir','reverse') % the x axis needs to be rotated to reflect the right direction

% Export the figure
saveas(gcf,[data.pdir, data.filename,'tracks.png'])

%% step 7: save the data======================================================

% all the data is arranged in one structure array data
save(fullfile([data.pdir, data.filename '_tracks.mat']),'data','-v7.3');