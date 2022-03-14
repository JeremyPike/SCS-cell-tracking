

% Matlab script to peform postporoessing of cell tracks and SCS segmentations
% Lok at al. (2022) 
% Author: Jeremy Pike, Image analyst for COMPARE, j.a.pike@bham.ac.uk


clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direcotry containing data
dataDir = "C:\Users\pikej\Documents\kai\data";
% minimum track length
minTrackLength = 5;
% number of channels in data
numChannels = 4;
% number of channels from ilastik
numIlastikChannels = 4;
% maximum distance from SCS to look at spots (microns)
scsMaxDistance = 100; 
% number of z slices to include (to all movies have same dims)
z_crop = 13;
% number of frames to include (to all movies have same dims)
t_crop = 32;
% frame binning
t_bin_width = 3;
% distance tolerance around SCS boundary for track to be defined as outside
buff_dist_out = 5;
% load csv file containing list of files to be processed
meta = readtable(strcat(dataDir, filesep, "meta_filt.csv"));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove files not suitable for analysis
meta = meta(meta.Use == 1, :);

% create variables to hold output measurments
meanExternalMapDisplacement =  zeros(size(meta,1), 1);
meanLeavingMapDisplacement = zeros(size(meta,1), 1);
meanEnteringMapDisplacement = zeros(size(meta,1), 1);
meanFrameMapDisplacement =  zeros(size(meta,1), 1);
numTracksEntering = zeros(size(meta,1), 1);
numTracksLeaving = zeros(size(meta,1), 1);
numTracksExternal = zeros(size(meta,1), 1);
numTracksInternal = zeros(size(meta,1), 1);
numTracks = zeros(size(meta,1), 1);
allExternalMapDisplacementWT = [];
allFrameMapDisplacementWT =[];
allExternalMapDisplacementKO = [];
allFrameMapDisplacementKO = [];

% loop through files
for i = 1: size(meta,1)
    
    row = meta(i, :);
    % retrieve path for ilastik pixel classification
    [~, filename, ext] = fileparts(row.ilastikPath);
    filepath = strcat(dataDir, filesep, filename, ext);
    % retrieve path for TrackMate results
    [~, filename, ext] = fileparts(row.trackingPath);
    trackingPath = strcat(dataDir, filesep, filename, ext);    
  
    % load ilastik classification probability map
    tif_info = imfinfo(filepath);
    ilastikProb = zeros(tif_info(1).Width, tif_info(1).Height, z_crop, t_crop+t_bin_width-1, numIlastikChannels);
    for t = 1 : t_crop+t_bin_width-1
        for z = 1 : z_crop
            for c = 1 : numIlastikChannels       
                 tifIndex =  c + (z - 1) * numIlastikChannels + (t - 1) * z_crop * numIlastikChannels;
                 ilastikProb(:, :, z, t, c) = imread(filepath, tifIndex);
            end
        end
    end
    
    % variables for processed SCS segmentation and distance maps
    scs = zeros(tif_info(1).Width, tif_info(1).Height, z_crop, t_crop);
    D = zeros(tif_info(1).Width, tif_info(1).Height, z_crop, t_crop);
    
    % for frame binning
    t_bin_pad = (t_bin_width - 1) / 2;
    
    % loop through timepoints
    for t = 1 : t_crop
        
        % mean frame binning of ilastik probability map
        ilastikProb_t = mean(ilastikProb(:, :, :, t:t+2*t_bin_pad, :), 4);
        % convert probailities to segmentation by assigning pixel class
        % with maximum probability
        [~, segmentation] = max(ilastikProb_t, [], 5);

        % to hold processed segmentation
        segmentation_processed_t = zeros(size(segmentation));

        % create coordinate matrices
        [X, Y] = meshgrid(1:size(segmentation, 1), 1:size(segmentation, 2));

        % loop through all classifications apart from background
         for c = 1 : numIlastikChannels - 1

            % class segmentation (plus one as first class is background) 
            comp = segmentation == c + 1;
            % fill any holes
            comp = imfill(comp,'holes');

            % if class segmentation is not empty
            if sum(comp(:))~=0
                % find connected components
                labels = bwlabeln(comp);
                % calculate volume of each cc
                stats = regionprops3(labels, 'Volume');
                [~, maxInd] = max(stats.Volume);
                % only keep largest cc (clear everything else)    
                comp = labels == maxInd;
            end

            % assign processed class channel
            segmentation_processed_t(comp) = c;

         end

        % SCS defined by combining all processed class segmentations
        scs_t = segmentation_processed_t >= 1;
        % fill any 2D holes in SCS 
        for z = 1 : z_crop
            scs_t(:,:,z) = imfill(scs_t(:,:,z),'holes');
        end

        % retrieve scaled voxel size (in microns)  
        voxelSizeScaled = [row.voxelX_scaled row.voxelY_scaled row.voxelZ_scaled];
        
        % create distance map representing distance from SCS in microns
        D_t = bwdistsc(scs_t, voxelSizeScaled);
        
        scs(:, :, :, t) = scs_t;
        D(:, :, :, t) = D_t;
    end

    spots_raw = readtable(trackingPath);
    
    % frame cropping
    spots_raw.FRAME = spots_raw.FRAME + t_bin_pad;
    spots_raw = spots_raw(spots_raw.FRAME >= 0, :);
    spots_raw = spots_raw(spots_raw.FRAME < t_crop, :);
    
    % convert spot spots to pixel positions in scaled images
    spots_raw.PIXEL_X_scaled = floor(spots_raw.POSITION_X / row.voxelX_scaled) + 1;
    spots_raw.PIXEL_Y_scaled = floor(spots_raw.POSITION_Y / row.voxelY_scaled) + 1;
    spots_raw.PIXEL_Z_scaled = floor(spots_raw.POSITION_Z / row.voxelZ_scaled) + 1;
    
    % remove spots greater than specifed distance from SCS
    check = [];
    for s = 1 : size(spots_raw, 1)
        if (D(spots_raw.PIXEL_Y_scaled(s), spots_raw.PIXEL_X_scaled(s), spots_raw.PIXEL_Z_scaled(s), spots_raw.FRAME(s)) <= scsMaxDistance)  
            check(s) = true;
        else
            check(s) = false;
        end
    end
    spots = spots_raw(logical(check), :);
    
    % retrieve track IDs
    trackIDUnique = unique(spots.TRACK_ID);
    trackIDUnique = trackIDUnique(trackIDUnique>0);
    
    track_type = zeros(size(trackIDUnique, 1), 1);
    % track count
    count = 0;
    % to hold number of spots in each track
    numSpots = [];
    % to hold normalised total track displacements
    externalMapDisplacement = [];
    leavingMapDisplacement = [];
    enteringMapDisplacement = [];
    % to hold framewise displacment towards SCS
    frameDisplacement = [];
    spots_processed = [];
    % loop through all tracks    
    for t = 1 : size(trackIDUnique, 1)
        % get spots in current track
        trackSpots = spots(spots.TRACK_ID == trackIDUnique(t),:);
        % track meets minumum length requirements
        if size(trackSpots, 1) >= minTrackLength
            % increase track count by one
            count = count + 1;
            % record number of spots in track
            numSpots(count) = size(trackSpots, 1);
            spots_processed = [spots_processed; trackSpots];
            % for each spot in track record distance to SCS
            trackD = [];
            for s = 1 : numSpots(count)
                trackD(s) = D(trackSpots.PIXEL_Y_scaled(s), trackSpots.PIXEL_X_scaled(s), trackSpots.PIXEL_Z_scaled(s), spots_raw.FRAME(s));
            end   
            % distance to SCS at start of track
            DStart = trackD(1);
            % distance to SCS at end of track
            DEnd = trackD(numSpots(count));

            % if track starts outside SCS
            if DStart > buff_dist_out
               
                % if track finished in SCS add one to entering count
                if DEnd == 0
                    numTracksEntering(i) = numTracksEntering(i) + 1;
                    track_type(t) = 1;
                    enteringFrame = find(trackD == 0, 1);
                    enteringMapDisplacement = [enteringMapDisplacement, DStart / (trackSpots.FRAME(enteringFrame) - trackSpots.FRAME(1))];
                else 
                % if track finished outside SCS add one to extrernal count    
                    numTracksExternal(i) = numTracksExternal(i) + 1;
                    track_type(t) = 4;
                    % framewise displacement towards SCS (postive is movement
                    % towards, negative movement away
                    trackFrameDisplacement = mean(- (trackD(2:numSpots(count)) - trackD(1:numSpots(count)-1)));
                    % record framewise displacment for this track
                    frameDisplacement = [frameDisplacement, trackFrameDisplacement];
                     % record normalised track displacment towards SCS
                    externalMapDisplacement = [externalMapDisplacement, -(DEnd - DStart) / (trackSpots.FRAME(numSpots(count)) - trackSpots.FRAME(1))];
                end
            end
            % if track starts in SCS
            if DStart == 0 
                % if track finishes in SCS
                if DEnd > buff_dist_out
                    % otherwise add one to leaving SCS count
                    numTracksLeaving(i) = numTracksLeaving(i) + 1;
                    track_type(t) = 2;
                    leavingMapDisplacement = [leavingMapDisplacement, - DEnd / (trackSpots.FRAME(numSpots(count)) - trackSpots.FRAME(find(trackD, 1) - 1))];
                elseif DEnd==0
                   % add one to internal SCS count
                    numTracksInternal(i) = numTracksInternal(i) + 1;
                    track_type(t) = 3;
                end
            end
        end
    end
    % total number of tracks for file
    numTracks(i) = count;
    % mean external map displacement for file
    meanExternalMapDisplacement(i) = mean(externalMapDisplacement);
    % mean external map displacement for file
    meanLeavingMapDisplacement(i) = mean(leavingMapDisplacement);
    meanEnteringMapDisplacement(i) = mean(enteringMapDisplacement);
    % mean external frame displacement for file
    meanFrameMapDisplacement(i) = mean(frameDisplacement);
   
end

% store results in table
meta.meanExternalMapDisplacement = meanExternalMapDisplacement;
meta.meanLeavingMapDisplacement = meanLeavingMapDisplacement;
meta.meanEnteringMapDisplacement = meanEnteringMapDisplacement;
meta.meanFrameMapDisplacement = meanFrameMapDisplacement;
meta.numTracks = numTracks;
meta.numTracksInternal = numTracksInternal;
meta.numTracksLeaving = numTracksLeaving;
meta.numTracksEntering = numTracksEntering;
meta.internalTrackRat = numTracksInternal ./ (numTracksInternal + numTracksLeaving + numTracksEntering);
meta.leavingTrackRat = numTracksLeaving ./ (numTracksInternal + numTracksLeaving + numTracksEntering);
meta.enteringTrackRat = numTracksEntering ./ (numTracksInternal + numTracksLeaving + numTracksEntering);

% save as csv file
writetable(meta, "scs_analysis.csv")
