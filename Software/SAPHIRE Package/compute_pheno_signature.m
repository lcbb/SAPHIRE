%%%%%%%%%%%%%%%%%%%%
% Computes dynamic features from the temporal state sequence annotations.
%
% Inputs:   SAPHIRE_Model - Cell array of SAPHIRE models for individual
%                           cell trajectories
%
% Outputs:  phenoComboDynamicFeats - structure array of composite dynamic
%                                    features derived from all cell 
%                                    trajectories in SAPHIRE_Model
%           phenoComboDynamicSig - composite phenotypic signature of
%                                  dynamics derived from all cell 
%                                  trajectories in SAPHIRE_Model
%               
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015 Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function [phenoComboDynamicFeats,phenoComboDynamicSig] = compute_pheno_signature(SAPHIRE_Model)

%Compute dynamics features from the inferred SAPHIRE annotations of most
%likely hidden states and mean values to generate signatures of time spent
%in different slices of PC space as well as directional transitions in PC
%space.

numTrajs = numel(SAPHIRE_Model);
SAPHIREIndivTrackDynamicFeats.numStates = cell2mat(cellfun(@(x) max(x.ML_states),SAPHIRE_Model,'uniformoutput',0));

numSlices = 12; %set number of slices to partition polar PC space into (360 degrees are divided into this number of slices)
sliceEdges = 0:360/numSlices:360; %edges of bins used

for i = 1:numTrajs
    
    stateSeq = SAPHIRE_Model{i}.ML_states;
    
    %STATE LOCATIONS IN PC SPACE:
    
    %Save frequencies of states with their respective means, converted to
    %polar coordinates (r, theta). For each trajectory we get a Kx3 matrix, where K is
    %the unique number of states (row) and for each row's columns the
    %following are stored : [theta, r, sigma].
    
    SAPHIREIndivTrackDynamicFeats.PC_state_params{i} = ...
        zeros(SAPHIREIndivTrackDynamicFeats.numStates(i),3); 
    
    for j = 1:numel(SAPHIRE_Model{i}.ML_params.sigma_emit)
       
        %Convert the current state's location in PC space to polar
        %coordinates
        [thetaTemp,rTemp] = cart2pol(SAPHIRE_Model{i}.ML_params.mu_emit(1,j),...
            SAPHIRE_Model{i}.ML_params.mu_emit(2,j));
        
        %Convert the displacement of angles in radians to degrees
        thetaTemp = rad2deg(thetaTemp);
        
        %Convert angles that are negative to be positive in counterclockwise
        %direction (i.e. -90 becomes 270)
        if thetaTemp < 0
            thetaTemp = 360 + thetaTemp;
        end
        
        SAPHIREIndivTrackDynamicFeats.PC_state_params{i}(j,:) = ...
            [thetaTemp,rTemp,SAPHIRE_Model{i}.ML_params.sigma_emit(j)];
    end
    
    %Calcuate the PC slice dwell frequencies (this is fraction of time a
    %cell's inferred states spend in a given slice for |numSlices| total number of slices).
    tempThetas = SAPHIREIndivTrackDynamicFeats.PC_state_params{i}(:,1); %thetas of states
    SAPHIREIndivTrackDynamicFeats.stateLocationFreqSlices{i} = ...
        histcounts(tempThetas(stateSeq),sliceEdges)./length(stateSeq);
    
    %Assign the state radii (from PC origin) to whichever slices they are
    %in. That is, for each slice we place the radii values of the states
    %that fall into that slice from the cell's trajectory.
    [~,idx] = histc(tempThetas,sliceEdges);
    SAPHIREIndivTrackDynamicFeats.stateLocationRadiiSlices{i} = cell(1,numSlices);
    uniqueLoc = unique(idx);
    for j = 1:numel(uniqueLoc)
        SAPHIREIndivTrackDynamicFeats.stateLocationRadiiSlices{i}{1,uniqueLoc(j)} = ...
            SAPHIREIndivTrackDynamicFeats.PC_state_params{i}(idx==uniqueLoc(j),2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %STATE TRANSITIONS IN PC SPACE (REGARDLESS OF LOCATION):
    
    %Compute the dwell frequencies in states that are TRANSITIONS from
    %other states, where the theta corresponds to the transition direction
    %and dwell frequency to the transitioned-to state. Also compute the
    %magnitude of the transitions between states in terms of distances in
    %PC space.
    locTrans = find(abs(stateSeq(2:end)-stateSeq(1:end-1))); %find locations where state transitions occur (point just before transition)
    temp = zeros(1,length(stateSeq));
    temp(locTrans+1) = 1;
    idx = cumsum(temp)+1;
    %Partition the states over time into cell array
    statePartitions = accumarray(idx',stateSeq',[],@(x) {x});
    transDirArray = cell(numel(statePartitions)-1,1); %store directions of transtions (first state isn't transitioned to so ignored)
    transMagArray = cell(numel(statePartitions)-1,1); %store magnitudes of transitions (first state isn't transitioned to so ignored)
    for j = 2:numel(statePartitions)
        transDirCart = SAPHIRE_Model{i}.ML_params.mu_emit(:,unique(statePartitions{j})) - ...
        SAPHIRE_Model{i}.ML_params.mu_emit(:,unique(statePartitions{j-1})); %cartesian coordinates of transition between states
        transDirTheta = rad2deg(cart2pol(transDirCart(1),transDirCart(2))); %direction of transition from (j-1)th to (j)th state.
        if transDirTheta < 0
            transDirTheta = 360 + transDirTheta;
        end
        transMag = norm(transDirCart); %magnitude of the transition
        transDirArray{j-1} = repmat(transDirTheta,numel(statePartitions{j}),1);
        transMagArray{j-1} = repmat(transMag,numel(statePartitions{j}),1);
    end
    
    %Some extra features to compute:
    %Store the frequency of transitions between ANY state for a given
    %trajectory (i.e. number of times a transition happens, normalized to
    %trajectory length).
    SAPHIREIndivTrackDynamicFeats.freqTransAnyStates(i,1) = (numel(statePartitions)-1)/length(stateSeq);
    %Store the variance of frequencies of state partitions (i.e.
    %periods between state transitions)
    SAPHIREIndivTrackDynamicFeats.varFreqStatePartitions(i,1) = var(cell2mat(cellfun(@numel,statePartitions,'uniformoutput',0))./length(stateSeq));
    
    %Calcuate the state transition direction frequencies (this is fraction
    %of time points transitions in direction of  given slice occur for
    %|numSlices| total number of slices). Note that longer residence states
    %get more frequency weighing from previous transition since we are more
    %confident this is a robust (longer) transition by the cell, as opposed
    %to transient, rapid, single transition lasting a single time point,
    %for example. Also store the displacements (magnitudes) of the
    %transitions, with same number of repeats based on time spend in
    %transitioned to state as for the frequencies.
    if SAPHIREIndivTrackDynamicFeats.numStates(i) == 1 %if 1 state only then no transitions
        SAPHIREIndivTrackDynamicFeats.stateTransitionFreqSlices{i} = zeros(1,numSlices);
        SAPHIREIndivTrackDynamicFeats.stateTransitionMagSlices{i} = repmat({0},1,numSlices);
    else
        concatTransDirArray = vertcat(transDirArray{:});
        concatTransMagArray = vertcat(transMagArray{:});
        SAPHIREIndivTrackDynamicFeats.stateTransitionFreqSlices{i} = histcounts(concatTransDirArray,sliceEdges)./length(stateSeq);
        %Assign the state radii (from PC origin) to whichever slices they are
        %in. That is, for each slice we place the radii values of the states
        %that fall into that slice from the cell's trajectory.
        [~,idx] = histc(concatTransDirArray,sliceEdges);
        SAPHIREIndivTrackDynamicFeats.stateTransitionMagSlices{i} = cell(1,numSlices);
        uniqueLoc = unique(idx);
        for j = 1:numel(uniqueLoc)
            SAPHIREIndivTrackDynamicFeats.stateTransitionMagSlices{i}{1,uniqueLoc(j)} = concatTransMagArray(idx==uniqueLoc(j));
        end
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Combine the histogram features for all the cells computed above into a
%dynamic phenotypic signature:

%State mean (center) locations:
tempAllRadii = vertcat(SAPHIREIndivTrackDynamicFeats.stateLocationRadiiSlices{:});
tempAllRadiiMeans = cell2mat(cellfun(@(x) mean(x),tempAllRadii,'uniformoutput',0));
tempAllRadiiMeans(isnan(tempAllRadiiMeans)) = 0;
tempAllFreq = vertcat(SAPHIREIndivTrackDynamicFeats.stateLocationFreqSlices{:});

if numTrajs  > 1
    %Mean across all cells of normalized state location dwell times
    phenoComboDynamicFeats.meanStateLocFreq = mean(tempAllFreq)';
    %Mean across all cells of state location radii
    phenoComboDynamicFeats.meanStateLocRadii = mean(tempAllRadiiMeans)';
else
    %Mean across all cells of normalized state location dwell times
    phenoComboDynamicFeats.meanStateLocFreq = tempAllFreq';
    %Mean across all cells of state location radii
    phenoComboDynamicFeats.meanStateLocRadii = tempAllRadiiMeans';
end
    
%State transitions:
tempAllTransDirMags = vertcat(SAPHIREIndivTrackDynamicFeats.stateTransitionMagSlices{:});
tempAllTransDirMagsMeans = cell2mat(cellfun(@(x) mean(x),tempAllTransDirMags,'uniformoutput',0));
tempAllTransDirMagsMeans(isnan(tempAllTransDirMagsMeans)) = 0;
tempAllTransDirFreqs = vertcat(SAPHIREIndivTrackDynamicFeats.stateTransitionFreqSlices{:});

if numTrajs  > 1
    %Mean across all cells of normalized state transition dwell times
    phenoComboDynamicFeats.meanTransDirFreq = mean(tempAllTransDirFreqs)';
    %Mean across all cells of state transition magnitudes
    phenoComboDynamicFeats.meanTransDirMag = mean(tempAllTransDirMagsMeans)';
else
    %Mean across all cells of normalized state transition dwell times
    phenoComboDynamicFeats.meanTransDirFreq = tempAllTransDirFreqs';
    %Mean across all cells of state transition magnitudes
    phenoComboDynamicFeats.meanTransDirMag = tempAllTransDirMagsMeans';
end

%Assemble the dynamic features into a composite signature
phenoComboDynamicSig = [phenoComboDynamicFeats.meanStateLocFreq;...
                        phenoComboDynamicFeats.meanStateLocRadii;...
                        phenoComboDynamicFeats.meanTransDirFreq;...
                        phenoComboDynamicFeats.meanTransDirMag];
