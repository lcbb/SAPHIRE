%%%%%%%%%%%%%%%%%%%%
% Demo of the SAPHIRE package for cellular processing and probabilistic time series modeling of
% live-cell, image-based phenotypic (e.g. shape) measurements.
%               
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% USER SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run full SAPHIRE 'analysis' demo or 'plot' demo with pre-loaded demo data
%or 'QC' GUI on a sample cell trajectory?
demoType = 'plot';

%Select which example cell trajectory to run demo for (1 - 4)
trajToAnalyze = 4; 

%Set pixel scale in imaging experiment (microns per pixel):
pxScale = 1.3;

%Set maximum number of possible hidden states to test for a cell trajectory
Kmax = 5; 

%Set parallel processing: 'on' for parallel fitting of models for the
%trajectories, 'off' otherwise. NOTE: set to 'off' for parallel processing
%of cell trajectories instead of models.
mcmc_params.parallel = 'off'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load live-cell images and masks
load(fullfile(pwd,'data','DemoData.mat'));

trajToAnalyze = min(trajToAnalyze,numel(CellData.actinMask));

%Read in from file the features to compute for each cell, set by user
[~,featLabels] = xlsread('Features.xlsx');

%Add scripts in Software folder to path
addpath(genpath(pwd));

if strcmp(demoType,'QC') %Cell image time series editing tool demo
  
    %Create temporary structure with all data needed for GUI
    numFrames = numel(CellData.Iactin{trajToAnalyze});
    Iarray.CB = CellData.Iactin{trajToAnalyze};
    Iarray.Nuc = CellData.Inuc{trajToAnalyze};
    Iarray.CBMask = cellfun(@(x) full(x),CellData.actinMask{trajToAnalyze},'uniformoutput',0);
    Iarray.NucCentMask = cellfun(@(x) full(x),CellData.nucCentMask{trajToAnalyze},'uniformoutput',0);
    Iarray.phenoLabels = repmat({'s'},numFrames,1);

    %Run GUI for editing cell trajectory (NOTE: all GUI data stored in
    %|dataGUI| variable in current workspace!)
    SAPHIRE_cellTrajQC_main(Iarray);
    
    %Update frame cell body masks
    Iarray.CBMaskEdited = cell(numFrames,1);
    for i = 1:numFrames
        [rLoc,cLoc] = find(dataGUI.rgbPaddedCellImgsTileIdx == i);
        %Get temp edited cell body mask for current frame padded
        IcbMaskEditTemp = dataGUI.editedCBMaskTile(unique(rLoc),unique(cLoc));
        %Get coordinates of cropped image border for current frame
        frameFillRegion = dataGUI.pad_frameFillTile(unique(rLoc),unique(cLoc));
        [rLoc,cLoc] = find(frameFillRegion);
        Iarray.CBMaskEdited{i} = IcbMaskEditTemp(min(rLoc):max(rLoc),min(cLoc):max(cLoc));
    end

    %Update edited phenotype labels. If cell labeled as apoptotic, make all
    %subsequent frames for given cell labeled as apoptotic as well.
    locApop = min(find(ismember(dataGUI.editedPhenoLabels,'a')));
    Iarray.phenoLabelsEdited = dataGUI.editedPhenoLabels;
    Iarray.phenoLabelsEdited(locApop:end) = cellstr('a');

    %Update which frames will be deleted 
    Iarray.deleteFrameBool = dataGUI.deleteFrameMask;
    
    %Save the edited cellular time series frames
    save(fullfile(pwd,'data','IarrayQCedited.mat'),'Iarray');
    
else
    
    if strcmp(demoType,'analysis')

        %Compute features from cell shape trajectories
        CellData.shapeFeatures = compute_features(CellData,pxScale,featLabels);
        %Save Demo data with features
        save(fullfile(pwd,'data','DemoData.mat'),'CellData');

    end

    %Perform PCA on normalized cell trajectory shape features
    feats = struct2cell(structfun(@(x) vertcat(x{:}),CellData.shapeFeatures,'uniformoutput',0));
    [~,score] = pca(zscore(horzcat(feats{:})));

    %Compute shape space cellular trajectories
    locIndivTraj = [];
    for trajNum = 1:numel(CellData.actinMask)
       locIndivTraj = [locIndivTraj;repmat(trajNum,numel(CellData.actinMask{trajNum}),1)];
    end
    for trajNum = 1:max(locIndivTraj)
        CellData.shapeSpaceTraj{trajNum,1} = score(locIndivTraj==trajNum,1:2)';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(demoType,'analysis')

        disp(sprintf('\n%s\n','Running SAPHIRE Analysis:'))

        %Run SAPHIRE analysis on the shape space trajectories - infer states,
        %parameters, and most likely temporal state sequence for each
        %trajectory
        SAPHIRE_Model = cell(numel(CellData.shapeSpaceTraj),1);
        if strcmp(mcmc_params.parallel,'off') 
            parpool();
            disp(sprintf('\n%s\n','Progress:'))
            parfor_progress(numel(CellData.shapeSpaceTraj));
            tic
            parfor trajNum = 1:numel(CellData.shapeSpaceTraj)
                [SAPHIRE_Model{trajNum}.PrM,...
                SAPHIRE_Model{trajNum}.ML_states,...
                SAPHIRE_Model{trajNum}.ML_params,...
                SAPHIRE_Model{trajNum}.full_results,...
                ~,...
                SAPHIRE_Model{trajNum}.logI]...
                    = hmm_process_cell_trajectory(CellData.shapeSpaceTraj{trajNum},Kmax,mcmc_params);
                parfor_progress;
            end
            toc
            delete(gcp)
            parfor_progress(0);
        else
           for trajNum = 1:numel(CellData.shapeSpaceTraj)
                tic
                [SAPHIRE_Model{trajNum}.PrM,...
                SAPHIRE_Model{trajNum}.ML_states,...
                SAPHIRE_Model{trajNum}.ML_params,...
                SAPHIRE_Model{trajNum}.full_results,...
                ~,...
                SAPHIRE_Model{trajNum}.logI]...
                    = hmm_process_cell_trajectory(CellData.shapeSpaceTraj{trajNum},Kmax,mcmc_params);
                toc
                disp(['Modeling cell trajectory ',num2str(trajNum),' of ',num2str(numel(CellData.shapeSpaceTraj)),' complete.']);
           end 
        end

        %Save SAPHIRE models
        save(fullfile(pwd,'data','Demo_SAPHIRE_models.mat'),'SAPHIRE_Model');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Load SAPHIRE cell trajectory models
    load(fullfile(pwd,'data','Demo_SAPHIRE_models.mat'));

    %Plot the temporal evolution of a cell in state space, inferred states with
    %most likely state sequence, state annotation of cell image time series
    %frames, and state transition diagram for annotated sequence.
    model = SAPHIRE_Model{min(trajToAnalyze,numel(SAPHIRE_Model))};
    traj = CellData.shapeSpaceTraj{min(trajToAnalyze,numel(SAPHIRE_Model))};
    I = CellData.Iactin{min(trajToAnalyze,numel(SAPHIRE_Model))};
    mask = CellData.actinMask{min(trajToAnalyze,numel(SAPHIRE_Model))};
    plot_model_results(model,traj,I,mask)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Generate a composite phenotype dynamics profile from the temporal state
    %annotations of the cells. Note that this works for MATLAB version 2015
    %onward, since it uses the function histcounts. Users may change histogram
    %functions that are compatible with previous versions (e.g. hist, histc).
    rel = version('-release');
    rel = str2double(rel(1:4));
    if rel >= 2015
        [phenoComboDynamicFeats,phenoComboDynamicSig] = compute_pheno_signature(SAPHIRE_Model);
    end

end








