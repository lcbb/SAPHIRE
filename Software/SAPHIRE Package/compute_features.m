%%%%%%%%%%%%%%%%%%%%
% Computes image-based features from time series cell objects
%
% Inputs:   Idata - structure array of temporal trajectories of individual
%                   cell body image and their binary masks
%           pxScale - scale of each pixel (e.g. microns per pixel from the
%           imaging experiment)
%           featLabels - cell array of names of features to compute       
% Outputs:  trajSeriesFeatures - structure of features computed for input
%                                cell objects
%               
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function trajSeriesFeatures = compute_features(Idata,pxScale,featLabels)

disp('Computing features...');

numTrajs = numel(Idata.actinMask);

for trajNum = 1:numTrajs
    for featNum = 1:numel(featLabels)
        eval(['trajSeriesFeatures.',featLabels{featNum},'{trajNum,1} = zeros(numel(Idata.Iactin{trajNum}),1);']);
    end
end

for trajNum = 1:numTrajs
    for frameNum = 1:numel(Idata.Iactin{trajNum})
        for featNum = 1:numel(featLabels)
            eval(['trajSeriesFeatures.',featLabels{featNum},'{trajNum,1}(frameNum) = wcf_',featLabels{featNum},...
                '(Idata.Iactin{trajNum}{frameNum},full(Idata.actinMask{trajNum}{frameNum}),pxScale);']);
        end
    end
    disp(['Feature calculations for trajectory ',num2str(trajNum),' of ',num2str(numTrajs),' complete.']);
end





