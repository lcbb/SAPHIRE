%%%%%%%%%%%%%%%%%%%%
% Plots results of SAPHIRE analysis on an individual cell trajectory
%
% Inputs:   model - SAPHIRE model strucutre of the cell trajectory
%           traj - dxT matrix of coordinates of the trajectory in
%                  shape-space, where d is the coordinate dimensionality
%                  (e.g. d = 2 for 2-D PC space) and T is the trajectory 
%                  length.
%           I - cell array of image frames of the fluorescence channel
%               images of the cell body
%           mask - cell array of image frames of the binary masks of the 
%                  cell body
%               
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function plot_model_results(model,traj,I,mask)

%Get the number of time points in the cell trajectory and number of states
numTimePts = size(traj,2);
numStates = size(model.ML_params.mu_emit,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figDims = get(0,'screensize');

%Plot the points in the cell's trajectory in state space in addition to the
%inferred states + sigma and 2*sigma.
figure;
set(gcf,'color','w');
set(gcf,'position',[295 figDims(end)-420 270 340]);

%Plot the maximum likelihood state annotations of time series points and
%underlying Gaussian states
cmap = distinguishable_colors(numel(unique(model.ML_states)));
for i = 1:numTimePts-1
    plot([traj(1,i),traj(1,i+1)],[traj(2,i),traj(2,i+1)],...
        '-','color',[0.3 0.3 0.3],'linewidth',1);
    hold on;
end
for i = 1:numTimePts
    plot(traj(1,i),traj(2,i),...
        'o','markerfacecolor',cmap(model.ML_states(i),:),'markeredgecolor',...
        'w','markersize',8);
    hold on;
end
for i = 1:max(model.ML_states)
    mu = model.ML_params.mu_emit(:,i);
    sigma = model.ML_params.sigma_emit(i);
    a = 0:2*pi/100:2*pi;
    line(sigma*cos(a)+mu(1), sigma*sin(a)+mu(2),'color','k','linewidth',1);
    hold on;
    line(2*sigma*cos(a)+mu(1),2*sigma*sin(a)+mu(2),'color','k','linewidth',1);
end
xLim = get(gca,'xlim');
yLim = get(gca,'ylim');
hold on;
xlabel('PC1')
ylabel('PC2')
set(gca,'xlim',[min(xLim(1),yLim(1)),max(xLim(2),yLim(2))],'ylim',[min(xLim(1),yLim(1)),max(xLim(2),yLim(2))]);
daspect([1 1 1])
colormap(cmap);
cA = colorbar('southoutside');
cA.TickLabelsMode = 'manual';
cA.TicksMode = 'manual';
cA.Ticks = 1/(2*numStates):1/numStates:1-(1/(2*numStates));
cA.TickLabels = 1:numStates;
cA.Label.String = 'State Color Label';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot temporal evolution of cell trajectory
figure;
set(gcf,'color','w');
set(gcf,'position',[15 figDims(end)-420 270 340]);

cmap = labmap([0,0],[10 90],[0 0]);
cmap = cmap(round(linspace(1,size(cmap,1),numTimePts-1)),:);
for i = 1:numTimePts-1
    plot([traj(1,i),traj(1,i+1)],[traj(2,i),traj(2,i+1)],...
        '-','color',cmap(i,:),'linewidth',2);
    hold on;
end
xlabel('PC1')
ylabel('PC2')
set(gca,'xlim',[min(xLim(1),yLim(1)),max(xLim(2),yLim(2))],'ylim',[min(xLim(1),yLim(1)),max(xLim(2),yLim(2))]);
daspect([1 1 1])
colormap(cmap);
cB = colorbar('southoutside');
cB.Ticks = [];
cB.TickLabels = '';
cB.Label.String = 'Time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot frame images of the cell across the trajectory colord as the most
%likely states

mask = cellfun(@(x) full(x),mask,'uniformoutput',0);

%First we must pad all cell images to be the same size (based on the
%biggest cell) so that montage doesn't skew the size of the cells during
%display.

%Get the largest image size:
maxVertSize = max(cell2mat(cellfun(@(x) size(x,1),mask,'uniformoutput',0))); 
maxHorizSize = max(cell2mat(cellfun(@(x) size(x,2),mask,'uniformoutput',0)));

%Set number of frames
numFrames = numel(I);

%Get the bit depth of the image
bitDepthCB = class(I{1});

%Loop through all cell images, pad them with zeros to be maxVertSize x
%maxHorizSize, and store them in an array.
pad_CBarray = cell(numFrames,1);
pad_CBMaskarray = cell(numFrames,1);
for i = 1:numFrames
    Icb = I{i};
    IcbMask = mask{i};
    %Pad the images to be the max size of all images
    [VertSize HorizSize] = size(Icb);
    vertPad = maxVertSize - VertSize;
    if mod(vertPad,2) ~= 0 %if odd number of pixels to pad
        padVertPre = ceil(vertPad/2);
        padVertPost = floor(vertPad/2);
    else %if even number of pixels to pad, make symmetric around the image
        padVertPre = vertPad/2;
        padVertPost = vertPad/2;
    end
    horizPad = maxHorizSize - HorizSize;
    if mod(horizPad,2) ~= 0 %if odd number of pixels to pad
        padHorizPre = ceil(horizPad/2);
        padHorizPost = floor(horizPad/2);
    else %if even number of pixels to pad, make symmetric around the image
        padHorizPre = horizPad/2;
        padHorizPost = horizPad/2;
    end
    %Pad the images:
    pad_CBarray{i} = padarray(Icb,[padVertPre padHorizPre],'pre');
    pad_CBarray{i} = padarray(pad_CBarray{i},[padVertPost padHorizPost],'post');
    pad_CBMaskarray{i} = padarray(IcbMask,[padVertPre padHorizPre],'pre');
    pad_CBMaskarray{i} = padarray(pad_CBMaskarray{i},[padVertPost padHorizPost],'post');
end

%Set the dimensions of the large montage/tile image of frames to make
%the tile array as square as possible
newDims = sqrt(numFrames);
if newDims - floor(newDims) ~= 0
    dim = ceil(newDims);
    newDims = [dim dim];
    residDim = prod(newDims) - numFrames;
    pad_CBMaskarray = [pad_CBMaskarray;repmat({false(maxVertSize,maxHorizSize)},residDim,1)];
    eval(['pad_CBarray = [pad_CBarray;repmat({',bitDepthCB,'(zeros(maxVertSize,maxHorizSize))},residDim,1)];']);
else
    newDims = [newDims newDims];
end

%Create a matrix that stores the indices of each frame's pixels in the
%larger montage matrix, and also reshape the images into a large
%matrix.
rgbPaddedCellImgsIdx = repmat({zeros(maxVertSize,maxHorizSize)},numel(pad_CBMaskarray),1);
for i = 1:numFrames
    rgbPaddedCellImgsIdx{i} = repmat(i,maxVertSize,maxHorizSize);
end

%Reshape all needed cell arrays into tile image
data.rgbPaddedCellImgsTileIdx = cell2mat(reshape(rgbPaddedCellImgsIdx,newDims(2),newDims(1))');
data.pad_CBMaskTile = cell2mat(reshape(pad_CBMaskarray,newDims(2),newDims(1))'); 
data.pad_CBTile = cell2mat(reshape(pad_CBarray,newDims(2),newDims(1))'); 

%Plotting
figure
cmap = distinguishable_colors(max(model.ML_states));
data.pad_CBTile(~data.pad_CBMaskTile) = 0;
imshow(imcomplement(data.pad_CBTile),[])
hold on;
boundPts = bwboundaries(data.pad_CBMaskTile);
for i = 1:numel(boundPts)
    plot(boundPts{i}(:,2),boundPts{i}(:,1),'color',cmap(model.ML_states(data.rgbPaddedCellImgsTileIdx(boundPts{i}(1,1),boundPts{i}(1,2))),:),'linewidth',2);
    hold on;
    [r,c] = find(data.rgbPaddedCellImgsTileIdx == i);
    r = min(r);c = min(c);
    text(c,r,num2str(model.ML_states(i)),'fontsize',8);
end
title(sprintf('%s\n%s','ML temporal state annotation sequence of cell trajectory','(time = left to right, then top to bottom)'),'fontsize',10)
set(gcf,'color','w');
tightfig();
set(gcf,'position',[figDims(3)-450,60,440,720]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot graph with nodes as states colored by state dwell frequency and edges
%state transitions, colored by transitioned-to state dwell time and
%thickness the transition magnitude.

if numStates > 1

    %State locations
    stateDwellFreqs = zeros(numStates,1);
    for i = 1:numStates
        stateDwellFreqs(i) = sum(model.ML_states==i)/numTimePts;
    end

    %State transitions
    stateTransDwellFreqs = zeros(numStates,numStates);
    stateTransMags = zeros(numStates,numStates);

    stateSeq = model.ML_states;
    locTrans = find(abs(stateSeq(2:end)-stateSeq(1:end-1))); %find locations where state transitions occur (time just before transition)
    temp = zeros(1,length(stateSeq));
    temp(locTrans+1) = 1;
    idx = cumsum(temp)+1;
    %Partition the states over time into cell array
    statePartitions = accumarray(idx',stateSeq',[],@(x) {x});
    statePartitionsUnique = cell2mat(cellfun(@(x) unique(x),statePartitions,'uniformoutput',0));
    statePartitionsNums = cell2mat(cellfun(@(x) numel(x),statePartitions,'uniformoutput',0));

    for i = 1:numStates
        for j = 1:numStates
            if i == j
                continue;
            else
                stateTransMags(i,j) = norm(model.ML_params.mu_emit(:,i)'-model.ML_params.mu_emit(:,j)');
                transVec = [i,j];
                loc = find(all(ismember([statePartitionsUnique(1:end-1),statePartitionsUnique(2:end)],transVec,'rows'),2))+1;
                stateTransDwellFreqs(i,j) = sum(statePartitionsNums(loc))/numel(stateSeq);
            end
        end
    end
    stateTransMags(stateTransDwellFreqs==0) = 0;

    %Make network diagram
    nodeLabels = cellfun(@(x) num2str(x),num2cell(1:numStates),'uniformoutput',0);
    nodeLabels = cellfun(@(x) ['State ',x],nodeLabels,'uniformoutput',0);
    stateNet = biograph(stateTransMags>0,nodeLabels,'Layouttype','hierarchical','NodeAutoSize','off','showweights','off','arrowsize',12);

    cmapDwellFreq = bilateralmap(0,0.07);
    cmapDwellFreq = cmapDwellFreq(160:256,:);
    colorRangeMaxDwell = max(stateDwellFreqs);
    colorScaleDwellFreq = linspace(0,colorRangeMaxDwell,size(cmapDwellFreq,1));
    % cmap = distinguishable_colors(numel(unique(model.ML_states)));
    for i = 1:numStates
        stateNet.nodes(i).Label = nodeLabels{i};
        stateNet.nodes(i).Shape = 'circle';
        stateNet.nodes(i).Size = [50 50];
        stateNet.nodes(i).FontSize = 12;
        stateNet.nodes(i).Color = interp1(colorScaleDwellFreq,cmapDwellFreq,stateDwellFreqs(i));
        stateNet.nodes(i).LineColor = [0 0 0];%cmap(i,:);
        stateNet.nodes(i).LineWidth = 1;
    end

    cmapTransFreq = bilateralmap(0,0.6);
    cmapTransFreq = cmapTransFreq(160:256,:);
    colorRangeMaxTrans = max(stateTransDwellFreqs(:));
    colorScaleTransDwellFreq = linspace(0,colorRangeMaxTrans,size(cmapTransFreq,1));
    stateTransDwellFreqsReshape = reshape(stateTransDwellFreqs',numStates^2,1);
    stateTransDwellFreqsReshape = stateTransDwellFreqsReshape(stateTransDwellFreqsReshape~=0);
    %     %Line widths based on magnitude of transitions
    %     stateTransMagsReshape = reshape(stateTransMags',numStates^2,1);
    %     stateTransMagsReshape = stateTransMagsReshape(stateTransMagsReshape~=0);
    %     scaleTransMags = linspace(min(stateTransMagsReshape),max(stateTransMagsReshape),60);
    %     lineWidthTransVec = linspace(1,6,60);
    for i = 1:numel(stateNet.edges)
        if numel(unique(stateTransDwellFreqsReshape)) == 1
            stateNet.edges(i).LineColor = cmapTransFreq(end,:);
        else
            stateNet.edges(i).LineColor = interp1(colorScaleTransDwellFreq,cmapTransFreq,stateTransDwellFreqsReshape(i));
        end
    %         if numel(unique(stateTransMagsReshape)) == 1
    %             stateNet.edges(i).LineWidth = lineWidthTransVec(1);
    %         else
    %             stateNet.edges(i).LineWidth = interp1(scaleTransMags,lineWidthTransVec,stateTransMagsReshape(i));
    %         end
    %     stateNet.edges(i).Weight = round(100*stateTransDwellFreqsReshape(i))/100;
    end


    hNet = view(stateNet);
    g = biograph.bggui(hNet);
    f = get(g.biograph.hgAxes,'Parent');
    set(f,'position',[35 60 345 280])
    set(f,'name','State Network')
    child_handles = allchild(0);
    names = get(child_handles,'Name');
    close(child_handles(find(strncmp('Biograph', names, 8))))
    figPos = get(f,'position');

    figure;
    axis off;
    colormap(cmapDwellFreq);
    caxis([0 colorRangeMaxDwell])
    c1 = colorbar('north');
    c1.Label.String = 'Normalized State Dwell Time';
    set(gcf,'position',[figPos(1)+figPos(3)+140 figPos(2)+figPos(4)-30 240 140]);
    set(gcf,'color','w');
    set(gcf,'name','Normalized State Dwell Time');
    
    figure;
    axis off;
    colormap(cmapTransFreq);
    caxis([0 colorRangeMaxTrans])
    c2 = colorbar('north');
    c2.Label.String = 'Normalized State Transition Dwell Time';
    set(gcf,'position',[0 0 240 140]);
    set(gcf,'position',[figPos(1)+figPos(3)+140 figPos(2)+figPos(4)-240 240 140]);
    set(gcf,'color','w');
    set(gcf,'name','Normalized State Transition Dwell Time');

end




