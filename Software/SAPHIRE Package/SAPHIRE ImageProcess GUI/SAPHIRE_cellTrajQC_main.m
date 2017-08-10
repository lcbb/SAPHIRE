%GUI for editing cellular time series images by the user. An individual
%cell trajectory image frames are shown to the user as a montage over time,
%and the user selects to delete certain frames, add/remove
%foreground/background regions, modify segmentation automatically with new
%parameters, re-label cell object phenotype (e.g. dividing or dying cells),
%to generate reliable cell trajectories for subsequent feature extraction
%and time series modeling.
%
% Input: Iarray.(CB,CBMask,nucCenters,Nuc) - structure containing cell
%       arrays of cell body images, cell body masks, nuclear centers, and
%       nuclear images.
%       Iarray.phenoLabels - cell array specifying phenotype labels of
%       cells as done originally by automated analysis/labeling.
%
% Outputs: dataGUI - GUI data of user-based modifications are stored in the
%                    workspace.
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_main(Iarray)

    close all;

    %Set parameters for cell body segmentation modifications
    Iarray.params.cbSeg.tophatDiskR = 40; %size of disk for background subtraction
    Iarray.params.cbSeg.otsuCorrectFrac = 0.01; %fraction of Otsu threshold to iterate by
    Iarray.params.cbSeg.edgeCorrectFrac = 0.01; %fraction of Canny edge detection thresholds to iterate by
    Iarray.params.cbSeg.minObjSize = 100; % minimum cell object size (those below this value are removed)
    
    fprintf('\n%s\n','STARTING CELLULAR IMAGE TIME SERIES EDITING TOOL')

    %First we must pad all cell images to be the same size (based on the
    %biggest cell object in the time series) so that montage doesn't skew
    %the size of the cells during display.

    %Get the largest image size:
    maxVertSize = max(cell2mat(cellfun(@(x) size(x,1),Iarray.CBMask,'uniformoutput',0))); 
    maxHorizSize = max(cell2mat(cellfun(@(x) size(x,2),Iarray.CBMask,'uniformoutput',0)));

    %Set number of frames
    numFrames = numel(Iarray.CB);
    
    %Get the bit depth of the images
    bitDepthNuc = class(Iarray.Nuc{1});
    bitDepthCB = class(Iarray.CB{1});
    
    %Loop through all cell images, mask the nuclear and cell body channels,
    %pad them with zeros to be maxVertSize x maxHorizSize, and store them
    %in an array.
    rgbPaddedCellImgs = cell(numFrames,1);
    pad_CBarray = cell(numFrames,1);
    pad_CBMaskarray = cell(numFrames,1);
    pad_Nucarray = cell(numFrames,1);
    pad_NucCenterMaskarray = cell(numFrames,1);
    pad_frameFill = cell(numFrames,1); %ones where each cropped frame is (rectangle of ones)
    for i = 1:numFrames
        Inuc = Iarray.Nuc{i}; 
        Icb = Iarray.CB{i};
        IcbMask = Iarray.CBMask{i};
        InucCentMask = Iarray.NucCentMask{i};
        %Pad the nucleus and cell body images to be the max size of all images
        [VertSize HorizSize] = size(Inuc);
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
        pad_Nucarray{i} = padarray(Inuc,[padVertPre padHorizPre],'pre');
        pad_Nucarray{i} = padarray(pad_Nucarray{i},[padVertPost padHorizPost],'post');
        pad_CBarray{i} = padarray(Icb,[padVertPre padHorizPre],'pre');
        pad_CBarray{i} = padarray(pad_CBarray{i},[padVertPost padHorizPost],'post');
        pad_CBMaskarray{i} = padarray(IcbMask,[padVertPre padHorizPre],'pre');
        pad_CBMaskarray{i} = padarray(pad_CBMaskarray{i},[padVertPost padHorizPost],'post');
        pad_NucCenterMaskarray{i} = padarray(InucCentMask,[padVertPre padHorizPre],'pre');
        pad_NucCenterMaskarray{i} = padarray(pad_NucCenterMaskarray{i},[padVertPost padHorizPost],'post');
        pad_frameFill{i} = padarray(true(size(Inuc)),[padVertPre padHorizPre],'pre');
        pad_frameFill{i} = padarray(pad_frameFill{i},[padVertPost padHorizPost],'post');
        %Combine images and make RGB:
        rgbPaddedCellImgs{i} = cat(3,zeros(maxVertSize,maxHorizSize),...
            pad_CBarray{i},imdilate(255*pad_NucCenterMaskarray{i},ones(10)));
    end

    %Set the dimensions of the large montage/tile image of frames to make
    %the tile array as square as possible
    newDims = sqrt(numFrames);
    if newDims - floor(newDims) ~= 0
        dim = ceil(newDims);
        newDims = [dim dim-1];
        residDim = prod(newDims) - numFrames;
        eval(['rgbPaddedCellImgs = [rgbPaddedCellImgs;repmat({',bitDepthCB,'(zeros(maxVertSize,maxHorizSize,3))},residDim,1)];']);
        pad_NucCenterMaskarray = [pad_NucCenterMaskarray;repmat({false(maxVertSize,maxHorizSize)},residDim,1)];
        pad_CBMaskarray = [pad_CBMaskarray;repmat({false(maxVertSize,maxHorizSize)},residDim,1)];
        eval(['pad_CBarray = [pad_CBarray;repmat({',bitDepthCB,'(zeros(maxVertSize,maxHorizSize))},residDim,1)];']);
        eval(['pad_Nucarray = [pad_Nucarray;repmat({',bitDepthNuc,'(zeros(maxVertSize,maxHorizSize))},residDim,1)];']);
        pad_frameFill = [pad_frameFill;repmat({false(maxVertSize,maxHorizSize)},residDim,1)];
    else
        newDims = [newDims newDims];
    end
    
    %Create a matrix that stores the indices of each frame's pixels in the
    %larger montage matrix, and also reshape the images into a large
    %matrix.
    rgbPaddedCellImgsIdx = repmat({zeros(maxVertSize,maxHorizSize)},numel(rgbPaddedCellImgs),1);
    for i = 1:numFrames
        rgbPaddedCellImgsIdx{i} = repmat(i,maxVertSize,maxHorizSize);
    end
    
    %Reshape all needed cell arrays into tile image
    rgbPaddedCellImgsTile = cell2mat(reshape(rgbPaddedCellImgs,newDims(2),newDims(1))');
    dataGUI.rgbPaddedCellImgsTileIdx = cell2mat(reshape(rgbPaddedCellImgsIdx,newDims(2),newDims(1))');
    dataGUI.pad_NucCenterMaskTile = cell2mat(reshape(pad_NucCenterMaskarray,newDims(2),newDims(1))'); 
    dataGUI.pad_CBMaskTile = cell2mat(reshape(pad_CBMaskarray,newDims(2),newDims(1))'); 
    dataGUI.pad_CBTile = cell2mat(reshape(pad_CBarray,newDims(2),newDims(1))'); 
    dataGUI.pad_NucTile = cell2mat(reshape(pad_Nucarray,newDims(2),newDims(1))'); 
    dataGUI.pad_frameFillTile = cell2mat(reshape(pad_frameFill,newDims(2),newDims(1))'); 
    
    %Get the coordinates (x,y) of each frame's corners
    dataGUI.rgbPaddedCellImgsTileCornerCoords = zeros(numFrames,4);
    for i = 1:numFrames
        [xTemp,yTemp] = find(dataGUI.rgbPaddedCellImgsTileIdx == i);
        dataGUI.rgbPaddedCellImgsTileCornerCoords(i,:) = [min(xTemp),max(xTemp),...
            min(yTemp),max(yTemp)];
    end
    
    %Get the coordinates (x,y) of each frame's image corners. Note that
    %this is different from above in that we here extract the corners of
    %the cell IMAGE within each time, not the padded tile's corners.
    dataGUI.rgbPaddedCellImgsTileImgCornerCoords = zeros(numFrames,4);
    for i = 1:numFrames
        [xTemp,yTemp] = find(dataGUI.rgbPaddedCellImgsTileIdx == i & dataGUI.pad_frameFillTile);
        dataGUI.rgbPaddedCellImgsTileImgCornerCoords(i,:) = [min(xTemp),max(xTemp),...
            min(yTemp),max(yTemp)];
    end
    
    %Get the cell body mask object boundary coordinates in the large tile
    %image
    rgbPaddedCellImgsTileBoundCoords = cell(numFrames,1);
    for i = 1:numFrames
       rgbPaddedCellImgsTileBoundCoords{i} = cell2mat(bwboundaries(pad_CBMaskarray{i}));
       %Add on the shift of coordinates in tile image
       rgbPaddedCellImgsTileBoundCoords{i}(:,1) = rgbPaddedCellImgsTileBoundCoords{i}(:,1) + dataGUI.rgbPaddedCellImgsTileCornerCoords(i,1) - 1;
       rgbPaddedCellImgsTileBoundCoords{i}(:,2) = rgbPaddedCellImgsTileBoundCoords{i}(:,2) + dataGUI.rgbPaddedCellImgsTileCornerCoords(i,3) - 1;
    end
       
    %%  IMAGE TIME SERIES EDITING GUI

    dataGUI.newDims = newDims;
    dataGUI.hfig = figure;
    dataGUI.deleteFrameMask = false(numFrames,1);
    dataGUI.editedCBMaskTile = dataGUI.pad_CBMaskTile;
    dataGUI.origPhenoLabels = Iarray.phenoLabels;
    dataGUI.editedPhenoLabels = dataGUI.origPhenoLabels;
    dataGUI.cbImgTile = rgbPaddedCellImgsTile(:,:,2);
    dataGUI.paramsMain = Iarray.params;

    %Plot all the cell outlines for shape trajectory and assign variables
    %to the segmentation outlines

    imshowpair(dataGUI.pad_CBTile,dataGUI.pad_NucTile)
    hold on;
    boundPts = bwboundaries(dataGUI.editedCBMaskTile);
    for i = 1:numel(boundPts)
        plot(boundPts{i}(:,2),boundPts{i}(:,1),'color','r','linewidth',0.5);
        hold on;
    end
    %Plot the text labels of cell phenotype at each frame
    for i = 1:numel(dataGUI.editedPhenoLabels)
        text(dataGUI.rgbPaddedCellImgsTileCornerCoords(i,3),...
            dataGUI.rgbPaddedCellImgsTileCornerCoords(i,1),...
            dataGUI.editedPhenoLabels{i},'color','w');
    end
    hold off;
    zoom on
    set(gcf,'toolbar','figure');
    set(gcf,'color','k')
    dataGUI.originalAxis = axis;

    rightShift = 600;
    set(gcf,'position',get(0,'screensize'))

    dataGUI.EditSegThreshBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Edit Seg Thresh',...
        'Position', [rightShift 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_editSegMaskReg(hObject,eventData));

    dataGUI.AddSegBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Add Seg',...
        'Position', [rightShift+(110*1) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_addCBMaskregion(hObject,eventData));

    dataGUI.RemSegBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Remove Seg',...
        'Position', [rightShift+(110*2) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_remMaskReg(hObject,eventData));

    dataGUI.DelFrameBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Delete Frame',...
        'Position', [rightShift+(110*3) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_delFrameinTile(hObject,eventData));

    dataGUI.DelAllFrames = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Delete All',...
        'Position', [rightShift+(110*4) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_delEntireTraj(hObject,eventData));

    dataGUI.phenoEditBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Pheno Edit',...
        'Position', [rightShift+(110*5) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_editFramePheno(hObject,eventData));

    dataGUI.finishedBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Finish',...
        'Position', [rightShift+(110*6) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_finished(hObject,eventData));

    dataGUI.resetBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Reset',...
        'Position', [rightShift+(110*7) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_reset(hObject,eventData));

    dataGUI.zoomOutBtn = uicontrol('Parent',dataGUI.hfig,'Style', 'pushbutton', 'String', 'Zoom Out',...
        'Position', [rightShift+(110*8) 50 100 20],...
        'Callback',@(hObject,eventData) SAPHIRE_cellTrajQC_zoomOut(hObject,eventData));

    %Array that stores all buttons
    dataGUI.btnHandleArray = [dataGUI.zoomOutBtn,dataGUI.EditSegThreshBtn,dataGUI.AddSegBtn,dataGUI.DelFrameBtn,dataGUI.RemSegBtn,dataGUI.finishedBtn,dataGUI.phenoEditBtn,dataGUI.resetBtn,dataGUI.DelAllFrames];

    guidata(dataGUI.hfig,dataGUI)

    waitfor(dataGUI.hfig) %wait for gui to close



   
    
    
         