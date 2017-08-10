%Plots modified cellular time series images in GUI window
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_plotMain(hObject,eventData)

dataGUI = guidata(hObject);  

% dataGUI.currAxis = axis;

imshowpair(dataGUI.pad_CBTile,dataGUI.pad_NucTile)
hold on;
boundPts = bwboundaries(dataGUI.editedCBMaskTile);
for i = 1:numel(boundPts)
    plot(boundPts{i}(:,2),boundPts{i}(:,1),'color','r','linewidth',0.5);
    hold on;
end

%Plot x's where frames are deleted
for i = 1:numel(dataGUI.deleteFrameMask)
    if dataGUI.deleteFrameMask(i) ~= 0
        [rowDel,colDel] = find(dataGUI.pad_NucCenterMaskTile & (dataGUI.rgbPaddedCellImgsTileIdx == i));
        plot(colDel,rowDel,'bx','markersize',30,'linewidth',3);
    end
end

%Plot the text labels of cell phenotype at each frame
for i = 1:numel(dataGUI.editedPhenoLabels)
    text(dataGUI.rgbPaddedCellImgsTileCornerCoords(i,3),...
        dataGUI.rgbPaddedCellImgsTileCornerCoords(i,1),...
        dataGUI.editedPhenoLabels{i},'color','w');
end

set(gcf,'position',get(0,'screensize'));
axis(dataGUI.currAxis)
hold off;   

pan on
 
%Update GUI data:
guidata(hObject,dataGUI);
