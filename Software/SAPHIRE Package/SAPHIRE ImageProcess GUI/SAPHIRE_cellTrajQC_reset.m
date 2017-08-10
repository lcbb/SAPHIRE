%Resets all cellular time series properties (phenotype labels, masks, etc.)
%to the original prior to user modifications.
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_reset(hObject,eventData)

dataGUI = guidata(hObject);  

dataGUI.deleteFrameMask = false(numel(dataGUI.deleteFrameMask),1);
dataGUI.editedCBMaskTile = dataGUI.pad_CBMaskTile;
dataGUI.editedPhenoLabels = dataGUI.origPhenoLabels;

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
 
%Update GUI data:
guidata(hObject,dataGUI);
