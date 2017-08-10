%User-based editing of cell segmentation masks
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_editSegMaskReg(hObject,eventData)

dataGUI = guidata(hObject);  

%Disable button when one is active
set(dataGUI.btnHandleArray,'enable','off');

dataGUI.currAxis = axis;

%User pressed button on location on editing frame's segmentation thresholds
while 1
    [xTemp,yTemp] = myginput(1,'hand');
    if xTemp < 1 || xTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,2) || yTemp < 1 || yTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,1) || (dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp))==0) 
        disp('MUST CLICK WITHIN A CELL FRAME!')
        continue;
    else
        break;
    end
end

%Get the frame image on which to update the segmentation parameters
clickedFrameIdx = dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp));
Icb = mat2gray(dataGUI.pad_CBTile(dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,1):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,2),...
    dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,3):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,4)));
Inuc = mat2gray(dataGUI.pad_NucTile(dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,1):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,2),...
    dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,3):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,4)));
InucMask = dataGUI.pad_NucCenterMaskTile(dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,1):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,2),...
    dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,3):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,4));
initialTotsuCB = graythresh(Icb)-0.4*(graythresh(Icb));
[~,initialTedgeCB] = edge(Icb,'canny');
                
%Load the interactive segmentation threshold setting GUI to get the mask
%for the cell in chosen frame.
cbMask = SAPHIRE_cellTrajQC_SegThreshSet(Icb,Inuc,InucMask,initialTotsuCB,initialTedgeCB,dataGUI.paramsMain);
dataGUI.editedCBMaskTile(dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,1):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,2),...
    dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,3):dataGUI.rgbPaddedCellImgsTileImgCornerCoords(clickedFrameIdx,4)) = cbMask;

%Make sure any added area outside each cropped frame's boundary isn't
%included.
dataGUI.editedCBMaskTile = dataGUI.pad_frameFillTile & dataGUI.editedCBMaskTile;

%Do final clean-up on masks:
dataGUI.editedCBMaskTile = imreconstruct(dataGUI.pad_NucCenterMaskTile,dataGUI.editedCBMaskTile);
dataGUI.editedCBMaskTile = imfill(dataGUI.editedCBMaskTile,'holes'); 

%Enable buttons again when done
set(dataGUI.btnHandleArray,'enable','on');

%Update GUI data:
guidata(hObject,dataGUI);

SAPHIRE_cellTrajQC_plotMain(hObject,eventData);

end


