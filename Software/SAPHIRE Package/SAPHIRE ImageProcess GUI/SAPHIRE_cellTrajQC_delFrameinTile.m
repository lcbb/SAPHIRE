%Deletes a user-specified image frame and places an "x" over the deleted
%frame.
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_delFrameinTile(hObject,eventData)

dataGUI = guidata(hObject);  

%Disable button when one is active
set(dataGUI.btnHandleArray,'enable','off');

dataGUI.currAxis = axis;

%User pressed left mouse on location on editing frame phenotype

while 1
    [xTemp,yTemp,~] = myginput(1,'hand');
    if xTemp < 1 || xTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,2) || yTemp < 1 || yTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,1) || (dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp))==0) 
        disp('MUST CLICK WITHIN A CELL FRAME!')
        continue;
    else
        break;
    end
end

%Get frame for which pheno label was edited
frameEditPheno = dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp));

if dataGUI.deleteFrameMask(frameEditPheno) == false
    dataGUI.deleteFrameMask(frameEditPheno) = true;
else
    dataGUI.deleteFrameMask(frameEditPheno) = false;
end

%Enable buttons again when done
set(dataGUI.btnHandleArray,'enable','on');

%Update GUI data:
guidata(hObject,dataGUI);

SAPHIRE_cellTrajQC_plotMain(hObject,eventData);
