%Modifies the phenotype label of a cell image frame based on the letter
%label specified by the user.
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_editFramePheno(hObject,eventData)

dataGUI = guidata(hObject);  

%Disable button when one is active
set(dataGUI.btnHandleArray,'enable','off');

dataGUI.currAxis = axis;

%User pressed button on location on editing frame phenotype
while 1
    [xTemp,yTemp,btnTemp] = myginput(1,'arrow');
    if xTemp < 1 || xTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,2) || yTemp < 1 || yTemp > size(dataGUI.rgbPaddedCellImgsTileIdx,1) || (dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp))==0) 
        disp('MUST CLICK WITHIN A CELL FRAME!')
        continue;
    else
        break;
    end
end

%Get frame for which pheno label was edited
frameEditPheno = dataGUI.rgbPaddedCellImgsTileIdx(round(yTemp),round(xTemp));

%Set the updated label
dataGUI.editedPhenoLabels{frameEditPheno} = char(btnTemp);

%Enable buttons again when done
set(dataGUI.btnHandleArray,'enable','on');

%Update GUI data:
guidata(hObject,dataGUI);

SAPHIRE_cellTrajQC_plotMain(hObject,eventData);
