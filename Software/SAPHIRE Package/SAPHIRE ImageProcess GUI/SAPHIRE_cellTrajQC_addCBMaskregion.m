%Adds user-defined region to the foreground of cell body mask
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_addCBMaskregion(hObject,eventData)

dataGUI = guidata(hObject);  

%Disable button when one is active
set(dataGUI.btnHandleArray,'enable','off');

dataGUI.currAxis = axis;

%Define matrix where to store manual modifications of regions in mask.
modifiedRegions = zeros(size(dataGUI.editedCBMaskTile));
  
fprintf('\n%s\n','Add closed freehand region to mask...')
modifiedRegions = SAPHIRE_cellTrajQC_modifyRegions(dataGUI,...
    'add',true,modifiedRegions);

%Apply modifications to regions of the mask.
dataGUI.editedCBMaskTile(modifiedRegions == 1) = 0;
dataGUI.editedCBMaskTile(modifiedRegions == 2) = 1;
dataGUI.editedCBMaskTile(modifiedRegions == 3) = 0;

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



