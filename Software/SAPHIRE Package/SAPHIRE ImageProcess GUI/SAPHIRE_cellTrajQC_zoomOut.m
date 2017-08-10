%Zooms out the cellular time series tile images to original scale in the
%GUI window.
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_zoomOut(hObject,eventData)

dataGUI = guidata(hObject);  

%Reset original axis
dataGUI.currAxis = dataGUI.originalAxis;

%Update GUI data:
guidata(hObject,dataGUI);

SAPHIRE_cellTrajQC_plotMain(hObject,eventData);

zoom on;
end


