%Deletes the entire cellular time series
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_delEntireTraj(hObject,eventData)

dataGUI = guidata(hObject);  

%Disable button when one is active
set(dataGUI.btnHandleArray,'enable','off');

dataGUI.currAxis = axis;

%Set to delete all frames
locNoDelete = ~dataGUI.deleteFrameMask;
dataGUI.deleteFrameMask(locNoDelete) = true;

%Check that user really wants to delete all frames
while 1
    fprintf('\n%s\n','Press ''spacebar'' to delete entire trajectory or ''u'' to undo.');
    pause
    currKey = get(dataGUI.hfig,'CurrentKey');
    if strcmp(currKey,'space') 
         fprintf('\n%s\n','Finished deleting entire trajectory')
         break;
    end
    if strcmp(currKey,'u') %undo
         dataGUI.deleteFrameMask(locNoDelete) = false;
         break;
    end
end

%Enable buttons again when done
set(dataGUI.btnHandleArray,'enable','on');

%Update GUI data:
guidata(hObject,dataGUI);

SAPHIRE_cellTrajQC_plotMain(hObject,eventData);
