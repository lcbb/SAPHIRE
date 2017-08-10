%Closes GUI and outputs the user modifications to the cellular time series
%as |dataGUI| variable to the workspace.
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function SAPHIRE_cellTrajQC_finished(hObject,eventData)

dataGUI = guidata(hObject);  

assignin('base','dataGUI',dataGUI) %update the variable data in the workspace

close all;

fprintf('\n%s\n\t','FINISHED EDITING IMAGE TIME SERIES OF CURRENT CELL');
    