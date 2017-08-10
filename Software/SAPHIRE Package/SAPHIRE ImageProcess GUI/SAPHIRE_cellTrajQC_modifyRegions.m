%Modifies user-defined cell mask region
%        
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function modifiedRegions = SAPHIRE_cellTrajQC_modifyRegions(dataGUI,modifyType,closedStatus,modifiedRegions)
        
     h_ModifyRegion = imfreehand('Closed',closedStatus);
     setColor(h_ModifyRegion,'b');

     if closedStatus == false
         tempAll = round(h_ModifyRegion.getPosition);
         [~,index] = unique(tempAll,'rows');
         P0 = tempAll(sort(index),:);
         D = [0; cumsum(sum(abs(diff(P0)),2))]; % Need the distance between points...
         P = interp1(D,P0,D(1):.5:D(end)); % ...to close the gaps
         P = unique(round(P),'rows');
         idx = sub2ind(size(mask),P(:,2),P(:,1));
     elseif closedStatus == true
         temp = createMask(h_ModifyRegion);
         idx = find(temp==1);
     end

     if strcmp(modifyType,'split')
         %Use temp when splitting to make connected components.
         temp = zeros(size(dataGUI.editedCBMaskTile));
         temp(idx) = 1;
         temp = bwmorph(temp,'thin');
         temp = bwmorph(temp,'diag');
         idx = find(temp==1);
         modifiedRegions(idx) = 1;
     elseif strcmp(modifyType,'add')
         if closedStatus == false
             temp = zeros(size(dataGUI.editedCBMaskTile));
             temp(idx) = 1;
             temp = bwmorph(temp,'thin');
             temp = bwmorph(temp,'diag');
             idx = find(temp==1);
         end
         modifiedRegions(idx) = 2; %Added pixels = 2 
     elseif strcmp(modifyType,'remove')
         if closedStatus == false
             temp = zeros(size(dataGUI.editedCBMaskTile));
             temp(idx) = 1;
             temp = bwmorph(temp,'thin');
             temp = bwmorph(temp,'diag');
             idx = find(temp==1);
         end
         modifiedRegions(idx) = 3; %Removed pixels = 3  
     end

     while 1
         fprintf('\n%s\n','Press ''spacebar'' to finish or ''u'' to undo region modification.');
         pause
         currKey = get(dataGUI.hfig,'CurrentKey');
         if strcmp(currKey,'space') 
             fprintf('\n%s\n','Finished modifying object')
             break;
         end
         if sum(modifiedRegions(:))>0 && strcmp(currKey,'u') %undo current modification.
             modifiedRegions(idx) = 0;
             delete(h_ModifyRegion);
             break;
         end
     end
     
end

