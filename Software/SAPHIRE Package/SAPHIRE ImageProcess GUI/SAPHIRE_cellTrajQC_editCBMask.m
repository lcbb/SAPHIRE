% Author: Simon Gordonov 
% Date: 04/09/15

function RegressModel_cellShapeTrajectoryQC_editCBMask(hObject,eventData)

dataGUI = guidata(hObject);  

dataGUI.currAxis = axis;

%Plot initial segmentation.
plot_segmentation(dataGUI);

%%
%UPDATE THE DATA (EDITED SEGMENTATION MASK)
guidata(hObject,dataGUI);

%Add keypress listener to the segmentation figure.
set(dataGUI.hfig,'KeyPressFcn', @keypress_segmentation);
    
%Define matrix where to store manual modifications of regions in mask.
modifiedRegions = zeros(size(dataGUI.editedCBMaskTile));

displaySegCmds();

while 1

    figure(1)
    k = waitforbuttonpress;
        
    if k ~= 0
            
        currKey = get(dataGUI.hfig,'CurrentKey');
        currModifier = get(dataGUI.hfig,'CurrentModifier');
            
            switch currKey
   
                case 'a' %add freehand regions to mask
                    
                    if length(currModifier) == 1 && strcmp(currModifier{:},'control') %closed region
                        fprintf('\n%s\n','Add closed freehand region to mask...')
                        modifiedRegions = modifyRegions(dataGUI,...
                            'add',true,modifiedRegions);
                    elseif isempty(currModifier) %non-closed region
                        fprintf('\n%s\n','Add open freehand region to mask...')
                        modifiedRegions = modifyRegions(dataGUI,...
                            'add',false,modifiedRegions);
                    end
                    
                case 's' %split freehand regions in mask
                    
                        fprintf('\n%s\n','Split objects in mask...')
                        modifiedRegions = modifyRegions(dataGUI,...
                            'split',false,modifiedRegions);
                    
                case 'r' %remove freehand regions in mask
                    
                    if length(currModifier) == 1 && strcmp(currModifier{:},'control') %closed region
                        fprintf('\n%s\n','Remove closed freehand region from mask...')
                        modifiedRegions = modifyRegions(dataGUI,...
                            'remove',true,modifiedRegions);
                    elseif isempty(currModifier) %non-closed region
                        fprintf('\n%s\n','Remove open freehand region from mask...')
                        modifiedRegions = modifyRegions(dataGUI,...
                            'remove',false,modifiedRegions);
                    end
                    
                case 'escape'
                    fprintf('\n%s\n','FINISHING MASK EDITING')
                    break;
                    
                otherwise
                    fprintf('\n%s\n','Invalid key pressed!')
                    displaySegCmds();
                    continue;
            end
            
            %Apply modifications to regions of the mask.
            dataGUI.editedCBMaskTile(modifiedRegions == 1) = 0;
            dataGUI.editedCBMaskTile(modifiedRegions == 2) = 1;
            dataGUI.editedCBMaskTile(modifiedRegions == 3) = 0;

            %Plot segmentation.
           plot_segmentation(dataGUI);
           displaySegCmds()
               
        end
    end
    
    %Do final clean-up on masks:
    dataGUI.editedCBMaskTile = imreconstruct(dataGUI.pad_NucCenterMaskTile,dataGUI.editedCBMaskTile);
    dataGUI.editedCBMaskTile = imfill(dataGUI.editedCBMaskTile,'holes'); 
    
    plot_segmentation(dataGUI);
    
    %Update GUI data:
    guidata(hObject,dataGUI);
    
    pan on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displaySegCmds()

    fprintf('\n%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n',...
        'Enter a command...',...
        'a - add open region to mask',...
        'ctrl+a - add closed region to mask',...
        'r - remove open region in mask',...
        'ctrl+r - remove closed region in mask',...
        'ESC - end segmentation module');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keypress_segmentation(~,~)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modifiedRegions = modifyRegions(dataGUI,modifyType,closedStatus,modifiedRegions)
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_segmentation(dataGUI)

        imshowpair(dataGUI.pad_CBTile,dataGUI.pad_NucTile)
        hold on;
        boundPts = bwboundaries(dataGUI.editedCBMaskTile);
        for i = 1:numel(boundPts)
            plot(boundPts{i}(:,2),boundPts{i}(:,1),'color','r','linewidth',0.5);
            hold on;
        end
        set(gcf,'position',get(0,'screensize'));
        axis(dataGUI.currAxis)
        hold off;   
        
end




