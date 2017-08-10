%Performs cell body segmentation based on input parameters
%
% Input: Icb - cell body fluorescence image
%        Inuc - cell nucleus flourescence image
%        InucMask - nucleus mask center
%        initialTotsuCB - initial Otsu threshold 
%        initialTedgeCM - initial, 2-element, vector of Canny edge thresholds
%        params - parameters for changing segmentation   
%
% Outputs: cbMask - updated cell body mask
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%

function cbMask = SAPHIRE_cellTrajQC_SegThreshSet(Icb,Inuc,InucMask,initialTotsuCB,initialTedgeCB,params)

    %Refine the nuclei centers
        maskNucShrink = bwmorph(InucMask,'shrink','inf');
    
    %Perform background subtraction using tophat:
        if params.cbSeg.tophatDiskR ~= 0
            Icb = imtophat(Icb,strel('disk',params.cbSeg.tophatDiskR));
        end

    %Initialize the segmentation global threshold and edge threshold:
        Totsu = initialTotsuCB; %Otsu modified threshold
        Tedge = initialTedgeCB; %Canny edge parameters
        
    %Initialize faster (non-final) segmentation.
        finalBool = false;
        
    %Define matrix where to store manual modifications of regions in mask.
        modifiedRegions = zeros(size(Icb));
    
    %Run initial segmentation.
        cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
        
    %Plot initial segmentation.
        hf_segment = figure(2);
        set(gcf,'color','k');
   
    %Get nuclei center coordinates for plotting and plot segmentation +
    %nuclei from nuclear segmentation.
        [tempR,tempC] = ind2sub(size(Icb),...
            find(full(maskNucShrink)));
        plot_segmentation(Icb,Inuc,cbMask,tempR,tempC);

    %Add keypress listener to the segmentation figure.
        set(hf_segment,'KeyPressFcn', @keypress_segmentation);

    while 1
        
        displaySegCmds()
        
        k = waitforbuttonpress;
        
        if k ~= 0
            
            currKey = get(hf_segment,'CurrentKey');
            currModifier = get(hf_segment,'CurrentModifier');
            
            switch currKey
                
                case 'uparrow'
                    
                    Totsu = adjust_otsuThresh(Totsu,params,'+');
                    disp(['New Totsu = ',num2str(Totsu)]);
                    %Run segmentation with new thresholds.
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);

                case 'downarrow'
                    
                    Totsu = adjust_otsuThresh(Totsu,params,'-');
                    disp(['New Totsu = ',num2str(Totsu)]);
                    %Run segmentation with new thresholds.
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                    
                case 'leftarrow'
                    
                    Tedge = adjust_edgeThresh(Tedge,params,'-');
                    disp(['New Tedge = ',num2str(Tedge)]);
                    %Run segmentation with new thresholds.
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                    
                case 'rightarrow' 
                    
                    Tedge = adjust_edgeThresh(Tedge,params,'+');
                    disp(['New Tedge = ',num2str(Tedge)]);
                    %Run segmentation with new thresholds.
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                    
                case 'a' %add freehand regions to mask
                    
                    if length(currModifier) == 1 && strcmp(currModifier{:},'control') %closed region
                        fprintf('\n%s\n','Add closed freehand region to mask...')
                        modifiedRegions = modifyRegions(cbMask,...
                            'add',true,modifiedRegions,hf_segment);
                    elseif isempty(currModifier) %non-closed region
                        fprintf('\n%s\n','Add open freehand region to mask...')
                        modifiedRegions = modifyRegions(cbMask,...
                            'add',false,modifiedRegions,hf_segment);
                    end
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                    
                case 'r' %remove freehand regions in mask
                    
                    if length(currModifier) == 1 && strcmp(currModifier{:},'control') %closed region
                        fprintf('\n%s\n','Remove closed freehand region from mask...')
                        modifiedRegions = modifyRegions(cbMask,...
                            'remove',true,modifiedRegions,hf_segment);
                    elseif isempty(currModifier) %non-closed region
                        fprintf('\n%s\n','Remove open freehand region from mask...')
                        modifiedRegions = modifyRegions(cbMask,...
                            'remove',false,modifiedRegions,hf_segment);
                    end
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                    
                case 'escape'
    
                    fprintf('\n%s\n','INITIATING FINAL SEGMENTATION')
                    finalBool = true;
                    %Run segmentation with final thresholds
                    cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions);
                                        
                otherwise
                    fprintf('\n%s\n','Invalid key pressed!')
                    displaySegCmds();
                    continue;
            end
            
            %Do final cleanup of mask after manual editing:
            cbMaskTemp = imfill(logical(cbMask),'holes');
            cbMaskTemp = imreconstruct(maskNucShrink,cbMaskTemp);
            cbMask = cbMask .* bwareaopen(cbMaskTemp,params.cbSeg.minObjSize);

            %Apply modifications to regions of the mask.
            cbMask(modifiedRegions == 2) = 2; %added regions get kept (cbMask pixels = 2) NOTE THAT YOU SHOULD ONLY ADD PIXELS TO OBJECT YOU WANT TO KEEP (i.e. cbMask = 2)
            
            %Plot segmentation.
            %Get nuclei center coordinates for plotting:
            [tempR,tempC] = ind2sub(size(Icb),...
                 find(full(maskNucShrink)));

            plot_segmentation(Icb,Inuc,cbMask,tempR,tempC);

            if finalBool == true
                break;
            end
        end
    end
    
    close(hf_segment);
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displaySegCmds()

    fprintf('\n%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n',...
        'Enter a command...',...
        'uparrow - increase global intensity threshold',...
        'downarrow - decrease global intensity threshold',...
        'leftarrow - decrease Canny thresholds',...
        'rightarrow - increase Canny thresholds',...
        '"a" - add thin region',...
        'ctrl + "a" - add closed region',...
        '"r" - remove thin region',...
        'ctrl + "r" - remove closed region',...
        'ESC - finish segmentation');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keypress_segmentation(~,~)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modifiedRegions = modifyRegions(cbMask,modifyType,closedStatus,modifiedRegions,hf_segment)
        
     h_ModifyRegion = imfreehand('Closed',closedStatus);
     setColor(h_ModifyRegion,'b');

     if closedStatus == false
         try
             P0 = h_ModifyRegion.getPosition;
             D = [0; cumsum(sum(abs(diff(P0)),2))]; % Need the distance between points...
             P = interp1(D,P0,D(1):.5:D(end)); % ...to close the gaps
             P = unique(round(P),'rows');
             idx = sub2ind(size(cbMask),P(:,2),P(:,1));
         catch
             temp = createMask(h_ModifyRegion);
             idx = find(temp==1);
         end
     elseif closedStatus == true
         temp = createMask(h_ModifyRegion);
         idx = find(temp==1);
     end

     if strcmp(modifyType,'split')
         %Use temp when splitting to make connected components.
         temp = zeros(size(cbMask));
         temp(idx) = 1;
         temp = bwmorph(temp,'thin');
         temp = bwmorph(temp,'diag');
         idx = find(temp==1);
         modifiedRegions(idx) = 1;
     elseif strcmp(modifyType,'add')
         if closedStatus == false
             temp = zeros(size(cbMask));
             temp(idx) = 1;
             temp = bwmorph(temp,'thin');
             temp = bwmorph(temp,'diag');
             idx = find(temp==1);
         end
         modifiedRegions(idx) = 2;
     elseif strcmp(modifyType,'remove')
         if closedStatus == false
             temp = zeros(size(cbMask));
             temp(idx) = 1;
             temp = bwmorph(temp,'thin');
             temp = bwmorph(temp,'diag');
             idx = find(temp==1);
         end
         modifiedRegions(idx) = 3;  
     end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function TotsuNew = adjust_otsuThresh(TotsuPrev,params,plus_minus_cb)
    eval(['TotsuNew = TotsuPrev',plus_minus_cb,'params.cbSeg.otsuCorrectFrac;']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TedgeNew = adjust_edgeThresh(TedgePrev,params,plus_minus_edge)
    eval(['TedgeNew = TedgePrev ',plus_minus_edge,'params.cbSeg.edgeCorrectFrac;']);     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_segmentation(Icb,Inuc,cbMask,tempR,tempC)
        
        %Get CC of cb mask:
        CC = bwconncomp(logical(cbMask));
        boundPts = bwboundaries(logical(cbMask)); 
        
        %Plot cb image
        imshow(cat(3,zeros(size(Icb)),imadjust(Icb),zeros(size(Icb))),...
            [],'initialmagnification','fit','border','tight'); 
        hold on;
        %Plot touching object boundaries in red and non-touching in white.
        for i = 1:CC.NumObjects
            plot(boundPts{i}(:,2),boundPts{i}(:,1),'color','b','linewidth',2);
            hold on;
        end

        hold on;
        plot(tempC,tempR,'b.');
        
        hold off;   
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cbMask = run_cb_segmentation(Icb,maskNucShrink,Totsu,Tedge,params,modifiedRegions)

    %Don't allow thresholds to go below zero
    Totsu(Totsu < 0) = 0;
    if Tedge(1) < 0
        Tedge(1) = 0;
    end
    if Tedge(2) < 0
        Tedge(2) = eps;
    end

    %Binary image using Totsu threshold (Otsu corrected).
    cbMask = im2bw(Icb,Totsu);

    %Find edges using Tedge threshold (Canny corrected).
    Iedge = edge(Icb,'canny',Tedge);

    %Combine the Otsu foreground with the Canny edges that touch the the
    %foreground:
    cbMask = cbMask | imreconstruct(cbMask,Iedge);

    %Perform additional tweaks as needed.
    cbMask = imclose(cbMask,strel('disk',2));

    %Remove objects that are smaller than minObjSize # of pixels:
    cbMask = bwareaopen(cbMask,params.cbSeg.minObjSize);
    
    %Create a filled holes cell body channel used for identifying
    %single nontouching cells later:
    cbMaskFillHoles = imfill(cbMask,'holes');

    %Remove objects that are smaller than minObjSize # of pixels:
    cbMaskFillHoles = bwareaopen(cbMaskFillHoles,params.cbSeg.minObjSize);
    
    %Make changes based on modify removed regions by user.
    cbMaskFillHoles(modifiedRegions==3) = 0; %removed regions.

    %Perform reconstruction of cell body channel objects based on cell nuclei points
    %used as markers. Use actin channel with filled holes to ensure the
    %Euler number is due to nuclei present and not cells creating
    %"closed rings":
    cbMaskNucRec = imreconstruct(maskNucShrink,cbMaskFillHoles);

    %Use nuclei markers as holes in the cell body channel objects, and
    %compute euler number per cell body channel object. If euler number is
    %not zero (i.e. cells are touching), then remove them. This step
    %retains only those cells that do not touch each other.
    cbMaskRecWithNucHoles = cbMaskNucRec.*imcomplement(maskNucShrink);
    CCcb = bwconncomp(cbMaskRecWithNucHoles);

    %Extract Euler number of each cell body channel object with nuclei
    %holes:
    CCcbEulerNum = regionprops(CCcb,'EulerNumber'); %structure
    CCcbEulerNum = [CCcbEulerNum.EulerNumber]'; %convert to array

    %Loop through the cell body objects and set pixels of cell body to
    %2 if cell DOES NOT touch others, and leave to 1 if it touches
    %others.
    for idx = 1:numel(CCcbEulerNum)
        if CCcbEulerNum(idx) == 0 
            cbMaskRecWithNucHoles(CCcb.PixelIdxList{idx}) = 2;
        end
    end
    
    %Fill holes in objects and output final cell body
    %mask (logical):
    cbMask = imfill(cbMaskRecWithNucHoles,'holes');
    
    %Remove all spur pixels
    maskOld = cbMask;
    maskNew = bwmorph(cbMask,'spur');
    while 1
        if all(maskOld(:) == maskNew(:))
            break
        else
            maskOld = maskNew;
            maskNew = bwmorph(maskOld,'spur');
        end
    end
    cbMask = maskNew;
    
    %Set the objects in cbMask that touch the border to be 1's (recall 1's
    %in cbMask mean objects touch each other, and now also touch border,
    %while 2's mean those cells are good to go for tracking i.e. don't
    %touch each other).
    cbMask = imclearborder(cbMask) + (cbMask & ~imclearborder(cbMask));
    
    fprintf('\n%s\n','SEGMENTATION ITERATION COMPLETE')

end



