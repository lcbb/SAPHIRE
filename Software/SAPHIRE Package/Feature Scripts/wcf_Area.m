%Feature: Area
%Type: Scalar
%Source: MATLAB regionprops
%Description: the actual number of pixels in the region. (This value may
%differ slightly from the value returned by bwarea, which weights different
%patterns of pixels differently.)
%Units: pixels^2 (or um^2 based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Area(I,Imask,unitConv)

        featVal = regionprops(Imask,'Area');
        featVal = featVal.Area;
        %Convert to distance units:
        featVal = featVal * unitConv^2;
