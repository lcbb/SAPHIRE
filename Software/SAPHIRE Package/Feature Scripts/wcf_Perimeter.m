%Feature: Perimeter
%Type: Scalar
%Source: MATLAB regionprops
%Description: the distance around the boundary of the region. regionprops
%computes the perimeter by calculating the distance between each adjoining
%pair of pixels around the border of the region. If the image contains
%discontiguous regions, regionprops returns unexpected results.
%Units: pixels (or um based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Perimeter(I,Imask,unitConv)

        featVal = regionprops(Imask,'Perimeter');
        featVal = featVal.Perimeter;
        %Convert to distance units:
        featVal = featVal * unitConv;
