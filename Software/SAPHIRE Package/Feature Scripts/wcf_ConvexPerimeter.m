%Feature: ConvexPerimeter
%Type: Scalar
%Source: derived from MATLAB regionprops measurements
%Description: perimeter of the convex hull of the object.
%Units: pixels(or nm based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_ConvexPerimeter(I,Imask,unitConv)

        feat = regionprops(Imask,'ConvexImage');
        featVal = regionprops(feat.ConvexImage,'Perimeter');
        %Convert to distance units:
        featVal = featVal.Perimeter * unitConv;

