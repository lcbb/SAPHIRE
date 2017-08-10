%Feature: ConvexArea
%Type: Scalar
%Source: MATLAB regionprops
%Description:  specifies the number of pixels in 'ConvexImage'. This
%property is supported only for 2-D input label matrices.
%Units: pixels^2 (or um^2 based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_ConvexArea(I,Imask,unitConv)

        featVal = regionprops(Imask,'ConvexArea');
        featVal = featVal.ConvexArea;
        %Convert to distance units:
        featVal = featVal * unitConv^2;
