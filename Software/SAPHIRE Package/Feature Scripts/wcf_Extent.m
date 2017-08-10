%Feature: Extent
%Type: Scalar
%Source: MATLAB regionprops
%Description: specifies the ratio of pixels in the region to pixels in the
%total bounding box. Computed as the Area divided by the area of the
%bounding box. This property is supported only for 2-D input label
%matrices.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Extent(I,Imask,unitConv)

        featVal = regionprops(Imask,'Extent');
        featVal = featVal.Extent;
