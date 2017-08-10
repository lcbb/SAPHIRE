%Feature: MajorAxisLength
%Type: Scalar
%Source: MATLAB regionprops
%Description: the length (in pixels) of the major axis of the ellipse that
%has the same normalized second central moments as the region. This
%property is supported only for 2-D input label matrices.
%Units: pixels (or um based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_MajorAxisLength(I,Imask,unitConv)

        featVal = regionprops(Imask,'MajorAxisLength');
        featVal = featVal.MajorAxisLength;
        %Convert to distance units:
        featVal = featVal * unitConv;
