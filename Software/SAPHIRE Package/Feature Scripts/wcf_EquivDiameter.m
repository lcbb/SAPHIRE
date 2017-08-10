%Feature: EquivDiameter
%Type: Scalar
%Source: MATLAB regionprops
%Description: specifies the diameter of a circle with the same area as the
%region. Computed as sqrt(4*Area/pi). This property is supported only for
%2-D input label matrices.
%Units: pixels (or um based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_EquivDiameter(I,Imask,unitConv)

        featVal = regionprops(Imask,'EquivDiameter');
        featVal = featVal.EquivDiameter;
        %Convert to distance units:
        featVal = featVal * unitConv;
