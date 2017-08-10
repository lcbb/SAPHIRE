%Feature: Eccentricity
%Type: Scalar
%Source: MATLAB regionprops
%Description: specifies the eccentricity of the ellipse that has the same
%second-moments as the region. The eccentricity is the ratio of the
%distance between the foci of the ellipse and its major axis length. The
%value is between 0 and 1. (0 and 1 are degenerate cases; an ellipse whose
%eccentricity is 0 is actually a circle, while an ellipse whose
%eccentricity is 1 is a line segment.) This property is supported only for
%2-D input label matrices.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Eccentricity(I,Imask,unitConv)

        featVal = regionprops(Imask,'Eccentricity');
        featVal = featVal.Eccentricity;
