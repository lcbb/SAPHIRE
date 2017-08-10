%Feature: Solidity
%Type: Scalar
%Source: MATLAB regionprops
%Description: the proportion of the pixels in the convex hull that are also
%in the region. Computed as Area/ConvexArea. This property is supported
%only for 2-D input label matrices.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Solidity(I,Imask,unitConv)

        featVal = regionprops(Imask,'Solidity');
        featVal = featVal.Solidity;
