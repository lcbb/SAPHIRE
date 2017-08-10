%Feature: Waviness (also called Roughness)
%Type: Scalar
%Source: derived from MATLAB regionprops measurements
%Description: ratio of the convex hull perimeter to the perimeter of the
%object. A common shape factor.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Waviness(I,Imask,unitConv)

        feats = regionprops(Imask,{'ConvexImage','Perimeter'});
        Pcvx = regionprops(feats.ConvexImage,'Perimeter');
        featVal = Pcvx.Perimeter / feats.Perimeter;
  
