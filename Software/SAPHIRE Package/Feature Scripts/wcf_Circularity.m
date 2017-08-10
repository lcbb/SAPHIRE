%Feature: Circularity (also called FormFactor)
%Type: Scalar
%Source: derived from MATLAB regionprops measurements
%Description: 4*pi*Area/Perimeter^2, equals 1 for a perfectly circular
%object.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_Circularity(I,Imask,unitConv)

        feats = regionprops(Imask,{'Area','Perimeter'});
        featVal = (4*pi*feats.Area)/feats.Perimeter^2;
  
