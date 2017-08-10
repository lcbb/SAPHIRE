%Feature: AxisRatio
%Type: Scalar
%Source: derived from MATLAB regionprops measurements
%Description: major axis length/minor axis length (axis lengths derived
%from best filling ellipse with same second moments as the object).
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_AxisRatio(I,Imask,unitConv)

        feats = regionprops(Imask,{'MajorAxisLength','MinorAxisLength'});
        featVal = feats.MajorAxisLength / feats.MinorAxisLength;
  
