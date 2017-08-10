%Feature: FeretMax
%Type: Scalar
%Source: imFeretDiameters.m by David Legland
%Description: caliper distances of a binary object at given set of angles
%(theta) computed 0 to 179 orientations (degrees). This function outputs
%the maximum value of caliper length from the orientations.
%Units: pixels (or um based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_FeretMax(I,Imask,unitConv)

        thetas = 0:179;
        diams = imFeretDiameter(Imask,thetas);
        featVal = max(diams);
        
        %Convert to distance units
        featVal = featVal .* unitConv;
