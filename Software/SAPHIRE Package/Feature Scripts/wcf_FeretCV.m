%Feature: FeretCV
%Type: Scalar
%Source: imFeretDiameters.m by David Legland
%Description: caliper distances of a binary object at given set of angles
%(theta) computed 0 to 179 orientations (degrees). This function outputs
%the CV value of caliper length from the orientations.
%Units: none
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_FeretCV(I,Imask,unitConv)

        thetas = 0:179;
        diams = imFeretDiameter(Imask,thetas);
        featVal = std(diams)/mean(diams);
        

