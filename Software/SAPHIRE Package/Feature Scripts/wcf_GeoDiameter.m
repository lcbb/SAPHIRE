%Feature: geodesic diameter
%Type: Scalar
%Source: David Legland 
%Description: length of the longest geodesic path within a particle. A
%geodesic path is the shortest path (series of neighbor pixels within
%particle) that connects two points in the particle. The geodesic diameter
%is a good alternative to other shape descriptors such as perimeter.
%Units: pixels (or um based on unitConv)
%
% Author: Simon Gordonov
% Date: 09/30/14

function featVal = wcf_GeoDiameter(I,Imask,unitConv)

        featVal = imGeodesicDiameter(Imask);
        %Convert to distance units:
        featVal = featVal * unitConv;
