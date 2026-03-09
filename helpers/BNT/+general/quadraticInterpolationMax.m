% Calculate maximum of given function using quadratic interpolation
%
% This function first interpolates the given function around it's maximum value,
% then the maximum is calculated as vertex of quadratic function. This function
% is intended to be used for obtaining a maximum of a binned function. After binning,
% the maximum value will be given with bin resolution. By using the interpolation,
% it is possible to increase precision.
% Units of vx/vy are determined by units of bins/y.
%
%  USAGE
%   [vx, vy] = quadraticInterpolationMax(bins, y, limits)
%   bins        Bin/x values for function y.
%   y           Function for which maximum value is obtained.
%   limits      [l1 l2] two element vector that defines area of interpolation.
%               The function y is interpolated on interval [max(y)-l1 max(y)+l2].
%   vx          X-value of maximum of y. It can be NaN if quadratic function has
%               it's minimum at that point instead of maximum, or if vx goes beyond
%               the range of bins.
%   vy          Y-value of maximum of y. It can be NaN if vx is NaN.
%
function [vx, vy] = quadraticInterpolationMax(bins, y, limits)
    if isempty(y)
        vx = nan;
        vy = nan;
        return;
    end
    limits = abs(limits);
    [~, i] = max(y);
    self = max(1, i-limits(1)):min(length(bins), i+limits(2));
    p = polyfit(bins(self), y(self), 2);
    vx = -p(2)/(2*p(1)); % x-value of quadratic vertex, standard equation

    % assign NaN to data with maxima beyond the bin range or with positive
    % second derivative (exhibiting a minimum instead of maximum)
    if p(1) > 0 || vx > bins(end) || vx < bins(1)
        vx = nan;
    end
    vy = polyval(p, vx);
end
