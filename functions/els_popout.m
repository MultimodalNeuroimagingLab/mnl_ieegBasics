function xyzs_pop = els_popout(xyzs, th, phi, varargin)
%
% Pops out xyz coordinates slightly in the viewing angle [th phi]
%
% INPUT:
%   xyzs    = nx3 matrix where each row is a different electrode. Columns are [xcoord ycoord zcoord]
%   th      = theta (azimuth) of viewing angle, in degrees
%   phi     = phi (elevation) of viewing angle, in degrees
%   mag     = (optional) magnitude of pop-out relative to max xyzs x-dim.
%             default = 0.05
%
% OUTPUT:
%   xyzs_pop = nx3 matrix of same dimension as input xyzs, but the coordinates are adjusted for pop-out
%
assert(nargin<=4, 'Exceeded number of inputs');

if nargin==4
    mag = varargin{1};
else
    mag = 0.05;
end

a_offset = mag*max(abs(xyzs(:,1)))*[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];
xyzs_pop = xyzs+repmat(a_offset, size(xyzs,1), 1);      