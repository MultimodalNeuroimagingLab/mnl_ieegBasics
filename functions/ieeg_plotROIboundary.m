
function ieeg_plotROIboundary(g,roi_labels,varargin)

% function to draw boundaries of selected ROIs
% 
% input:
%     g: gifti (pial or inflated) for slected hemisphere (l/r)
%     roi_labels: integer labels corresponding to ROIs labels. vx1 matrix v-->#of vertices
%     varargin{1}: boundary color, [r g b], default is white
%     
% output:
%     renders the ROI boundaries on the loaded brain figure
%     
%
% Example usage:
%   g = gifti( fullfile( freesurferPath, ['inflated.' hemi '.surf.gii']));
%   sulcal_labels = read_curv( fullfile( freesurferPath, 'surf', [hemi 'h.sulc']));
%   roi_labels = ieeg_loadAtlas(freesurferPath, 'r', 'Wang')
%   ieeg_RenderGifti( g, sulcal_labels);
%   plotROIboundary( g, roi_labels, [0 0.99 0])
% 
% ZQ 2025
%

    if nargin == 2
        bClr = [0.99 0.99 0.99];
    elseif nargin == 3
        bClr = varargin{1};
    else
        assert(1==0, 'Error: Check the number of parameters')
    end

    faces = g.faces;
    vertices = g.vertices;
    
    % Find edges between vertices of different ROIs (boundary edges)
    edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])]; % All edges in surface
    edges = sort(edges, 2); % Sort vertex indices for consistency
    
    % Identify edges where the two vertices belong to different ROIs
    roi_edges = roi_labels(edges(:,1)) ~= roi_labels(edges(:,2));
    
    % Boundary edges are those connecting vertices with different ROI labels
    boundary_edges = edges(roi_edges, :);
    
    % Plot boundary edges in bClr
    for i = 1:size(boundary_edges,1)
        v1 = boundary_edges(i,1);
        v2 = boundary_edges(i,2);
        plot3(vertices([v1 v2],1), vertices([v1 v2],2), vertices([v1 v2],3), ...
            'Color', bClr, 'LineWidth', 2);
    end


end
