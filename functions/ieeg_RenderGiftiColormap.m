function tH = ieeg_RenderGiftiColormap(g,c)
% function to render a gifti 
% 
% input:
%     g: gifti file with faces and vertices
%     c: color for every vertex
%
% output:
%   th: returns trimesh handle so you can change it
%       for example, tH.FaceAlpha = 0.5 will make the rendering transparent 
%
% Viewing Angle: can be changed with ecog_ViewLight(90,0), changes both
% angle and light accordingly
%
% Example usage:
%   gii = gifti('/path/to/gifti');
%   [~, label, colortable] = read_annotation('path/to/.aparc.a2009s.annot'); % Destrieux labels
% 
%   label = changem(label, 1:size(colortable.table(:,5)), colortable.table(:,5)); % change labels to range from 1:n
%   label_cmap = colortable.table(:, 1:3)./256; % normalize between 0 and 1
%   label_cmap(1, :) = 0.7*[1 1 1]; % change unlabeled from black to gray
%   label_cmap = brighten(label_cmap, 0.7); % make more pastel-colored
% 
%   figure; ieeg_RenderGiftiColormap(gii, label, label_cmap, colortable.struct_names);
%
%
% DH 2017

tH = trimesh(g.faces, g.vertices(:,1), g.vertices(:,2), g.vertices(:,3), c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
l1 = light;
lighting gouraud
material([.3 .9 .2 50 1]); 
axis off
set(gcf,'Renderer', 'zbuffer')
view(270, 0);
set(l1,'Position',[-1 0 1])


