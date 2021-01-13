%
% This function takes an nx1 char cell array of electrode sites and expands hyphen-contracted items into their non-contracted
% counterparts.
% Helper function to ieeg_labelSeizure.m
%
% expandedSites = expandHyphens(sites)
%       sites =         nx1 cell array of electrode sites, where hyphens mean all sites in between the 2 numbers
%
% Returns:
%       expandedSites = mx1 cell array of electrode sites, where all hyphenated entries in <sites> have been replaced by
%       expanded entries
%
%   E.g. sites =         {'LIT2'}
%                        {'LP2-4'}
%                        {'LZ3'}
%
%        expandedSites = {'LIT2'}
%                        {'LP2'}
%                        {'LP3'}
%                        {'LP4'}
%                        {'LZ3'}
%
% HH 2020
%

function expandedSites = expandHyphens(sites)
    expandedSites = {};
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    
    for i=1:length(sites)
        if ~contains(sites{i}, '-')
            expandedSites = [expandedSites; sites{i}]; % do not modify site if it has no hyphen
            continue 
        end
        
        elements = strtrim(split(sites{i}, '-'));
        assert(numel(elements) == 2, sprintf('%s does not contain exactly 1 hyphen', sites{i}));
        
        siteChars = elements{1}(~isDigit(elements{1})); % the characters of the site (e.g. 'LA')
        numStart = str2double(elements{1}(isDigit(elements{1}))); % the first number of the site
        numEnd = str2double(elements{2}(isDigit(elements{2}))); % the last number of the site
        
        thisExpandedSite = arrayfun(@(x) sprintf('%s%d', siteChars, x), numStart:numEnd, 'UniformOutput', false)';
        expandedSites = [expandedSites; thisExpandedSite(:)];
    end
end