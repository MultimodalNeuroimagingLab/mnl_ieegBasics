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
% Note: if the first element of a site pair starts with '0', the output will be padded in the tens position
%       E.g. input {'LP03-LP05'} or {'LP03-LP5'} returns {'LP03'; 'LP04'; 'LP05'}
%           but input {'LP3-LP05'} will just return {'LP3'; 'LP4'; 'LP5'}.
%       There is no issue when crossing to double-digit elements. E.g. input {'LP09-LP11'} returns {'LP09'; 'LP10'; 'LP11'} as expected
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
        strNumStart = elements{1}(isDigit(elements{1})); % before converting to double, used to check if starts with 0
        numStart = str2double(strNumStart); % first number of the site
        numEnd = str2double(elements{2}(isDigit(elements{2}))); % last number of the site
        
        if startsWith(strNumStart, '0') % pad single digits with '0' in tens place
            thisExpandedSite = arrayfun(@(x) sprintf('%s%02d', siteChars, x), numStart:numEnd, 'UniformOutput', false)';
        else
            thisExpandedSite = arrayfun(@(x) sprintf('%s%d', siteChars, x), numStart:numEnd, 'UniformOutput', false)';
        end
        expandedSites = [expandedSites; thisExpandedSite(:)];
    end
end