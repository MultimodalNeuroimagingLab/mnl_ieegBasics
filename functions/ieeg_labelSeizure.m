%
% This script appends a column to the input electrodes table, with the header "seizure_label". It overwrites the
% existing <seizure_label> column if one already exists.
%
% Levels are: "SOZ", "IrritativeZone", "EarlyPropagationZone", "Resected", "ResectedEdge", "n/a"
% Mutually exlusive (warning will appear for electrodes with mutually exclusive labels):
%       "SOZ" & "IrritativeZone"
%       "SOZ" & "EarlyPropagationZone"
%       "Resected" & "ResectedEdge"
%
% "n/a" is added by default to all unlabeled electrodes
%
% HH 2020
%

while true % load electrodes table
    try
        elecPath = input("Path to electrodes table:\n", 's');
        T = readtable(elecPath, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    catch
        fprintf("ERROR: Invalid input path\n\n");
        continue
    end
    break
end

if ismember('seizure_label', T.Properties.VariableNames)
    warning('Input file is already labelled; you will be overwriting the <seizure_label> column.');
    input('Press any key to continue...\n', 's');
end

goToEnd = 'n';
while strcmp(goToEnd, 'n')
    
    namesPadded = T.name; % print all electrode names
    namesPadded(end+1:ceil(length(namesPadded)/5)*5) = {''};
    fprintf("All electrodes (n=%d):\n", length(T.name));
    disp(reshape(namesPadded, [], 5));

    seizLabs = repmat({''}, size(T, 1), 1); % string labels, ultimately appended to electrodes table
    seizVals = cell([size(T, 1), 1]); % numerical representation of labels

    for i = 1:5 % loop through the different types of levels
        switch i
            case 1, lev = 'SOZ';
            case 2, lev = 'IrritativeZone';
            case 3, lev = 'EarlyPropagationZone';
            case 4, lev = 'Resected';
            case 5, lev = 'ResectedEdge';
        end

        while true
            sitesStr = input(sprintf('\n%s sites (space/comma-separated, not case-sensitive):\n', lev), 's'); % input sites for current level
            if isempty(sitesStr) % if nothing inputed
                ixes = [];
                break;
            end
                
            sites = strtrim(split(upper(sitesStr), {' ', ','}));
            sites(cellfun(@isempty, sites)) = []; % remove empty sites due to comma, space splitting
            ixes = find(ismember(T.name, sites)); % indices in T.name where input sites are found
            lib = ismember(sites, T.name); % 0 for input sites that aren't found in T.name

            if all(lib)
                break
            else
                fprintf("ERROR: %s not found in electrodes table\n", string(join(sites(~lib), ', ')));
            end
        end

        for j = 1:length(ixes)
            if isempty(seizVals{ixes(j)}) % first label added for this site
                seizVals{ixes(j)} = i; % append numerical label
                seizLabs{ixes(j)} = lev; % append string label
            else
                seizVals{ixes(j)} = [seizVals{ixes(j)}, i];
                seizLabs{ixes(j)} = [seizLabs{ixes(j)}, ', ', lev];
            end
        end 
    end

    exc1 = cellfun(@(x) all(ismember([1, 2], x)), seizVals); % indices that are labeled both "SOZ" and "IrritativeZone"
    exc2 = cellfun(@(x) all(ismember([1, 3], x)), seizVals); % indices that are labeled both "SOZ" and "EarlyPropagationZone"
    exc3 = cellfun(@(x) all(ismember([4, 5], x)), seizVals); % indices that are labeled both "Resected" and "ResectedEdge"

    if any(exc1), warning('Electrodes %s have mutually exclusive labels: "SOZ" & "IrritativeZone".', string(join(T.name(exc1), ', '))); end
    if any(exc2), warning('Electrodes %s have mutually exclusive labels: "SOZ" & "EarlyPropagationZone".', string(join(T.name(exc2), ', '))); end
    if any(exc3), warning('Electrodes %s have mutually exclusive labels: "Resected" & "ResectedEdge".', string(join(T.name(exc3), ', '))); end

    for s = 1:length(seizLabs) % print labels for all 
        if ~isempty(seizLabs{s})
            fprintf('%s: %s\n', T.name{s}, seizLabs{s});
        end
    end
    fprintf('(Rest are ''n/a'')\n\n');

    while true
        goToEnd = lower(input('Confirm? [y/n], (''n'' will take you back to reinput levels):\n', 's'));
        if any(strcmp({'y', 'n'}, goToEnd)), break, end
    end
end

seizLabs(cellfun(@isempty, seizLabs)) = {'n/a'};
T.seizure_label = seizLabs; % will overwrite if column already exists
writetable(T, elecPath, 'FileType','text','Delimiter','\t');
fprintf('\nLabelled table successfully saved to:\n%s\n', elecPath);
