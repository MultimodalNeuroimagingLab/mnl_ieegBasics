%
% This script appends a column to the input electrodes table, with the header "seizure_zone". It overwrites the
% existing <seizure_zone> column if one already exists.
%
% Levels are:
%     "SOZ": "Seizure Onset Zone, the region where the recorded clinical seizures originated during the recording period.",
%     "IrritativeZone": "Region of cortex that generates interictal epileptiform discharges, but not seizures",
%     "EarlyPropagationZone": "Region of cortex that generates the initial seizure symptoms. Not seizure onset, but the propagation of seizure from SOZ into this region within first 3 seconds from seizure onset.",
%     "Resected": "Region of cortex that was resected",
%     "ResectedEdge": "Region of cortex that is within 1 cm of the edge of the resected area."
% Mutually exlusive (warning will appear for electrodes with mutually exclusive labels):
%       "SOZ" & "IrritativeZone"
%       "SOZ" & "EarlyPropagationZone"
%       "Resected" & "ResectedEdge"
%
% "n/a" is added by default to all unlabeled electrodes
%
% HH 2020

disp("Enter path to electrodes table:", 's');
[file,path] = uigetfile('*.tsv');
if isequal(file, 0)
    disp('User selected Cancel');
    return
else
    elecPath = fullfile(path,file);
    disp(['User selected ', elecPath]);
end

T = readtable(elecPath, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A', 'n/a'});

% check if seizure_zone has already been filled (default state is all-NaN column, which is numeric)
if ismember('seizure_zone', T.Properties.VariableNames) && ~isa(T.seizure_zone, 'double') 
    answer = questdlg('Input file is already labelled! You will overwrite the <seizure_zone> column... Do you want to continue?', ...
        'SOZ Labels', ...
        'Yes, overwrite', 'No, cancel', 'No, cancel');
    
    switch answer
        case 'No, cancel'
            return
    end  
end

goToEnd = 'n';
while strcmp(goToEnd, 'n')
    
    namesPadded = T.name; % print all electrode names
    namesPadded(end+1:ceil(length(namesPadded)/5)*5) = {''};
    fprintf("\nAll electrodes (n=%d):\n", length(T.name));
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
            sitesStr = input(sprintf('\nEnter all %s sites (separate entries with space/comma, use ''-'' for multiple (e.g. LA3-8), case-insensitive):\n', lev), 's'); % input sites for current level
            if isempty(sitesStr) % if nothing inputed
                ixes = [];
                break;
            end
            
            sites = strtrim(split(upper(sitesStr), {' ', ','}));
            sites(cellfun(@isempty, sites)) = []; % remove empty sites due to comma, space splitting
            sites = expandHyphens(sites); % expands hyphenated entries (e.g. 'LIT3-10' -> 'LIT3', 'LIT4', ... 'LIT10')
            
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
    if any(exc1)
        answer = questdlg(sprintf('Electrodes %s have mutually exclusive labels:\n"SOZ" & "IrritativeZone".\nKeep only:', string(join(T.name(exc1), ', '))), ...
        'SOZ vs. IrritativeZone', ...
        'SOZ', 'IrritativeZone', 'Keep both', 'SOZ');
    
        switch answer
            case 'SOZ'
                seizVals(exc1) = cellfun(@(x) x(x ~= 2), seizVals(exc1), 'UniformOutput', false); % remove IrritativeZone entries
                seizLabs(exc1) = strrep(seizLabs(exc1), 'SOZ, IrritativeZone', 'SOZ');
            case 'IrritativeZone'
                seizVals(exc1) = cellfun(@(x) x(x ~= 1), seizVals(exc1), 'UniformOutput', false); % remove SOZ entries
                seizLabs(exc1) = strrep(seizLabs(exc1), 'SOZ, IrritativeZone', 'IrritativeZone');
        end
    end
    
    exc2 = cellfun(@(x) all(ismember([1, 3], x)), seizVals); % indices that are labeled both "SOZ" and "EarlyPropagationZone (after SOZ vs irritative decision)"
    if any(exc2)
        answer = questdlg(sprintf('Electrodes %s have mutually exclusive labels:\n"SOZ" & "EarlyPropagationZone".\nKeep only:', string(join(T.name(exc2), ', '))), ...
        'SOZ vs. EarlyPropagationZone', ...
        'SOZ', 'EarlyPropagationZone', 'Keep both', 'SOZ');
    
        switch answer
            case 'SOZ'
                seizVals(exc2) = cellfun(@(x) x(x ~= 3), seizVals(exc2), 'UniformOutput', false); % remove EarlyPropagationZone entries
                seizLabs(exc2) = strrep(seizLabs(exc2), 'SOZ, EarlyPropagationZone', 'SOZ');
            case 'EarlyPropagationZone'
                seizVals(exc2) = cellfun(@(x) x(x ~= 1), seizVals(exc2), 'UniformOutput', false); % remove SOZ entries
                seizLabs(exc2) = strrep(seizLabs(exc2), 'SOZ, EarlyPropagationZone', 'EarlyPropagationZone');
        end
    end
    
    exc3 = cellfun(@(x) all(ismember([4, 5], x)), seizVals); % indices that are labeled both "Resected" and "ResectedEdge"
    if any(exc3)
        answer = questdlg(sprintf('Electrodes %s have mutually exclusive labels:\n"Resected" & "ResectedEdge".\nKeep only:', string(join(T.name(exc3), ', '))), ...
        'Resected vs. ResectedEdge', ...
        'Resected', 'ResectedEdge', 'Keep both', 'ResectedEdge');
    
        switch answer
            case 'Resected'
                seizVals(exc3) = cellfun(@(x) x(x ~= 5), seizVals(exc3), 'UniformOutput', false); % remove IrritativeZone entries
                seizLabs(exc3) = strrep(seizLabs(exc3), 'Resected, ResectedEdge', 'Resected');
            case 'ResectedEdge'
                seizVals(exc3) = cellfun(@(x) x(x ~= 4), seizVals(exc3), 'UniformOutput', false); % remove IrritativeZone entries
                seizLabs(exc3) = strrep(seizLabs(exc3), 'Resected, ResectedEdge', 'ResectedEdge');
        end
    end
    
    fprintf('\nSummary:\n')
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
T.seizure_zone = seizLabs; % will overwrite if column already exists
writetable(T, elecPath, 'FileType','text','Delimiter','\t');
fprintf('\nLabeled table successfully saved to:\n%s\n', elecPath);
