%
%   Reorders an electrode table or file according to the row order of
%   bipolar channel names from a table or file
%
%   It ignores upper or lower case in the channels or electrode names.
%   
%   elecsOut = ieeg_sortElectrodes(electrodes, channels, saveFile)
%
%       electrodes      = the electrodes to sort. Either as a table or as a
%       filepath to electrodes.tsv file channels        = the channels that
%       the electrodes will be sorted to (by name). Either as a table or as
%       a filepath to electrodes.tsv file
%
%   Returns:
%       elecsOut        = Table with the rows of electrodes.tsv sorted in
%       order of channel names in channels.tsv
%                         Channel names not found in electrodes.tsv (i.e.
%                         ends of leads without coordinates) will be NaN
%                         rows in elecsOut, in order to align overall table
%                         structure with channels NaN rows are saved to
%                         file as 'n/a' to match BIDS formatting
%
%
%   Harvey Huang, Multimodal Neuroimaging Lab, Mayo Clinic, 2021
%   MM Updating name for saving the sorted electrodes file, January 2023
%   DH Changed code to deal with bipolar channels, August 2026

function [elecsOut1, elecsOut2, elecsMean] = ieeg_sortBipElectrodes(elecs, bip_channels)

    if ~istable(elecs)
        elecs = readtable(elecs, 'FileType', 'text', 'Delimiter', '\t'); % keep hyphens for filename
    end
    
    if ~istable(bip_channels)
        channs = readtable(bip_channels, 'FileType', 'text', 'Delimiter', '\t');
    else
        channs = bip_channels;
    end
    
    % preallocate elecsSave1 (first electrode in pair)
    varTypes = cellfun(@class, table2cell(elecs(1, :)), 'UniformOutput', false);
    varTypes(strcmp('char', varTypes)) = {'string'}; % so that Matlab doesn't scream at me for preallocating with 'char'
    varTypes{1} = 'name'; % ignore first col for indexing purposes
    
    elecsSave1 = table('Size', [height(bip_channels), length(elecs.Properties.VariableNames)], ...
                            'VariableNames', elecs.Properties.VariableNames, ...
                            'VariableTypes', repmat({'string'}, [1, length(elecs.Properties.VariableNames)]));
    elecsSave1.name = channs.name; % save with original channel names

    % preallocate elecsSave1 (second electrode in pair)
    elecsSave2 = elecsSave1; 

    % preallocate average electrode position for plotting
    mean_xyz = zeros(height(elecsSave1),3);

    % now assign values to bipolar electrodes table
    for kk = 1:height(elecsSave1)
        el1 = extractBefore(elecsSave1.name{kk},'-');
        el2 = extractAfter(elecsSave1.name{kk},'-');
        if ~isempty(el1) && ~isempty(el2) && ~isempty(find(ismember(upper(elecs.name),upper(el1)),1)) && ~isempty(find(ismember(upper(elecs.name),upper(el2)),1))
            el1_ind = find(ismember(upper(elecs.name),upper(el1)));
            el2_ind = find(ismember(upper(elecs.name),upper(el2)));

            mean_xyz(kk,:) = mean([elecs.x([el1_ind el2_ind]) elecs.y([el1_ind el2_ind]) elecs.z([el1_ind el2_ind])]);
            
            elecsSave1(kk,:) = elecs(el1_ind,:);
            elecsSave2(kk,:) = elecs(el2_ind,:);
            elecsSave1.name{kk} = el1;
            elecsSave2.name{kk} = el2;
        elseif ~isempty(el1) && ~isempty(el2) 
            elecsSave1.name{kk} = el1;
            elecsSave2.name{kk} = el2;
            mean_xyz(kk,:) = NaN;
        else
            mean_xyz(kk,:) = NaN;
        end
    end
    
    % make sure to keep numbers as doubles
    elecsSave1 = convertvars(elecsSave1, "x", "double");
    elecsSave1 = convertvars(elecsSave1, "y", "double");
    elecsSave1 = convertvars(elecsSave1, "z", "double");
    elecsSave1 = convertvars(elecsSave1, "Destrieux_label", "double");
    elecsSave2 = convertvars(elecsSave2, "x", "double");
    elecsSave2 = convertvars(elecsSave2, "y", "double");
    elecsSave2 = convertvars(elecsSave2, "z", "double");
    elecsSave2 = convertvars(elecsSave2, "Destrieux_label", "double");
    
    % replace missing stuff for cells or strings with 'n/a'
    for kk = 1:length(elecsSave1.Properties.VariableNames)
        if iscell(elecsSave1.(elecsSave1.Properties.VariableNames{kk})) || isstring(elecsSave1.(elecsSave1.Properties.VariableNames{kk}))
            thisColumn = elecsSave1.(elecsSave1.Properties.VariableNames{kk});
            thisColumn(ismissing(thisColumn)) = {'n/a'};
            elecsSave1.(elecsSave1.Properties.VariableNames{kk}) = thisColumn;
        end
    end
    % replace missing stuff for cells or strings with 'n/a'
    for kk = 1:length(elecsSave2.Properties.VariableNames)
        if iscell(elecsSave2.(elecsSave2.Properties.VariableNames{kk})) || isstring(elecsSave2.(elecsSave2.Properties.VariableNames{kk}))
            thisColumn = elecsSave2.(elecsSave2.Properties.VariableNames{kk});
            thisColumn(ismissing(thisColumn)) = {'n/a'};
            elecsSave2.(elecsSave2.Properties.VariableNames{kk}) = thisColumn;
        end
    end

    elecsOut1 = elecsSave1; 
    elecsOut2 = elecsSave2; 

    name = channs.name;
    x = mean_xyz(:,1);
    y = mean_xyz(:,2);
    z = mean_xyz(:,3);
    hemisphere = elecsSave1.hemisphere;
    hemisphere(ismissing(hemisphere)) = {'n/a'};
    elecsMean = table(name,x,y,z,hemisphere);


end
 