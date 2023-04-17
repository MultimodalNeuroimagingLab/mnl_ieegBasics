%% Preprocessing class for CCEP mef data
%   Mef data is loaded using the matmef package (Max van den Boom, 2020) and then subjected to modular preprocessing steps at the user's discretion.
%   Available preprocessing steps are highpass filtering, common average rereferencing, removing line noise, and baseline subtraction.
%   Data can be loaded all at once (ch x time points), have some preprocessing steps performed, and THEN converted to trial structure
%   (ch x time points x trials), or data can be loaded directly in trial structure via readMef3 (ch x time points x trials) to save memory.
%   Data can also be pruned down to only the desired channels.
%   Multiple mef objects, corresponding to multiple trials, can be concatenated together. This is particularly useful when different stim pairs have been split up
%   across multiple recording runs.
%   Incoming and outgoing CCEPs can be plotted for any/all channels and saved to a directory for preliminary viewing.
%
% DEPENDENCIES:
%   matmef, by Max van den Boom
%   mnl_ieegBasics
%
% USAGE:
%   
%   --CONSTRUCTOR--
%   Creates a new ccep_PreprocessMef object. Loads the channels table and (optionally) the events table. Hyphens are removed from channel names and
%   stim pair names (mef data will be loaded using original names, and hyphens will then be removed from metadata channel names).
%   >> mefObj = ccep_PreprocessMef(mefPath, channelsPath);
%   >> mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);
%       mefPath =               char. Path to .mefd folder to load
%       channelsPath =          char. Path to channels.tsv file, matching the mef data.
%       eventsPath =            (optional) char. Path to events.tsv file indicating where trials are in the mef data. If not given here, it must be given when
%                                   calling loadMefTrials
%      Returns:
%       mefObj =                ccep_PreprocessMef object. Contains get-fields: sub, channels, evts, metadata, srate, dataAll, data, tt, progress
%
%   --USEFUL FIELDS--
%   Manually set subject name. Subject name is automatically extrapolated from channelsPath in constructor, but will be randomly assigned if that fails. This is
%   used in file names when saving input/output plots via plotInputs/plotOutputs
%   >> mefObj.sub = 'sub-Name';
%   
%   Display preprocessing progress
%   >> mefObj.progress;
%      Returns: char array of preprocessing steps performed, user-readable
%  
%   Other readable properties are: channels, evts, metadata, srate, tt, and dataAll (after loadMefAll) or data (after loadMefTrials)
%
%   --OBJECT METHODS--
%   Detailed documentation is described below each method header, and accessed by typing "help ccep_PreprocessMef.methodName" or "help mefObj.methodName" if mefObj exists.
%   >> getpaths
%   >> filterEvents
%   >> loadMefAll
%   >> highpass
%   >> loadMefTrials
%   >> car
%   >> bipolar
%   >> pruneChannels
%   >> removeLN
%   >> subtractBaseline
%   >> plotInputs
%   >> plotOutputs
%
%   --COPYING AN INSTANCE--
%   Make a (shallow) copy of the mefObj. This means that all of the object properties (tables, arrays, etc.) are copied and independent of the original.
%   >> mefObj2 = copy(mefObj);
%       mefObj =                ccep_PreprocessMef object
%      Returns:
%       mefObj2 =               Shallow copy of mefObj
%
%   --CONCATENATING MULTIPLE INSTANCES--
%   See: help ccep_PreprocessMef.catMef
%
% EXAMPLES:
%
%   Ex. 1. Create a ccep_PreprocessMef object by loading directly as trials. Keep only 4.0mA and 6.0mA trials.
%   Apply common average referencing by 64-ch blocks, prune down to first 20 channels, remove line noise by spectrum interpolation, and subtract mean baseline.
%   Save all input and output CCEP plots to an output directory for preliminary viewing. Extract channels and data as a separate variables   
%   >> mefPath = '/path/to/file.mefd';                              % paths to inputs
%   >> channelsPath = '/path/to/channels.tsv';
%   >> eventsPath = '/path/to/events.tsv';
%   >> mefObj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath);   % construct ccep_PreprocessMef object, with eventsPath
%   >> mefObj.filterEvents('electrical_stimulation_current', {'4.0 mA', '6.0 mA'}); % keep only events with stim current == 4.0 or 6.0 mA
%   >> mefObj.loadMefTrials([-1, 2]);                               % load mef data for each trial from -1 to 2s around onset
%   >> mefObj.car(true);                                            % apply common average reference by 64-ch block
%   >> mefObj.pruneChannels(1:20);                                  % keep only the first 20 channels
%   >> mefObj.removeLN('SpectrumEstimation');                       % remove line noise using spectrum interpolation method
%   >> mefObj.subtractBaseline([-0.5, -0.05]);                      % subtract baseline from each trial calculated as the mean from -0.5s to -0.05s
%   >> mefObj.plotInputs([], [], [], '/path/to/output/folder');         % save plots of incoming CCEPs to all (20 pruned) channels to output folder, on [-1, 2]s time interval
%   >> mefObj.plotOutputs([], [], [], [], 'path/to/output/folder');         % save plots of outgoing CCEPs from all stim pairs to output folder, on [-1, 2]s time interval
%   >> channels = mefObj.channels;                                  % extract channels table and data as separate variables to use for later processing
%   >> sigdata = mefObj.data;
%
%   Ex. 2. Create a ccep_PreprocessMef object by loading all data (ch x timepoints). Apply high-pass filter, remove line noise using notch filter, then assemble into
%   trial structure. Apply common average rereferencing without restriction to 64-ch blocks and subtract median baseline. Display preprocessing progress, and open
%   plots for a few incoming CCEPs and a few outgoing CCEPs without automatically saving.
%   (Assume same paths as in ex. 1)
%   >> mefObj = ccep_PreprocessMef(mefPath, channelsPath);               % construct ccep_PreprocessMef object, this time not yet giving eventsPath
%   >> mefObj.loadMefAll;                                           % load all mef data in ch x time points structure
%   >> mefObj.highpass;                                             % apply high-pass filter to each channel
%   >> mefObj.removeLN('notch');                                    % remove line noise from each channel using notch filter
%   >> mefObj.loadMefTrials([-1.5, 2.5], eventsPath);               % convert dataAll to trial structure, from -1.5 to 2.5s around onset. eventsPath needed because it wasn't given at construction
%   >> mefObj.car;                                                  % apply common average reference without restricting to 64-ch blocks
%   >> mefObj.subtractBaseline([-1, -0.005], 'median');             % subtract baseline from each trial calculated as the median from -1s to -0.005s
%   >> mefObj.progress;                                             % display preprocessing progress
%   >> mefObj.plotInputs({'RA1', 'RA2', 'RC4'}, [-0.1, 1]);         % open plots of incoming CCEPs to 3 channels for inspection, from -0.1s to 1s around trial onset
%   >> mefObj.plotOutputs(5:10, [-0.1, 1]);                         % open plots of outgoing CCEPs from 6 stim sites for inspection, from -0.1s to 1s around trial onset
%
%   Ex. 3. Load 3 ccep_PreprocessMef objects and concatenate them
%   (Assume mefPaths, channelsPaths, and eventsPaths are length-3 cell array of char paths to their respective files)
%   >> mefObjs = cell(1, 3);                                        % initialize empty cell array to store mefObjs
%   >> for ii = 1:3
%   >>      mefObjs{ii} = ccep_PreprocessMef(mefPaths{ii}, channelsPaths{ii}, eventsPaths{ii}); % construct ccep_PreprocessMef object with each set of paths
%   >>      mefObjs{ii}.loadMefTrials([-1, 2]);                     % load mef data for each trial
%   >>      ...                                                     % perform the necessary preprocessing steps
%   >> end
%   >> mefObjCat = ccep_PreprocessMef.catMef(mefObjs{:});              % unpack cell array of ccep_PreprocessMef objects and concatenate to a single object
%
% Harvey Huang 2021
% 2023/04/17
%   - plotOutputs now only plots SEEG and ECOG channels
%   - implemented bipolar re-referencing method based on ieeg_bipolarSEEG.m from mnl_ieegBasics
%
%  	References discussing whether high pass filtering and baseline subtraction should be done:
%       Delorme, A. (2023). EEG is better left alone. Scientific Reports, 13(1), 2372.
%       van Driel, J., Olivers, C. N., & Fahrenfort, J. J. (2021). High-pass filtering artifacts in multivariate classification of neural time series data.
%           Journal of Neuroscience Methods, 352, 109080.
%
classdef ccep_PreprocessMef < matlab.mixin.Copyable % allow shallow copies
    
    properties % modifiable properties
        sub char % subject name
    end
    
    properties (SetAccess = private) % can read but not set directly
        channels table
        evts table
        metadata struct
        srate double
        dataAll double % "all" data, ch x time points
        data double % ch x time points x trials, either restructured from dataAll or loaded directly from mefPath
        tt double
        progress char
    end
    
    properties (Access = private) % cant read or set directly, used for class consistency purposes
        mefPath char
        channelsPath char
        eventsPath char
        chsPruned cell % stores info on channels that have been pruned
        concatenated logical = false; % true if this obj was a concatenated output from multiple objs
    end
    
    methods
        
        function obj = ccep_PreprocessMef(mefPath, channelsPath, eventsPath) % Specify paths and load tables
            obj.mefPath = mefPath;
            obj.channelsPath = channelsPath;
            obj.setSub;
            obj.channels = readtableRmHyphens(channelsPath);
            obj.progress = 'Constructed object';
            
            if nargin < 3 % no events. Load later when putting into trial structure
                disp('No events file loaded');
                return
            end
            
            obj.eventsPath = eventsPath;
            try
                obj.evts = readtableRmHyphens(eventsPath, 'electrical_stimulation_site', 1);
            catch
                obj.evts = readtable(eventsPath, 'FileType', 'text', 'Delimiter', '\t'); % if not CCEP events
            end
            obj.progress = sprintf('%s\n> Loaded events', obj.progress);
            
            try
                stimType = join(unique(obj.evts.electrical_stimulation_type), ', ');
                dur = unique(obj.evts.duration); assert(length(dur) == 1, 'More than 1 stim duration in events');
                current = join(unique(obj.evts.electrical_stimulation_current), ', ');
                fprintf('Electrical stimulation type: %s; duration: %.02e; current: %s\n', stimType{1}, dur, current{1});
            catch
                %pass
            end
            
            try % remove bad trials
                fprintf('Removing %d not good events\n', sum(~strcmpi(obj.evts.status, 'good')));
                obj.evts(~strcmpi(obj.evts.status, 'good'), :) = [];
                obj.progress = sprintf('%s\n> Removed non-''good'' events', obj.progress);
            catch
                warning('Could not remove non-''good'' events from obj.evts.status');
            end
        end
        
        function paths = getPaths(obj) % returns paths used to load data
        %   Returns the file paths that were used to load mef data, channels, and events
        %   >> paths = mefObj.getPaths;
        %      Returns:
        %           paths =                 struct. Contains fields mef, channels, events. Each field stores the char array path used to load that data. If 
        %                                   mefObj is a concatenated ccep_PreprocessMef object (catMef, below), then each field will contain multiple lines,
        %                                   corresponding to each concatenated ccep_PreprocessMef object, in order of concatenation.
            paths = struct();
            paths.mef = obj.mefPath;
            paths.channels = obj.channelsPath;
            paths.events = obj.eventsPath;
        end
        
        function filterEvents(obj, varargin)
        %   Filters events by column values. Keeps only the events where events.<Name> corresponds to elements in events.<Value>. Deletes the other events
        %   Can input multiple Name, Value pairs in a row (Name, Value, Name, Value, ...)
        %   >> mefObj.filterEvents(Name, Value, ...);
        %       Name =                  char array. Column name in mefObj.evts to indicate which column to look for matches. E.g. 'electrical_stimulation_current'
        %       Value =                 char array or cell array of chars to match in mefObj.evts.<Name>. E.g. '6.0 mA' or {'4.0 mA', '6.0 mA'}
            for ii = 1:2:length(varargin)-1 % iterate thru and remove all lines that do not satisfy input requirements
                key = varargin{ii}; val = varargin{ii+1};
                if ischar(val), val = {val}; end
                valStr = join(val, ', ');
                if ~sum(ismember(obj.evts.(key), val))
                    error('Cannot find any matches for %s in obj.%s', valStr{1}, key);
                end
                obj.evts(~ismember(obj.evts.(key), val), :) = [];
                obj.progress = sprintf('%s\n> Filtered obj.%s for %s', obj.progress, key, valStr{1});
            end
        end
        
        function loadMefAll(obj) % load ch x time points matrix of entire mef data
        %   Load all time series data from mefPath in a ch x timepoints matrix. Data is stored as mefObj.dataAll. Also configures mefObj.metadata and mefObj.srate fields
        %   >> mefObj.loadMefAll;
            disp('Loading FULL mef data');
            channelsMef = readtable(obj.channelsPath, 'FileType', 'text', 'Delimiter', '\t'); % preserve names to load mef
            [obj.metadata, obj.dataAll] = readMef3(obj.mefPath, [], channelsMef.name);
            obj.srate = obj.metadata.time_series_metadata.section_2.sampling_frequency;
            %obj.dataAll = obj.applyConversionFactor(obj.dataAll); % apply conversion upon loading. 2021/10/05 - commented out because this is being done in matmef
            obj.changeName(); % remove hyphenated names
            obj.progress = sprintf('%s\n> Loaded all ch x time points mef data', obj.progress);
        end
        
        function highpass(obj)
        %   Applies high-pass filtering to each channel in mefObj.dataAll (after mefObj.loadMefAll). Uses ieeg_highpass.m from mnl_ieegBasics
        %   Butterworth filter: Stop band (30 dB attenuation) < 0.05 Hz, Pass band (-3dB cutoff freq) > 0.50 Hz. Forward-reverse filtering using filtfilt
        %   >> mefObj.highpass;
            assert(~isempty(obj.dataAll), 'highpass can only be applied to dataAll (channels x timepoints)');
            disp('Applying highpass filter');
            obj.dataAll = ieeg_highpass(obj.dataAll', obj.srate)';
            obj.progress = sprintf('%s\n> High-pass filter', obj.progress);
        end
        
        function loadMefTrials(obj, trange, eventsPath)
        %   Loads or converts mef data into trial structure (ch x time points x trials), as mefObj.data
        %   If dataAll is already loaded with loadMefAll, then the trial structure is pulled directly from dataAll (fast). mefObj.dataAll IS THEN DELETED.
        %   If dataAll hasn't been loaded, then data is loaded fresh from mefPath using readMef3 (slow).
        %   Moving forward, only mefObj.data OR mefObj.dataAll will exist, not both!
        %   Configures mefObj.tt, which is a 1xn double time point vector
        %   >> mefObj.loadMefTrials(trange);
        %   >> mefObj.loadMefTrials(trange, eventsPath);
        %       trange =                1x2 double. [start, stop] time, in seconds, around each event onset to load as trials
        %       eventsPath =            (optional) char. Only needed if not given during object construction. Path to events.tsv file. If mefObj.events already exists
        %                                   and this input is given, this will overwrite the existing mefObj.events
            if nargin >= 3
                if ~isempty(obj.evts), warning('Overwriting existing events'); end
                obj.eventsPath = eventsPath;
                try
                    obj.evts = readtableRmHyphens(eventsPath, 'electrical_stimulation_site', 1);
                catch
                    obj.evts = readtable(eventsPath, 'FileType', 'text', 'Delimiter', '\t'); % if not CCEP events
                end
                obj.progress = sprintf('%s\n> Loaded events', obj.progress);
                
                try % print some info electrical stim info about the events
                    stimType = join(unique(obj.evts.electrical_stimulation_type), ', ');
                    dur = unique(obj.evts.duration); assert(length(dur) == 1, 'More than 1 stim duration in events');
                    current = join(unique(obj.evts.electrical_stimulation_current), ', ');
                    fprintf('Electrical stimulation type: %s; duration: %.02e; current: %s\n', stimType{1}, dur, current{1});
                catch
                    pass
                end
            elseif isempty(obj.evts)
                error('Path to events file must be input as 3rd argument');
            end
            
            if isempty(obj.metadata) % did not get from loading dataAll
                obj.metadata = readMef3(obj.mefPath); % get sampling frequency
                obj.srate = obj.metadata.time_series_metadata.section_2.sampling_frequency;
            end
            
            obj.tt = (0:(trange(end)-trange(1))*obj.srate-1)/obj.srate + trange(1); % samples around stim to return
            ranges = [obj.evts.sample_start + obj.tt(1)*obj.srate, ...
                      obj.evts.sample_start + obj.tt(1)*obj.srate + length(obj.tt)]; % [start end+1] samples for each trial. start is inclusive and corresponds to time 0. end is exclusive
            
            if ~isempty(obj.dataAll)
                disp('Converting dataAll into trial structure');
                obj.data = nan([height(obj.channels), length(obj.tt), height(obj.evts)]); % chs x time x trials
                for ii = 1:height(obj.evts)
                    obj.data(:, :, ii) = obj.dataAll(:, ranges(ii, 1)+1 : ranges(ii, 2)); % Conversion from 0-index argument to Matlab 1-index. Equivalent to readMef3 on ranges(tr, 1):ranges(tr, 2).
                end
                obj.dataAll = []; % clear memory
                obj.progress = sprintf('%s\n> Converted data to trials', obj.progress);
            else
                disp('Loading data trials from mefd file');
                channelsMef = readtable(obj.channelsPath, 'FileType', 'text', 'Delimiter', '\t'); % preserve names to load mef
                [~, obj.data] = readMef3(obj.mefPath, [], channelsMef.name, 'samples', ranges); % readmef3 is 0 index! So this loads up to ranges(end)-1
                obj.changeName(); % remove hyphenated names
                %obj.data = obj.applyConversionFactor(obj.data); % apply conversion upon loading. 2021/10/05 - commented out because it is being done in matmef
                obj.progress = sprintf('%s\n> Loaded data in trials', obj.progress);
            end
        end
        
        function car(obj, by64)
        %   Common average rereferencing. If mefObj.dataAll exists, this call applies ieeg_car from mnl_ieegBasics, across all channels with 'status'=='good'.
        %   If mefObj.data exists, this call applies ccep_CAR or ccep_CAR64blocks from mnl_ieegBasics on each trial separately.
        %   This CANNOT be performed on a concatenated ccep_PreprocessMef object (ccep_PreprocessMef.catMef below), and SHOULD NOT be performed after mefObj.pruneChannels is called.
        %   >> mefObj.car;
        %   >> mefObj.car(by64);
        %       by64 =                  (optional) logical, default = false. If true, ccep_CAR64blocks is used instead of ccep_CAR, which applies a separate rereference
        %                                   to each contiguous block of 64 SEEG channels. Both functions threshold channels that form the reference by early variance and 
        %                                   late variance criteria. This input is ignored if trial structure data (mefObj.data) doesn't exist (since ieeg_car is used)
            if nargin < 2, by64 = false; end
            assert(~obj.concatenated, 'Common average rereferencing cannot be performed after concatenating multiple ccep_PreprocessMef objects');
            if ~isempty(obj.chsPruned), warning('Applying CAR after pruning channels may be inaccurate!!'); end
            
            if ~isempty(obj.dataAll) && isempty(obj.data)
                warning('Applying ieeg_car to dataAll because trial data hasn''t been created. by64 input ignored');
                chans2incl = find(strcmpi(obj.channels.type, 'SEEG') & strcmpi(obj.channels.status, 'good'));
                obj.dataAll = ieeg_car(obj.dataAll', chans2incl)';
                obj.progress = sprintf('%s\n> Common average rereference', obj.progress);
            elseif isempty(obj.dataAll) && ~isempty(obj.data)
                if by64
                    assert(isempty(obj.chsPruned), 'Cannot apply 64-SEEG-ch block CAR after pruning channels');
                    disp('Applying CAR by 64-SEEG-ch block to trial data');
                    obj.data = ccep_CAR64blocks(obj.data, obj.tt, obj.channels, obj.evts.electrical_stimulation_site);
                    obj.progress = sprintf('%s\n> Common average rereference by 64-SEEG-ch block', obj.progress);
                else
                    disp('Applying CAR to trial data');
                    obj.data = ccep_CAR(obj.data, obj.tt, obj.channels, obj.evts.electrical_stimulation_site);
                    obj.progress = sprintf('%s\n> Common average rereference', obj.progress);
                end
            else, error('Either dataAll or data (not both) needs to exist');
            end
        end
        
        function bipolar(obj, seg5, seg6)
        % Bipolar re-referencing using ieeg_bipolarSEEG.m from mnl_ieegbasics. Operates on mefObj.dataAll or mefObj.data, whichever exists.
        % Only rereferences the channels with type == 'SEEG'. All non-SEEG channels are concatenated below the re-referenced SEEG channels
        % Channels table is modified to match the bipolar re-referenced channels. THIS MEANS THE NAMES, ORDER, AND NUMBER OF CHANNELS CHANGE.
        % Bipolar channel status = 'bad' if EITHER input channel in the bipolar pair was bad
        % >> mefObj.bipolar
        % >> mefObj.bipolar(seg5, seg6);
        %       seg5 =              (optional) char or cell array of char. List of leads that are segmented in groups of 5. E.g., 'LA'.
        %                               Bipolar pairs that cross segments will not be included in the output channels
        %       seg6 =              (optional) char or cell array of char. List of leads that are segmented in groups of 6, like seg5 above.
        %
            if nargin < 3 || isempty(seg6), seg6 = {}; end % no 6-segmented leads
            if nargin < 2 || isempty(seg5), seg5 = {}; end % no 5-segmented leads
            if ~isempty(obj.chsPruned), warning('Channels were pruned. Bipolar re-referencing will fail if there are non-consecutive positions on a lead'); end
            
            % apply only to SEEG channels
            chsSEEG = obj.channels(strcmp(obj.channels.type, 'SEEG'), :);
            badChans = find(~strcmpi(chsSEEG.status, 'good'));
            
            % bipolar rereference across all data
            if ~isempty(obj.dataAll) && isempty(obj.data)
                fprintf('Applying Bipolar re-referencing to dataAll\n');
                
                % apply only to SEEG data
                dataIn = obj.dataAll(strcmp(obj.channels.type, 'SEEG'), :);
                [dataOut, bipolarChans, badChans] = ieeg_bipolarSEEG(dataIn', chsSEEG.name, badChans, seg5, seg6);
                
                % vertically concatenate bipolar-rereferenced SEEG data with the non-SEEG data
                obj.dataAll = [dataOut'; obj.dataAll(~strcmp(obj.channels.type, 'SEEG'), :)];
            
            % bipolar rereference for each trial separately
            elseif isempty(obj.dataAll) && ~isempty(obj.data)
                fprintf('Applying Bipolar re-referencing to data\n');
                
                dataBip = [];
                for ii = 1:height(obj.evts) % iterate across trials
                    
                    % Apply only to SEEG data
                    dataIn = obj.data(strcmp(obj.channels.type, 'SEEG'), :, ii);
                    if ii == 1, verbose = true; else, verbose = false; end
                    [dataOut, bipolarChans, badChans] = ieeg_bipolarSEEG(dataIn', chsSEEG.name, badChans, seg5, seg6, verbose);
                        
                    dataBip = cat(3, dataBip, dataOut'); % concatenate trials
                end
                
                % concatenate with non-SEEG data (in 1st dimension)
                obj.data = cat(1, dataBip, obj.data(~strcmp(obj.channels.type, 'SEEG'), :, :));
                
            end
            
            % modify the channels table to accommodate bipolar-rereferenced data
            eles = split(bipolarChans, '-', 2); % get the minuend of each bipolar pair, to find matching rows in channels
            minuends = eles(:, 1);
            chsSEEG = chsSEEG(ismember(chsSEEG.name, minuends), :); % get the subset of chsSEEG that match th eminuends
            chsSEEG.name = bipolarChans; chsSEEG.status(badChans) = {'bad'}; % cases where the bipolar ch is bad because subtrahend was bad
            obj.channels = [chsSEEG; obj.channels(~strcmp(obj.channels.type, 'SEEG'), :)]; % concatenate with non-SEEG channels
            
            obj.progress = sprintf('%s\n> Bipolar rereference', obj.progress);
            
        end
        
        function pruneChannels(obj, chList) % keep only relevant channels
        %   Prune channels to keep only channels of interest (to save memory). Modifies mefObj.channels, mefObj.dataAll (if exists), mefObj.data (if exists)
        %   >> mefObj.pruneChannels(chList);
        %       chList =                cell array of char, or integer. If cell array, each element is a channel name to keep. E.g. {'RA1', 'RA2', 'RC4'}.
        %                                   If integer array, each element is a channel index to keep. E.g. 21:30 keeps the 21st - 30th channels
            if isnumeric(chList), chList = obj.channels.name(chList); end % convert indices to names
            if ischar(chList), chList = {chList}; end % if single ch given
            
            assert(sum(ismember(obj.channels.name, chList)) > 0, 'No channels will be returned');
            obj.chsPruned = setdiff(obj.channels.name, chList); % keep track which channels have been removed
            
            if ~isempty(obj.dataAll) && isempty(obj.data)
                obj.dataAll(~ismember(obj.channels.name, chList), :) = [];
            elseif isempty(obj.dataAll) && ~isempty(obj.data)
                obj.data(~ismember(obj.channels.name, chList), :, :) = [];
            end
            obj.channels(~ismember(obj.channels.name, chList), :) = [];
            obj.progress = sprintf('%s\n> Pruned channels', obj.progress);
        end
        
        function removeLN(obj, method, opts)
        %   Line noise removal, using either a notch filter (mnl_ieegBasics ieeg_notch.m) or the spectrum interpolation method implemented by Mewett, Nazeran, and Reynolds.
        %   >> mefObj.removeLN('notch');
        %   >> mefObj.removeLN('notch', f);
        %   >> mefObj.removeLN('SpectrumEstimation');
        %   >> mefObj.removeLN('SpectrumEstimation', opts);
        %       'notch' =               (Faster) Removes line noise with notch filter. Specify f to be either 60 (default) or 50 to remove line noise and harmonics at
        %                                   [60, 120, 180] Hz or [50, 100, 150] Hz.
        %       'SpectrumEstimation'    (Slower) Removes line noise with the spectrum interpolation method (mnl_ieegBasics/external/removeLineNoise_SpectrumEstimation.m).
        %                                   This is advantageous when there are high-amplitude transient artifacts (e.g. from stimulation) in the data.
        %                                   opts gets passed to the opts input of that function. Default = 'LF = 60, NH = 3, HW = 3' (line noise fundamental freq =
        %                                   60 Hz, remove 3 harmonics, half-width = 3 Hz). E.g. also set window to 4096 samps: 'LF = 60, NH = 3, HW = 3, M = 4096'.
            switch lower(method)
                case 'notch' % probably code a more sophisticated notch filter later
                    if nargin < 3, opts = 60; end % base notch frequency
                    if ischar(opts), opts = str2double(opts); end % ensure numerical
                    if ~isempty(obj.dataAll) && isempty(obj.data)
                        disp('Removing line noise on dataAll with ieeg_notch');
                        obj.dataAll = ieeg_notch(obj.dataAll', obj.srate, opts, true)';
                    elseif isempty(obj.dataAll) && ~isempty(obj.data)
                        disp('Removing line noise on trial data with ieeg_notch');
                        for ii = 1:size(obj.data, 1)
                            fprintf('.');
                            obj.data(ii, :, :) = ieeg_notch(squeeze(obj.data(ii, :, :)), obj.srate, opts, 2, false);
                        end
                        fprintf('\n');
                    else, error('Either dataAll or data (not both) needs to exist');
                    end
                    obj.progress = sprintf('%s\n> Removed line noise at %d Hz with notch filter', obj.progress, opts);
                    
                case 'spectrumestimation'
                    if nargin < 3, opts = 'LF = 60, NH = 3, HW = 3'; end % default params
                    if ~isempty(obj.dataAll) && isempty(obj.data)
                        assert(isempty(obj.data), 'trial data shouldn''t exist if dataAll still exists');
                        disp('Removing line noise on dataAll with removeLineNoise_SpectrumEstimation');
                        obj.dataAll = removeLineNoise_SpectrumEstimation(obj.dataAll, obj.srate, opts, true);
                        if ~isreal(obj.dataAll)
                            warning('Line noise removal yielded complex outputs for some channels (e.g. digitalinputs). Keeping real part only');
                            obj.dataAll = real(obj.dataAll);
                        end
                    elseif isempty(obj.dataAll) && ~isempty(obj.data)
                        disp('Removing line noise on trial data with removeLineNoise_SpectrumEstimation');
                        for ii = 1:size(obj.data, 1)
                            fprintf('.');
                            try
                                obj.data(ii, :, :) = removeLineNoise_SpectrumEstimation(squeeze(obj.data(ii, :, :))', obj.srate, opts, false)';
                            catch
                                warning('Could not remove line noise for channel %s (likely contains nan values)', obj.channels.name{ii});
                            end
                        end
                        fprintf('\n');
                        if ~isreal(obj.data)
                            warning('Line noise removal yielded complex outputs for some channels (e.g. digitalinputs). Keeping real part only');
                            obj.data = real(obj.data);
                        end
                    else, error('Either dataAll or data (not both) needs to exist');
                    end
                    obj.progress = sprintf('%s\n> Removed line noise with spectrum estimation; settings: %s', obj.progress, opts);
                    
                otherwise
                    error("Unexpected method arg. Options = 'Notch' or 'SpectrumEstimation'");
            end
            
        end
        
        function subtractBaseline(obj, baseRange, method)
        %   Subtracts mean or median baseline voltage from each trial in mefObj.data, calculated on an input time range
        %   >> mefObj.subtractBaseline(baseRange);
        %   >> mefObj.subtractBaseline(baseRange, method);
        %       baseRange =             1x2 double. [start, stop] time, in seconds, to calculate the baseline voltage on for each trial
        %       method =                (optional) char, default = 'mean'. 'mean' or 'median', method to calculate baseline voltage
            assert(~isempty(obj.data), 'Baseline can only be subtracted from trial data');
            assert(baseRange(1)>=obj.tt(1) && baseRange(end)<=obj.tt(end), 'baseline time range not within tt time range');
            if nargin < 3, method = 'mean'; end
            switch lower(method)
                case 'mean'
                    obj.data = obj.data - mean(obj.data(:, obj.tt>=baseRange(1) & obj.tt<=baseRange(end), :), 2);
                case 'median'
                    obj.data = obj.data - median(obj.data(:, obj.tt>=baseRange(1) & obj.tt<=baseRange(end), :), 2);
                otherwise
                    error("method must be given as 'mean' or 'median'");
            end
            obj.progress = sprintf('%s\n> Subtracted %s baseline on [%.3f, %.3f] s', obj.progress, method, baseRange(1), baseRange(end));
        end
        
        function plotInputs(obj, chs, trange, yspace, dir, fmt)
        %   Generate plots of incoming CCEPs to all channels or channels of interest. Can only be run after mefObj.loadMefTrials creates the trial structure data
        %   Plots are opened if no <dir> input given, or saved without opening if <dir> input is given.
        %   If <dir> is given, a metadata file named 'metadata_<mefObj.sub>_incomingCCEP.txt' is also written to <dir>, containing info on the date/time the command was
        %   executed, the paths used to to load the data, preprocessing steps executed, and channels saved.
        %   >> mefObj.plotInputs;
        %   >> mefObj.plotInputs(chs, trange, yspace);
        %   >> mefObj.plotInputs(__, dir);
        %   >> mefObj.plotInputs(__, dir, fmt);
        %       chs =                   (optional) cell array or integer array. If not given, incoming CCEPs to all channels are plotted. If given, can be cell array of channel
        %                                   names, e.g. {'RA1', 'RA2', 'RC4'}; or integer array of channel indices, e.g. 21:30, to plot incoming CCEPs for.
        %       trange =                (optional) 1x2 double. [start, stop] time, in seconds, around each trial to plot the incoming CCEPs. If not given, the entire
        %                                   mefObj.tt time interval is plotted.
        %       yspace =                (optional) positive double. The spacing between adjacent stim site lines. Smaller value corresponds to more gain in the plotted signals. Default = 500
        %       dir =                   (optional) char. Path to directory (folder) to save plots to. If not given, the plot at each channel will be opened without saving.
        %                                   If given, the plots will not be open but rather directly saved in <dir> instead.
        %       fmt =                   (optional) char. String formatting used to name the plots if saved to <dir>. Must contain 2 '%s', where the first one corresponds
        %                                   to the subject name (mefObj.sub) and the second one corresponds to the channel saved. Default = '%s_incomingCCEP_%s'.
        %                                   E.g. using default fmt, if mefObj.sub == 'sub-harbo' and the current channel plotted is 'RA1', the file saved would be at
        %                                   '<dir>/sub-harbo_incomingCCEP_RA1.png'.
            assert(~isempty(obj.data), 'Plots can only be generated from trial data');
            
            if nargin < 6, fmt = '%s_incomingCCEP_%s'; end % string formatting to save images
            assert(count(fmt, '%s') == 2, 'fmt needs to contain exactly 2 ''%s''s');
            if nargin < 5, dir = []; end
            if nargin < 4 || isempty(yspace), yspace = 500; end
            if nargin < 3 || isempty(trange), trange = [min(obj.tt), max(obj.tt)]; end
            if nargin < 2 || isempty(chs), chs = obj.channels.name; end % plot all channels
            if isnumeric(chs), chs = obj.channels.name(chs); end
            if ischar(chs), chs = {chs}; end % if single input
            if length(chs) > 10 && isempty(dir)
                answer = questdlg(sprintf('Are you sure you want to open %d figures', length(chs)), ...
                                          'Plot inputs', 'Yes, I have SO MUCH RAM to spare', 'No, cancel', 'No, cancel');
                if strcmp(answer, 'No, cancel'), return; end
            end  
            
            stimSites = unique(obj.evts.electrical_stimulation_site, 'stable');
            
            for ii = 1:length(chs)
                dataCh = squeeze(obj.data(strcmpi(obj.channels.name, chs{ii}), :, :));
                assert(~isempty(dataCh), 'No input data to channel %s', chs{ii});
                
                if ~isempty(dir) % save to dir, don't show plot
                    f = figure('Position', [200, 200, 600, 800], 'visible', 'off');
                else
                    figure('Position', [200, 200, 600, 800]);
                end
                hold on;
                for jj = 1:size(stimSites, 1)
                    yline(-(jj-1)*yspace, 'Color', 0.5*[1, 1, 1]);
                    if ismember(chs{ii}, split(stimSites{jj}, '-')), continue; end % don't plot stimulated channel
                    
                    meanTrial = mean(dataCh(:, strcmp(obj.evts.electrical_stimulation_site, stimSites{jj})), 2);
                    plot(obj.tt, meanTrial-(jj-1)*yspace, 'LineWidth', 1);
                end
                xline(0, 'Color', 'r');
                plot([0.05 0.05], [yspace*0.5 yspace*1.5], 'k', 'LineWidth', 2); % scale
                text(0.055, yspace, sprintf('%d \\muV', yspace));
                hold off
                
                xlim([min(trange), max(trange)]);
                ylim([-(jj+1)*yspace, 2*yspace]); % [-1000 from bottom to 1000 from top]
                xlabel('Time (s)');
                ylabel('Stimulated Electrodes');
                set(gca, 'YTick', (-yspace)*(size(stimSites, 1)-1:-1:0), 'YTickLabel', flip(stimSites), 'FontSize', 10);
                title(sprintf('Input to: %s', chs{ii}));
                
                if ~isempty(dir) % save to dir
                    fprintf('Saving INPUT plot for %s\n', chs{ii});
                    saveas(f, fullfile(dir, sprintf(fmt, obj.sub, chs{ii})), 'png');
                    close(f);
                end
            end
            
            if ~isempty(dir), obj.writeMeta(fullfile(dir, sprintf('metadata_%s.txt', sprintf(fmt, obj.sub, 'ch'))), chs); end % metadata about parameters
        end
        
        function plotOutputs(obj, sites, trange, yspace, maxCh, dir, fmt)
        %   Generate plots of outgoing CCEPs from all stimulated electrode pairs or stim pairs of interest.
        %   Plots are opened if no <dir> input given, or saved without opening if <dir> input given.
        %   If <dir> is given, a metadata file named 'metadata_<mefObj.sub>_outgoingCCEP.txt' is also written to <dir>, containing info on
        %   the date/time the command was executed, the paths used to to load the data, preprocessing steps executed, and stim pairs saved.
        %   >> mefObj.plotOutputs;
        %   >> mefObj.plotOutputs(sites, trange, yspace, maxCh);
        %   >> mefObj.plotOutputs(__, dir);
        %   >> mefObj.plotOutputs(__, dir, fmt);
        %       sites =                 (optional) cell array or integer array. If not given, outgoing CCEPs from all stim electrode pairs are plotted. If given, can be
        %                                   cell array of stim pair names, e.g. {'RA1-RA2', 'RA2-RA3', 'RC4-RC4'}; or integer array of stim pair indices, in the order
        %                                   that they are uniquely encountered going down the events file, e.g. 11:20 for the 11th - 20th unique stim pairs, to plot
        %                                   outgoing CCEPs for.
        %       trange =                (optional) 1x2 double. [start, stop] time, in seconds, around each trial to plot the outgoing CCEPs. If not given, the entire
        %                                   mefObj.tt time interval is plotted.
        %       yspace =                (optional) positive double. The spacing between adjacent channel lines. Smaller value corresponds to more gain in the plotted signals. Default = 500
        %       maxCh =                 (optional) positive integer. The maximum number of channels to plot on one plot. Channels will be divided ~evenly on multiple plots. Default = 128.
        %       dir =                   (optional) char. Path to directory (folder) to save plots to. If not given, the plot for each stim pair will be opened without saving.
        %                                   If given, the plots will not be open but rather directly saved in <dir> instead.
        %       fmt =                   (optional) char. String formatting used to name the plots if saved to <dir>. Must contain 2 '%s', where the first one corresponds
        %                                   to the subject name (mefObj.sub) and the second one corresponds to the stim pair saved. Default = '%s_outgoingCCEP_%s'.
        %                                   E.g. using default fmt, if mefObj.sub == 'sub-harbo' and the current stim pair plotted is 'RA1-RA2', the file saved would be at
        %                                   '<dir>/sub-harbo_outgoingCCEP_RA1-RA2.png'.
            assert(~isempty(obj.data), 'Plots can only be generated from trial data');
            
            if nargin < 7, fmt = '%s_outgoingCCEP_%s'; end % string formatting to save images
            assert(count(fmt, '%s') == 2, 'fmt needs to contain exactly 2 ''%s''s');
            if nargin < 6, dir = []; end
            if nargin < 5 || isempty(maxCh), maxCh = 128; end
            if nargin < 4 || isempty(yspace), yspace = 500; end
            if nargin < 3 || isempty(trange), trange = [min(obj.tt), max(obj.tt)]; end
            if nargin < 2 || isempty(sites), sites = unique(obj.evts.electrical_stimulation_site, 'stable'); end % plot all stim sites
            if isnumeric(sites)
                allSites = unique(obj.evts.electrical_stimulation_site, 'stable'); % sort in order of first encounter
                sites = allSites(sites);
            end % indices
            if ischar(sites), sites = {sites}; end % if single input
            if length(sites) > 10 && isempty(dir)
                answer = questdlg(sprintf('Are you sure you want to open %d figures', length(sites)), ...
                                          'Plot outputs', 'Yes, I have SO MUCH RAM to spare', 'No, cancel', 'No, cancel');
                if strcmp(answer, 'No, cancel'), return; end
            end
            
            chs = obj.channels.name; % all channels
            status = obj.channels.status; % status of all channels
            dataChs = obj.data; % pull out to remove non-ephys channels
            
            % Only plot the SEEG and ECOG channels
            try
                chs(~ismember(upper(obj.channels.type), {'SEEG', 'ECOG'})) = [];
                status(~ismember(upper(obj.channels.type), {'SEEG', 'ECOG'})) = [];
                dataChs(~ismember(upper(obj.channels.type), {'SEEG', 'ECOG'}), :, :) = [];
            catch
                warning('Could not isolate only SEEG and ECoG channels to plot');
            end
            
            for ii = 1:length(sites)
                dataStim = dataChs(:, :, strcmpi(obj.evts.electrical_stimulation_site, sites{ii})); % transpose to match plotInput structure
                assert(~isempty(dataStim), 'No output data from stim site %s', sites{ii});
                dataStim = mean(dataStim, 3)'; % mean across stim trials
                
                % how many figures to divide channels into
                nBlocks = ceil(length(chs)/maxCh);
                chBlockSize = ceil(length(chs)/nBlocks);
                
                for nn = 1:nBlocks
                    
                    chStart = (nn-1)*chBlockSize + 1;
                    chEnd = min(nn*chBlockSize, length(chs));
                    chsCurr = chs(chStart:chEnd);
                    statusCurr = status(chStart:chEnd);
                
                    if ~isempty(dir) % save to dir, don't show plot
                        f = figure('Position', [200, 200, 600, 800], 'visible', 'off');
                    else
                        figure('Position', [200, 200, 600, 800]);
                    end
                    hold on;

                    for jj = 1:length(chsCurr)
                        yline(-(jj-1)*yspace, 'Color', 0.5*[1, 1, 1]);
                        if ismember(chsCurr{jj}, split(sites{ii}, '-')), continue; end % skip stimulated channel
                        if ~strcmpi(statusCurr{jj}, 'good'), continue; end % 
                        
                        plot(obj.tt, dataStim(:, chStart+jj-1)-(jj-1)*yspace, 'LineWidth', 1); % channel
                    end
                    xline(0, 'Color', 'r');
                    plot([0.05 0.05], [yspace*0.5 yspace*1.5], 'k', 'LineWidth', 2); % scale
                    text(0.055, yspace, sprintf('%d \\muV', yspace));
                    hold off

                    xlim([min(trange), max(trange)]);
                    ylim([-(jj+1)*yspace, 2*yspace]); % [-1000 from bottom to 1000 from top]
                    xlabel('Time (s)');
                    ylabel('Recording Electrodes');
                    set(gca, 'YTick', (-yspace)*(size(chsCurr, 1)-1:-1:0), 'YTickLabel', flip(chsCurr), 'FontSize', 10);
                    
                    if nBlocks == 1
                        title(sprintf('Output from: %s', sites{ii}));
                        fmtCurr = fmt;
                    else % indicate in title and file name which block we are in
                        title(sprintf('Output from: %s (%d/%d)', sites{ii}, nn, nBlocks));
                        fmtCurr = sprintf('%s_%dof%d', fmt, nn, nBlocks);
                    end

                    if ~isempty(dir) % save to dir
                        fprintf('Saving OUTPUT plot (%d) for %s\n', nn, sites{ii});
                        saveas(f, fullfile(dir, sprintf(fmtCurr, obj.sub, sites{ii})), 'png');
                        close(f);
                    end
                
                end
                
            end
            
            if ~isempty(dir), obj.writeMeta(fullfile(dir, sprintf('metadata_%s.txt', sprintf(fmt, obj.sub, 'stimSite'))), sites); end % metadata about parameters
        end
        
    end
    
    methods(Static)
        
        function objCat = catMef(varargin)
        %   STATIC function.
        %   Concatenates across trials from multiple ccep_PreprocessMef objects into a single (new) ccep_PreprocessMef object. Input ccep_PreprocessMef objects are unaffected.
        %   Each ccep_PreprocessMef object must be in trial structure (possess mefObj.data, and must have the exact same channel names, srate, and time points (mefObj.tt)
        %   to be concatenated. Additionally, the ccep_PreprocessMef objects SHOULD have the same subject name (mefObj.sub) preprocess steps executed (mefObj.progress);
        %   warnings will be given if these 2 matches are not met.
        %
        %   >> mefObjCat = ccep_PreprocessMef.catMef(mefObj1, mefObj2, ... mefObjN);
        %       mefObj1 ... mefObjN =   ccep_PreprocessMef objects whose trials are to be concatenated, in the order that they are given as input.
        %      Returns:
        %       mefObjCat =             ccep_PreprocessMef object. The fields: sub, metadata, progress are copied from mefObj1.
        %                                   mefObjCat.evts are vertically concatenated across all mefObjs. mefObjCat.data are concatenated across dim 3 (trials) across all mefObjs.
        %                                   mefObjCat.mefPath, channelsPath, eventsPath are concatenated across all mefObjs
            assert(isa(varargin{1}, 'ccep_PreprocessMef'), 'Inputs must be ccep_PreprocessMef objects');
            objCat = copy(varargin{1}); % copy object so as to modify first obj
            assert(~isempty(objCat.data), 'Object 1 not in trial structure -- cannot be concatenated');
            
            if length(varargin) == 1
                fprintf('Only one mefObj input given. Nothing is concatenated.\n');
                return
            end
            
            for ii = 2:length(varargin)
                obj2 = varargin{ii}; % next object in line to be concatenated
                
                % Check for necessary/recommended matches
                assert(~isempty(obj2.data), 'Object 1 not in trial structure -- cannot be concatenated');
                assert(objCat.srate == obj2.srate, 'Sampling frequency mismatch with obj%d', ii);
                assert(all(strcmp(objCat.channels.name, obj2.channels.name)), 'Channels mismatch with obj%d', ii);
                assert(all(objCat.tt == obj2.tt), 'Time vector mismatch with obj%d', ii);
                if ~strcmp(objCat.sub, obj2.sub), warning('Subject name mismatch with obj%d', ii); end
                if ~strcmp(objCat.progress, obj2.progress), warning('Preprocessing mismatch with obj%d. Output progress will not capture this', ii); end
                
                % Concatenate events and data
                objCat.evts = [objCat.evts; obj2.evts];
                objCat.data = cat(3, objCat.data, obj2.data);
                objCat.mefPath = sprintf('%s\n%s', objCat.mefPath, obj2.mefPath);
                objCat.channelsPath = sprintf('%s\n%s', objCat.channelsPath, obj2.channelsPath);
                objCat.eventsPath = sprintf('%s\n%s', objCat.eventsPath, obj2.eventsPath);
            end
            fprintf('Merged across %d ccep_PreprocessMef objects\n', ii);
            objCat.concatenated = true;
            objCat.progress = sprintf('%s\n> Merged across %d ccep_PreprocessMef objects', objCat.progress, ii);
        end
        
    end
    
    
    methods (Access = private) % helper functions
        
        function setSub(obj) % set subject name upon construction
            try
                pathInfo = dir(obj.channelsPath);
                eles = split(pathInfo.name, '_');
                obj.sub = eles{1};
                fprintf('Subject name: %s\n', obj.sub);
            catch
                obj.sub = sprintf('random%s', char(randi([65, 90], 1, 4)));
                warning('Could not extract subject id from mefPath. Setting to ''%s\n''', obj.sub);
            end
        end
        
        function sig = applyConversionFactor(obj, sig) % apply conversion factor to loaded mef data
            convFac = obj.metadata.time_series_metadata.section_2.units_conversion_factor;
            if convFac < 0 % TODO: change this after correcting xltec natus conversion factor
                warning('Ignoring negative sign in conversion factor');
                convFac = -convFac;
            end
            sig = sig*convFac;
        end
        
        function changeName(obj) % replace hyphenated names in metadata with hyphen-removed names
            for ii = 1:length(obj.channels.name)
                obj.metadata.time_series_channels(ii).name = obj.channels.name(ii);
            end
        end
        
        function writeMeta(obj, path, chs) % write a metadata file for incoming/outgoing CCEP plots
            f = fopen(path, 'w');
            fprintf(f, '#%s\n', datetime);
            fprintf(f, 'MefPath: %s\nChannelsPath: %s\nEventsPath: %s\n', obj.mefPath, obj.channelsPath, obj.eventsPath);
            fprintf(f, '#\nProgress:\n%s\n', obj.progress);
            fprintf(f, '#\nSaved:');
            for ii = 1:size(chs, 1)
                fprintf(f, sprintf(' %s', chs{ii}));
            end
            fclose(f);
        end
        
    end
    
end
