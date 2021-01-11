function [config, path_table] = NM_process(config, step, use_adjustments)
%--------------------------------------------------------------------------
% NM_process NuMorph processing pipeline designed to perform channel 
% alignment, intensity adjustment, and stitching on multi-channel 
% light-sheet images.
%--------------------------------------------------------------------------
% Usage:  
% [config, path_table] = NM_process(config, step, use_adjustments)
% 
% Inputs:
%   config - config structure for processing.
%   step - 'stitch', 'align', 'intensity'. Perform only one of the 3 steps.
%   use_adjustments - (default: true) apply intensity adjustments if
%   performing 1 step.
%
% Output:
%   config - config structure after processing
%   path_table - table of filenames used and additional image information
%--------------------------------------------------------------------------

% Load configuration from .mat file, if not provided
if nargin<1
    config = load(fullfile('templates','NM_variables.mat'));
elseif isstring(config)
    config = load(fullfile(config,'NM_variables.mat'));
elseif ~isstruct(config)
    error("Invalid configuration input")
end

% Check config structure to make sure it's correct
if ~isfield(config,'adjust_intensity') || ~isfield(config,'align_channels')
    error("config structure is not from processing template.")
end

% Default to run full pipeline
if nargin<2
    step = 'process';
end

% Apply intensity adjustments if only aligning or stitching
if nargin<3 
    if ~isequal(config.adjust_intensity,"false")
        use_adjustments = true;
    else
        use_adjustments = false;
    end
end

% Check step input
if ~ismember({step},{'process','align','stitch','intensity'})
    error("Invalid processing step specified.")
end

% Make an output directory
if exist(config.output_directory,'dir') ~= 7
    mkdir(config.output_directory);
end

% Make a variables directory
config.var_directory = fullfile(config.output_directory,'variables');
if exist(config.var_directory,'dir') ~= 7
    mkdir(config.var_directory);
end

fprintf("%s\t Working on sample %s \n",datetime('now'),config.sample_name)

%% Create directories
% Update image directory if using processed images
if ~isequal(config.use_processed_images,"false")
    config.img_directory = fullfile(config.output_directory,config.use_processed_images);
    if ~exist(config.img_directory,'dir')
        error("Could not locate processed image directory %s\n",config.img_directory)
    end
    if isequal(config.use_processed_images,"aligned")
        warning("Images already aligned. Skipping channel alignment")
        config.channel_alignment = "false";
    end
end

%% Read image filename information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate table containing image information
path_table = path_to_table(config);

% Count number of x,y tiles for each channel
nchannels = length(unique(path_table.channel_num));
ncols = zeros(1,nchannels); nrows = zeros(1,nchannels);
for i = 1:nchannels
    path_sub = path_table(path_table.channel_num == i,:);
    ncols(i) = length(unique(path_sub.x));
    nrows(i) = length(unique(path_sub.y));
end

% Check if all resolutions are equal
equal_res = all(cellfun(@(s) config.resolution{1}(3) == s(1,3),config.resolution));

%% Intensity Adjustment
% Configure intensity adjustments
if use_adjustments || isequal(step,'intensity')
    % If single step, will only be loading intensity adjustments
    if isequal(step,'stitch') || isequal(step,'align')
        config.adjust_intensity = "true";
    end
    [config, path_table] = perform_intensity_adjustment(config, path_table, nrows, ncols);
end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note no specific checks on tiles/markers included
if nargin>1 && isequal(step,'stitch')    
    [config, path_table] = perform_stitching(config, path_table);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'align')
    [config, path_table] = perform_channel_alignment(config, path_table, equal_res);
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
elseif nargin>1 && isequal(step,'intensity')
    if nargout == 1; path_table = []; end
    if nargout < 1; clear config; end
    return
end

%% Run multiple steps
% If equal numbers of column and row tiles, perform alignment on each tile,
% then do stitching on aligned images. This assumes the same resolution for
% each channel. If the channels are imaged at different resolutions, stitch
% each channel if multi-tile. Then run alignment
if all(ncols == max(ncols)) && all(nrows == max(nrows))
    % Equal number of tiles
    % Channel Alignment
    [config, path_table] = perform_channel_alignment(config,path_table,equal_res);
    
    % Stitching
    if ncols(1)*nrows(1) == 1
        fprintf("%s\t 1 tile detected for marker %s \n",datetime('now'),config.markers(1));
    else
        [config, path_table] = perform_stitching(config, path_table);
    end
else
    % Different number of tiles
    fprintf("%s\t Different number of tiles between channels \n",datetime('now'))

    % Check resolutions are present
    assert(length(config.resolution) == length(config.markers),"Specify resolution "+...
        "for each channel as cell array");
    
    % Stitching
    for i = 1:length(config.markers)
        % Continue if only 1 tile
        if ncols(i)*nrows(i) == 1
            fprintf("%s\t 1 tile detected for marker %s \n",datetime('now'),config.markers(i));
            continue
        else
            config.stitch_sub_channel = i;
        end
        [config, path_table] = perform_stitching(config, path_table);
    end
    
    % Channel Alignment
    % First make sure there's only 1 tile in each channel
    assert(max(path_table.x) == 1,...
        "More than 1 tile detected for marker %s. Images "+...
        "must be stitched first when multiple resolutions are present.",...
        config.markers(ncols == max(ncols)))
    assert(max(path_table.y) == 1,...
        "More than 1 tile detected for marker %s. Images "+...
        "must be stitched first when multiple resolutions are present.",...
        config.markers(nrows == max(nrows)))
    
    [config, path_table] = perform_channel_alignment(config, path_table, equal_res);
end

%% Return outputs
if nargout == 1; path_table = []; end
if nargout < 1; config = []; end
    
end


function [config, path_table] = perform_intensity_adjustment(config, path_table, nrows, ncols)
% Intensity adjustments

if length(config.adjust_tile_position) == 1
    config.adjust_tile_position = repmat(config.adjust_tile_position,1,length(config.markers));
end

if length(config.adjust_tile_shading) == 1
    config.adjust_tile_shading = repmat(config.adjust_tile_shading,1,length(config.markers));
end

% Check selection
if isequal(config.adjust_intensity,"false")
    % No intensity adjustments
    fprintf("%s\t No intensity adjustments selected \n",datetime('now'));
    config.adj_params = [];
    [~, config] = check_adj_parameters([],config);
    return
elseif ~isequal(config.adjust_intensity, "true") && ~isequal(config.adjust_intensity, "update")
    error("Unrecognized selection for adjust_intensity. "+...
    "Valid options are ""true"", ""update"", or ""false"".")
end

adj_file = fullfile(config.output_directory,"variables",'adj_params.mat');
loaded_params = false;

if isempty(config.update_intensity_channels) && isequal(config.adjust_intensity,"update")
    config.update_intensity_channels = 1:length(config.markers);
end

% Check for previously calculated parameters
if exist(adj_file,"file") == 2 && length(config.update_intensity_channels) ~= length(config.markers)
    fprintf("%s\t Loading adjustment parameters \n",datetime('now'));
    load(fullfile(config.output_directory,'variables','adj_params.mat'),'adj_params')
    [adj_params, config] = check_adj_parameters(adj_params,config,nrows,ncols);
    
    % Attach to config structure
    config.adj_params = adj_params;
    
    % If all are loaded can return 
    if isequal(config.adjust_intensity, "true")
        return  
    else
        loaded_params = true;
    end
end

if ~loaded_params || isequal(config.adjust_intensity,"update")
    % Create new adj_params structure for all channels
    config.update_intensity_channels = 1:length(config.markers);

    % Define intensity adjustment parameters from scratch. Intensity
    % adjustment measurements should be made on raw images.
    fprintf("%s\t Defining new adjustment parameters \n",datetime('now'));        
    adj_fields = {'adjust_tile_shading','adjust_tile_position',...
        'lowerThresh','upperThresh','signalThresh'};

    % Update parameters from config structure
    for i = 1:length(adj_fields)
        params.(adj_fields{i}) = config.(adj_fields{i});
    end
end

% Calculate intensity adjustments
 for k = 1:length(config.markers)            
    if ~ismember(k,config.update_intensity_channels)
        continue
    else
        marker = config.markers(k);
    end
    fprintf("%s\t Measuring intensity for %s \n",datetime('now'),config.markers(k));
    stack = path_table(path_table.markers == config.markers(k),:);
    
    % Save marker and img_location
    if length(config.img_directory) == length(config.markers)
        params.img_directory = config.img_directory(k);
    else
        params.img_directory = config.img_directory;
    end

    % Take measurements    
    [lowerThresh_measured, upperThresh_measured, signalThresh_measured, t_adj, y_adj, flatfield, darkfield] = ...
        measure_images(config, stack, k);

    % Save adjustments to parameter structure
    params.lowerThresh = lowerThresh_measured;
    params.upperThresh = upperThresh_measured;
    params.signalThresh = signalThresh_measured;
    
    params.y_adj = y_adj;
    params.t_adj = t_adj;
    params.flatfield = flatfield;
    params.darkfield = darkfield;
    params.darkfield_intensity = config.darkfield_intensity(k);
    
    % Save params to adj params structure
    adj_params.(marker) = params;
end

% Check for user-defined intensity threshold values
[adj_params, config] = check_adj_parameters(adj_params, config, nrows, ncols);
config.adj_params = adj_params;

% Reorganize adjustment parameters and save to output directory
fprintf("%s\t Saving adjustment parameters \n", datetime('now'));
save(fullfile(config.output_directory,'variables','adj_params.mat'), 'adj_params')
config.adjust_intensity = "true";

% Save flatfield and darkfield as seperate variables
if isequal(config.adjust_tile_shading,"basic")
    fprintf("%s\t Saving flatfield and darkfield images \n",datetime('now'));
    save(fullfile(config.output_directory,'variables','flatfield.mat'), 'flatfield')
    save(fullfile(config.output_directory,'variables','darkfield.mat'), 'darkfield')
end

% Save samples
if isequal(config.save_samples,"true")
    fprintf('%s\t Saving samples \n',datetime('now'));    
    save_samples(config,'intensity',path_table)
end

end

function [config, path_table] = perform_channel_alignment(config,path_table,equal_res)
% Channel Alignment

% Check selection
if isequal(config.channel_alignment,"false")
    fprintf("%s\t No channel alignment selected. \n",datetime('now'));
    return
elseif ~isequal(config.channel_alignment, "true") && ~isequal(config.channel_alignment, "update")
    error("Unrecognized selection for channel_alignment. "+...
    "Valid options are ""true"", ""update"", or ""false"".")
end

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
nb_tiles = ncols * nrows;
position_mat = reshape(1:nb_tiles,[ncols,nrows])';

% Check if z resolution is consistent for all channels
if ~equal_res
    error("NuMorph currently does not support multi-channel alignment at "+...
            "multiple z resolutions")
end

if isequal(config.align_method,'translation')
    fprintf("%s\t Aligning channels by translation \n",datetime('now'));

    % Determine z displacement from reference channel
    % Check first if .mat file exists in output directory
    if exist(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'file') == 2 && ...
            ~isequal(config.update_z_adjustment,"true")
        fprintf("%s\t Loading z displacement matrix \n",datetime('now'));
        load(fullfile(config.output_directory,'variables','z_displacement_align.mat'),'z_displacement_align')
    else
        for k = 2:length(config.markers)
            fprintf("%s\t Performing channel alignment z calculation for %s/%s \n",...
                datetime('now'),config.markers{k},config.markers{1});
            z_tile = zeros(1,numel(position_mat));
            ave_signal = zeros(1,numel(position_mat));
            for idx = 1:numel(position_mat)
                [y,x] = find(position_mat==idx);
                path_ref = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers(1),:);
                path_mov = path_table(path_table.x==x & path_table.y==y & path_table.markers == config.markers(k),:);
                % This measures the displacement in z for a given channel to the reference
                [z_tile(idx),ave_signal(idx)] = z_align_channel(config,path_mov,path_ref,k);
                fprintf("%s\t Predicted z displacement of %d for tile %d x %d \n",...
                    datetime('now'),z_tile(idx),y,x);
            end
            z_matrix = reshape(z_tile,[ncols, nrows])';
            z_displacement_align.(config.markers(k)) = z_matrix;
        end
        % Save displacement variable to output directory
        fprintf("%s\t Saving z displacement matrix \n",datetime('now'));
        save(fullfile(config.output_directory,'variables','z_displacement_align.mat'), 'z_displacement_align')
    end

    % Create variable
    alignment_table = cell(ncols,nrows);
    save_path = fullfile(config.output_directory,'variables','alignment_table.mat');

    % Define which tiles to align if only doing subset
    tiles_to_align = 1:nb_tiles;
    if ~isempty(config.align_tiles) && all(ismember(config.align_tiles,tiles_to_align))
        tiles_to_align = config.align_tiles;
    elseif ~isempty(config.align_tiles) && ~all(ismember(config.align_tiles,tiles_to_align))
        error("Selected subset of tiles to align outside range of all tiles")
    end

    % Perform channel alignment
    for idx = tiles_to_align
        [y,x] = find(position_mat==idx);            
        fprintf("%s\t Aligning channels to %s for tile %d x %d \n",...
                    datetime('now'),config.markers{1},y,x);
        path_align = path_table(path_table.x==x & path_table.y==y,:);
        alignment_table{y,x} = align_by_translation(config,path_align,z_displacement_align);
        save(save_path,'alignment_table')
    end
    
    fprintf("%s\t Alignment completed! \n",datetime('now'));

    % Change image directory to aligned directory so that subsequent
    % steps load these images. Otherwise point config to alignment params
    if isequal(config.save_images,"true")
        config.img_directory = fullfile(config.output_directory,"aligned");
        path_table = path_to_table(config,"aligned");
    else
        config.load_alignment_params = "true";
    end
    
elseif isequal(config.align_method,'elastix')
    fprintf("%s\t Aligning channels using B-splines \n",datetime('now'));

    % Warn user if save_images flag is off. Can't apply params during
    % stitching
    if isequal(config.save_images,"false") && isequal(config.load_alignment_params,"false")
        warning("save_image flag is set to ""false"". Images need to be saved "+...
            "for downstream processing after elastix alignment.")
        pause(5)            
    end

    % Create variable
    alignment_params = cell(nrows,ncols);
    save_path = fullfile(config.output_directory,'variables','alignment_params.mat');
    if ~exist(save_path,'file')
        save(save_path,'alignment_params','-v7.3')
    else
        load(save_path,'alignment_params')
        assert(size(alignment_params,1) == nrows, "Number of row positions do not "+...
            "match loaded alignment parameters")
        assert(size(alignment_params,2) == ncols, "Number of column positions do not "+...
            "match loaded alignment parameters")
    end

    % Define which tiles to align if only doing subset
    tiles_to_align = 1:nb_tiles;
    if ~isempty(config.align_tiles) && all(ismember(config.align_tiles,tiles_to_align))
        tiles_to_align = config.align_tiles;
    elseif ~isempty(config.align_tiles) && ~all(ismember(config.align_tiles,tiles_to_align))
        error("Selected subset of tiles to align outside range of all tiles")
    end

    % Check for intensity thresholds which are required for alignment
    if isempty(config.lowerThresh) || isempty(config.upperThresh) || isempty(config.signalThresh)
        config = check_for_thresholds(config,path_table);
    end

    % If pre-aligning, check for and load calculated translation parameters
    alignment_table = [];
    if isequal(config.pre_align,"true")
        if exist(fullfile(config.var_directory,'alignment_table.mat'),'file') == 2
            load(fullfile(config.var_directory,'alignment_table.mat'),'alignment_table')
            if any(any(cellfun(@(s) isempty(s), alignment_table(position_mat)),'all'))
                idx = find(cellfun(@(s) isempty(s), alignment_table(position_mat)));
                warning("Missing alignment tables at tiles %s.",...
                    string(regexprep(num2str(idx'),' ,','')))
                pause(1)
            end
        else
            fprintf("%s\t Could not locate alignment table for pre-alignment. Generating "+...\
                "new table for pre-alignment.\n",datetime('now'));
            config_pre = config;
            config_pre.align_method = "translation";
            config_pre.align_tiles = [];
            config_pre.save_images = "false";
            config_pre.only_pc = "true";
            perform_channel_alignment(config_pre,path_table,equal_res);    
            load(fullfile(config.var_directory,'alignment_table.mat'),'alignment_table')
        end
    end

    % Perform alignment
    for i = tiles_to_align
        [y,x] = find(position_mat==i);
        path_align = path_table(path_table.x==x & path_table.y==y,:);
        if ~isempty(alignment_table)
            config.alignment_table = alignment_table{y,x};
        end
        fprintf("%s\t Aligning channels to %s for tile 0%dx0%d \n",...
            datetime('now'),config.markers{1},y,x);

        alignment_results = elastix_channel_alignment(config,path_align,true);
        
        % Save with the exception of only aligning slices
        if ~isequal(config.channel_alignment,"update") && ~isempty(config.align_slices)
            alignment_params{y,x} = alignment_results;
            save(save_path,'alignment_params','-v7.3')
        end

        % Save samples
        if isequal(config.save_samples,"true")
            fprintf('%s\t Saving samples \n',datetime('now'));
            save_samples(config,'alignment',path_align)
        end
    end

    % Change image directory to aligned directory so that subsequent
    % steps load these images
    config.img_directory = fullfile(config.output_directory,'aligned');
    path_table = path_to_table(config,"aligned",false,true);

    % Update tile intensity adjustments using newly aligned images.
    % Also, set light sheet width adjustments + flatfield adjustments
    % to false as these were applied during the alignment step
    config.adjust_tile_shading = repmat("false",1,length(config.markers));
    config.adj_params.adjust_tile_shading = config.adjust_tile_shading;

   % In case applying tile adjustments, re-calculate thresholds
    if isequal(config.adjust_intensity,"true") && isequal(config.adjust_tile_position,"true")
        for k = 1:length(config.markers)            
            fprintf("%s\t Updating intensity measurements for marker %s using "+...
                "newly aligned images \n",datetime('now'),config.markers(k));
            stack = path_table(path_table.markers == config.markers(k),:);
            [~,~,~,t_adj] = measure_images(config, stack, k);
            config.adj_params.t_adj{k} = t_adj;
        end
    end
else
    error("%s\t Unrecognized selection for align_method. "+...
        "Please select ""translation"", ""elastix"".\n",string(datetime('now')))
end

% Check number of loaded tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
if nb_tiles ~= ncols*nrows
    warning("Number of aligned tiles and raw image tiles are not equal")
    pause(5)
end

end

function [config, path_table] = perform_stitching(config,path_table)
% Stitching

% Check selection
if isequal(config.stitch_images,"false")
    fprintf("%s\t No image stitching selected. \n",datetime('now'));
    return
elseif ~isequal(config.stitch_images, "true") && ~isequal(config.stitch_images, "update")
    error("Unrecognized selection for stitch_images. "+...
    "Valid options are ""true"", ""update"", or ""false"".")
end

% Count number of x,y tiles
ncols = length(unique(path_table.x));
nrows = length(unique(path_table.y));
fprintf("%s\t Perfoming iterative 2D image stitching \n",datetime('now'));

% Check for aligned images
if isequal(config.img_directory, fullfile(config.output_directory,'aligned'))
    fprintf("%s\t Using images from aligned directory \n",datetime('now'));
    config.alignment_table = [];
elseif isequal(config.load_alignment_params,"true")
    fprintf("%s\t Loading channel alignment translations \n",datetime('now'));
    load(fullfile(config.output_directory,'variables','alignment_table.mat'),'alignment_table')

    % Check tiles
    empty_alignment_tiles = find(cellfun(@(s) isempty(s),alignment_table'));
    if empty_alignment_tiles
        error("No alignment parameters to apply for tiles %s",...
            sprintf("%s",num2str(empty_alignment_tiles)))
    end
    config.alignment_table = alignment_table;
else
    config.alignment_table = [];
end

% Check for intensity thresholds which are required for alignment
if isempty(config.lowerThresh) || isempty(config.upperThresh) || isempty(config.signalThresh)
    config = check_for_thresholds(config,path_table);
end

% Trim markers to only the ones specified
if ~isempty(config.stitch_sub_channel)
   markers = config.markers(config.stitch_sub_channel);
   path_table = path_table(ismember(path_table.markers,markers),:);
end

% Trim z positions that aren't present in all tiles
min_z_mat = zeros(nrows,ncols);
max_z_mat = zeros(nrows,ncols);
for i = 1:ncols
    for j = 1:nrows
        max_z_mat(j,i) = max(path_table.z(path_table.x==i & path_table.y==j));
        min_z_mat(j,i) = min(path_table.z(path_table.x==i & path_table.y==j));
    end
end
max_z = min(max_z_mat(:));
min_z = max(min_z_mat(:));
path_table(path_table.z>max_z,:)=[];
path_table(path_table.z<min_z,:)=[];

% Calculate adjustments in z
if isequal(config.update_z_adjustment,"true") || exist(fullfile(config.var_directory,'z_disp_matrix.mat'),'file') ~= 2
    % Calculate z displacements
    fprintf("%s\t Calculating z displacement matrix for stitching \n",datetime('now'));
    z_adj = calculate_adjusted_z(path_table,nrows,ncols,config.markers,...
            config.overlap,config.z_positions,config.z_window,config.lowerThresh,...
            config.output_directory);
    [~,z_idx] = setdiff(path_table.file,z_adj.file);
    path_table(z_idx,:) = [];
    path_table.z_adj = z_adj.z_adj;
else
    % Load z displacement info from a matrix
    fprintf("%s\t Loading z displacement matrix for stitching \n",datetime('now'));
    load(fullfile(config.var_directory, 'z_disp_matrix.mat'), 'z_disp_matrix')
    assert(nrows == size(z_disp_matrix,1) && ncols == size(z_disp_matrix,2), "Loaded z adjustment parameters do not match number of tiles "+...
        "detected in the input image directory. Check input image directory or update z adjustment")
    z_adj = apply_adjusted_z(path_table, z_disp_matrix);
    [~,z_idx] = setdiff(path_table.file,z_adj.file);
    path_table(z_idx,:) = []; 
    path_table.z_adj = z_adj.z_adj;
end

% Check whether to stitch from previously calculated transforms.
nb_images = length(min(path_table.z_adj):max(path_table.z_adj));
var_file = fullfile(config.output_directory,'variables','stitch_tforms.mat');
h_stitch_tforms = []; v_stitch_tforms = [];
if ~isequal(config.stitch_images,"update")
    if exist(var_file,'file') == 2
        fprintf("%s\t Loading previously calculated stitching parameters \n",datetime('now'));
        load(var_file, 'h_stitch_tforms','v_stitch_tforms')
        
        % Check if sizes match up
        if size(h_stitch_tforms,2) ~= nb_images
           error("%s\t Number of slices z slices in stitching parameters does not match loaded images \n",string(datetime('now')))
        elseif size(h_stitch_tforms,1) ~= (ncols-1)*nrows*2
           error("%s\t Number of column positions does not match loaded stitching parameters \n",string(datetime('now')))
        elseif size(v_stitch_tforms,1) ~= (nrows-1)*2
           error("%s\t Number of row positions does not match loaded stitching parameters \n",string(datetime('now')))
        end
    end
end

if isempty(h_stitch_tforms)
    % Perform stitching
    stitch_iterative(config,path_table)
else
    % Perform stitching using loaded parameters
    stitch_from_loaded_parameters(path_table, h_stitch_tforms, v_stitch_tforms, config)
end

fprintf("%s\t Stitching completed! \n",datetime('now'));

end