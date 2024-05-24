function preprocessing = numorph_preprocessing(varargin)
    %Required: --input_dir', @ischar);
    %Required: --output_dir', @ischar);
    %Required: --parameter_file', @ischar);
    %Required: --sample_name', @ischar);
    %Required: --stage', @ischar);



    % step: possible steps to run for preprocessing 
    % 'process' : runs complete preprocessing pipeline
    % 'intensity' : runs intensity adjustment step
    % 'align' : runs alignment step
    % 'stitch' : runs stitching step

    %call the parser function to parse the input arguments
    params = parseInputArgs(varargin{:});

    % Get the values of the parameters
    input_dir = params.input_dir;
    output_dir = params.output_dir;
    parameter_file = params.parameter_file;
    sample_name = params.sample_name;
    stage = params.stage;

    % Run the preprocessing pipeline

    %addpath(genpath(pwd)); % do i need this ?
    NM_setup; % setup numorph

    % read in the config file which is a csv file
    % convert config file to a structure like in original numorph
    home_path = fileparts(which('NM_config'));
    tmp_folder = fullfile(home_path,'data','tmp', 'NM_variables.mat');

    csv_to_mat_file(parameter_file, tmp_folder, input_dir, output_dir);

    % old output directory path from the prevoius version of the wrapper 
    %output_dir = load(tmp_folder, 'output_directory');
    

    switch stage
        case 'process'
            config = NM_config('process', sample_name);
            NM_process(config, "process");
        case 'intensity'
            config = NM_config('process', sample_name);
            NM_process(config, "intensity");
        case 'align'
            config = NM_config('process', sample_name);
            NM_process(config, "align");
        case 'stitch'
            config = NM_config('process', sample_name);
            NM_process(config, "stitch");
        otherwise
            error('Invalid step. Please choose from: process, intensity, align, stitch');
    end

    % Specify the directory where the .mat files are located/ old 
    %directoryPath = output_dir.output_directory; / old
    

    matFilePaths = getAllMatFiles(output_dir);
    

    % Loop through each .mat file
    for i = 1:length(matFilePaths)
        % Construct the full path to the .mat file
        matFilePath = matFilePaths{i};
        disp(matFilePath);
    
        % Load the .mat file
        loadedData = load(matFilePath);
    
        % Assuming the .mat file contains a single variable or a structure.
        % You might need to adjust this part if your .mat files are different.
        dataFieldNames = fieldnames(loadedData);
        if length(dataFieldNames) == 1
            data = loadedData.(dataFieldNames{1});
        
            % Check if the data is a table, if not, you might need to convert or handle differently
            if istable(data)
                % Construct the path for the .csv file
                % It uses the same base filename as the .mat file
                [filePath, fileName, ~] = fileparts(matFilePath);
                csvFilePath = fullfile(filePath, [fileName, '.csv']);
            
                % Write the table to a .csv file
                writetable(data, csvFilePath);
                fprintf('Converted %s to %s\n', matFilePaths{i}, [fileName, '.csv']);
            elseif isstruct(data)
                % Convert structure to JSON string
                jsonStr = jsonencode(data);
                % Construct the path for the output JSON file
                [filePath, fileName, ~] = fileparts(matFilePath);
                jsonFilePath = fullfile(filePath, [fileName, '.json']);
                
                fid = fopen(jsonFilePath, 'w');
                if fid == -1
                    error('Cannot create JSON file');
                end
                fwrite(fid, jsonStr, 'char');
                fclose(fid);
                fprintf('Converted %s to %s\n', matFilePaths{i}, [fileName, '.json']);
            else
                % Convert structure to JSON string
                jsonStr = jsonencode(data);
                % Construct the path for the output JSON file
                [filePath, fileName, ~] = fileparts(matFilePath);
                jsonFilePath = fullfile(filePath, [fileName, '.json']);
                
                fid = fopen(jsonFilePath, 'w');
                if fid == -1
                    error('Cannot create JSON file');
                end
                fwrite(fid, jsonStr, 'char');
                fclose(fid);
                fprintf('Converted %s to %s\n', matFilePaths{i}, [fileName, '.json']);
            end
        % Assuming the the .mat file contains a structure with multiple fields    
        elseif length(dataFieldNames) >= 1
            % Loop through each field if you have multiple fields in your structure
            allData = struct();
            for idx = 1:length(dataFieldNames)
                fieldName = dataFieldNames{idx};
                allData.(fieldName) = loadedData.(fieldName);
            end
            % Convert the structure to a JSON string
            jsonData = jsonencode(allData);
            
            % Define the JSON file path
            [filePath, fileName, ~] = fileparts(matFilePath);
            jsonFilePath = fullfile(filePath, [fileName, '.json']);
            %jsonFilePath = 'converted_data.json';
            
            % Write the JSON string to a file
            fileId = fopen(jsonFilePath, 'w');
            if fileId == -1
                error('Failed to create JSON file: %s', jsonFilePath);
            end
            fwrite(fileId, jsonData, 'char');
            fclose(fileId);
            disp(['Data converted and saved to JSON: ', jsonFilePath]);
        else
            warning('%s contains no variables and was not converted.', matFilePaths{i});
        end
    end

end



function csv_to_mat_file(csvFilePath, matFilePath, input_dir, output_dir)
    % Read the CSV file into a table
    tbl = readtable(csvFilePath, 'ReadVariableNames', true, 'Delimiter', ',');

    % Create a structure to hold the variables
    S = struct();
    S.img_directory = string(input_dir);
    S.output_directory = string(output_dir);

    % Iterate over the rows of the table and add to the structure
    for i = 1:height(tbl)
        parameter = tbl.Parameter{i};
        value = tbl.Value{i};
        % save vars to structure according to original numorph 

        % "adjust_intensity": true, update, false; Whether to calculate and apply any of the following intensity adjustments. Intensity adjustment measurements should typically be performed on raw images
        if (parameter == "adjust_intensity")
            S.adjust_intensity = string(value);

        % true, false; Can be 1xn_channels. Normalize tile intensities by position using overlapping regions    
        elseif (parameter == "adjust_tile_position")
            S.adjust_tile_position = string(value);

        % basic, manual, false; Can be 1xn_channels. Perform shading correction using BaSIC algorithm or using manual measurements from UMII microscope    
        elseif (parameter == "adjust_tile_shading")
            S.adjust_tile_shading = string(value);

        % Option to align only certain channels (set to >1)    
        elseif (parameter == "align_channels")
            if isnan(str2double(value))
                S.align_channels = double([]);
            else
                S.align_channels = str2double(value);
            end

        % Only for alignment by elastix. Option to align only certain chunks    
        elseif (parameter == "align_chunks")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.align_chunks = double([]);
            else
                S.align_chunks = str2double(value);
            end

        % elastix, translation; Channel alignment by rigid, 2D translation or non-rigid B-splines using elastix    
        elseif (parameter == "align_method")
            S.align_method = string(value);

        % Option to align only certain slice ranges. Set as cell array for non-continuous ranges (i.e. {1:100,200:300})    
        elseif (parameter == "align_slices")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.align_slices = double([]);
            else
                S.align_slices = cell(value);
            end

        % interger; Only for alignment by translation. Number of images sampled for determining translations. Images in between are interpolated    
        elseif (parameter == "align_stepsize")
            S.align_stepsize = str2double(value);

        % Option to align only certain stacks and not all stacks. Row-major order
        elseif (parameter == "align_tiles")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.align_tiles = double([]);
            else
                S.align_tiles = str2double(value);
            end

        % sigmoid, linear, max    
        elseif (parameter == "blending_method")
            S.blending_method = string(value); 

        % integer >= 0; Crops borders during stitching. Increase if images shift significantly between channels to prevent zeros values from entering stitched image    
        elseif (parameter == "border_pad")
            S.border_pad = str2double(value);

        % "channel_alignment": true, update, false; 
        elseif (parameter == "channel_alignment")
            S.channel_alignment = string(value);

        % Channel id    
        elseif (parameter == "channel_num")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.channel_num = string([]);
            else
                S.channel_num = string(value);
            end

        % integer; Padding around chunks. Should be set to value greater than the maximum expected translation in z    
        elseif (parameter == "chunk_pad")
            S.chunk_pad = str2double(value);

        % 1xn_channels; Constant darkfield intensity value (i.e. average intensity of image with nothing present)        
        elseif (parameter == "darkfield_intensity")
            S.darkfield_intensity = str2double(value);

        % [0,1]; Factor controlling amount of adjustment to apply. Set to 1 for absolute DoG    
        elseif (parameter == "DoG_factor")
            S.DoG_factor = str2double(value);

        % true,false; Apply difference of gaussian enhancement of blobs    
        elseif (parameter == "DoG_img")
            S.DoG_img = string(value);

        % 1x2 numeric; Min/max sigma values to take differene from.    
        elseif (parameter == "DoG_minmax")
            value = split(value, ';');
            S.DoG_minmax = str2double(value);

        % 1xn_channels-1 string; Name of folders containing elastix registration parameters. Place in /supplementary_data/elastix_parameter_files/channel_alignment    
        elseif (parameter == "elastix_params")
            S.elastix_params = string(value);


         % "none", "horizontal", "vertical", "both"; Flip image along horizontal or vertical axis    
        elseif (parameter == "flip_axis")
            S.flip_axis = string(value);

         % 1xn_channels numeric; Gamma intensity adjustment    
        elseif (parameter == "Gamma")
            S.Gamma = str2double(value);

        % Group name/id    
        elseif (parameter == "group")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.group = string([]);
            else
                S.group = string(value);
            end

        % "left","right","both","none"    
        elseif (parameter == "hemisphere")
            S.hemisphere = string(value);

        % 1xn_channels-1 interger; Match histogram bins to reference channel? If so, specify number of bins. Otherwise leave empty or set to 0. This can be useful for low contrast images    
        elseif (parameter == "hist_match")
            S.hist_match = str2double(value);

        % completely ignore marker from processing steps.    
        elseif (parameter == "ignore_markers")
            S.ignore_markers = string(value);
        
        % Input image directory path
        %elseif (parameter == "img_directory")
            %S.img_directory = string(value);

        % [-0.5,0.5]; Displacement of light-sheet along y axis. Value of 0.5 means light-sheet center is positioned at the top of the image    
        elseif (parameter == "laser_y_displacement")
            value = split(value, ';');
            S.laser_y_displacement = str2double(value);

        % true, false; Apply channel alignment translations during stitching    
        elseif (parameter == "load_alignment_params")
            S.load_alignment_params = string(value);

         % 1xn_channels numeric; Lower intensity for rescaling    
        elseif (parameter == "lowerThresh")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.lowerThresh = double([]);
            else
                S.lowerThresh = str2double(value);
            end

        % 1xn_channels interger; Light sheet width setting for UltraMicroscope II as percentage    
        elseif (parameter == "ls_width")
            S.ls_width = str2double(value);

        % Name of markers present    
        elseif (parameter == "markers")
            value = split(value, ';');
            
            if (strlength(value) == 0)
                S.markers = string([]);
            else
                S.markers = string(value);
            end

        % numeric; Mask intensity threshold for choosing signal pixels in elastix channel alignment. Leave empty to calculate automatically    
        elseif (parameter == "mask_int_threshold")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.mask_int_threshold = double([]);
            else
                S.mask_int_threshold = str2double(value);
            end

        % integer; Chunk size for elastix alignment. Decreasing may improve precision but can give spurious results    
        elseif (parameter == "max_chunk_size")
            S.max_chunk_size = str2double(value);

        % numeric >= 1; Max radius of cell nuclei along x/y in pixels. Required also for DoG filtering    
        elseif (parameter == "nuc_radius")
            S.nuc_radius = str2double(value);

        % true, false; Use only phase correlation for registration. This gives only a quick estimate for channel alignment.     
        elseif (parameter == "only_pc")
            S.only_pc = string(value);

        % 1x3 string specifying sample orientation
        % Orientation key: anterior(a)/posterior(p), superior(s)/inferior(i), left(l)/right(r)
        elseif (parameter == "orientation")
            S.orientation = string(value);

        % Directory to save results    
        %elseif (parameter == "output_directory")
            %S.output_directory = string(value);

         % 0:1; overlap between tiles as fraction    
        elseif (parameter == "overlap")
            S.overlap = str2double(value);

        % 1x3 string of regular expression specifying image row(y), column(x), slice(z)    
        elseif (parameter == "position_exp")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.position_exp = string([]);
            else
                S.position_exp = string(value);
            end

        % true, false; (Experimental) Option to pre-align using translation method prior to non-linear registration    
        elseif (parameter == "pre_align")
            S.pre_align = string(value);

        % 1x3 integer. Amount of downsampling along each axis. Some downsampling, ideally close to isotropic resolution, is recommended    
        elseif (parameter == "resample_s")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.resample_s = double([]);
            else
                S.resample_s = str2double(value);
            end

        % true, false; Rescaling intensities and applying gamma    
        elseif (parameter == "rescale_intensities")
            S.rescale_intensities = string(value);

        % Image reolution in um/voxel    
        elseif (parameter == "resolution")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.resolution = double([]);
            else
            value = str2double(value);
            value = reshape(value,[1,3]);
            S.resolution = num2cell(value,2);
            end

        % 90 or -90; Rotate image    
        elseif (parameter == "rotate_axis")
            S.rotate_axis = str2double(value);

        elseif (parameter == "sample")
            S.sample = string(value);

        elseif (parameter == "sample_id")
            S.sample_id = string(value);

        % [0,1]; Fraction of images to read and sample from. Setting to 1 means use all images    
        elseif (parameter == "sampling_frequency")
            S.sampling_frequency = str2double(value);

        elseif (parameter == "save_flag")
            S.save_flag = str2double(value);

        % true or false; Save images during processing. Otherwise only parameters will be calculated and saved    
        elseif (parameter == "save_images")
            S.save_images = string(value);

        % true, false; Save sample results for each major step    
        elseif (parameter == "save_samples")
            S.save_samples = string(value);

        % 0:1; Recommended: ~0.05. Steepness of sigmoid-based blending. Larger values give more block-like blending    
        elseif (parameter == "sd")
            S.sd = str2double(value);

        % integer vector; Subset tile positions for calculating shading correction (row major order). It's recommended that bright regions are avoid    
        elseif (parameter == "shading_correction_tiles")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.shading_correction_tiles = double([]);
            else
                S.shading_correction_tiles = str2double(value);
            end

        % numeric >= 1; Factor for adjusting the total effect of shading correction. Greater values lead to a smaller overall adjustment    
        elseif (parameter == "shading_intensity")
            S.shading_intensity = str2double(value);

        %  numeric >= 1; Factor for adjusting smoothness of shading correction. Greater values lead to a smoother flatfield image    
        elseif (parameter == "shading_smoothness")
            S.shading_smoothness = str2double(value);

        % true, false; Refine stitching using SIFT algorithm (requires vl_fleat toolbox)    
        elseif (parameter == "sift_refinement")
            S.sift_refinement = string(value);

        % 1xn_channels numeric; Rough estimate for minimal intensity for features of interest    
        elseif (parameter == "signalThresh")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.signalThresh = double([]);
            else
                S.signalThresh = str2double(value);
            end

        % true, false; Whether a single sheet was used for acquisition    
        elseif (parameter == "single_sheet")
            S.single_sheet = string(value);

        % 1xn_channels, "gaussian", "median", "guided". Apply a smoothing filter    
        elseif (parameter == "smooth_img")
            S.smooth_img = string(value);

            % 1xn_channels numeric; Size of smoothing kernel. For median and guided filters, it is the dimension of the kernel size
        elseif (parameter == "smooth_sigma")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.smooth_sigma = double([]);
            else
                S.smooth_sigma = str2double(value);
            end

        % "stitch_images": true, update, false; 2D iterative stitching
        elseif (parameter == "stitch_images")
            S.stitch_images = string(value);

        % z index; Start stitching from specific position. Otherwise this will be optimized    
        elseif (parameter == "stitch_start_slice")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.stitch_start_slice = double([]);
            else
                S.stitch_start_slice = str2double(value);
            end

        % channel index; If only stitching certain channels    
        elseif (parameter == "stitch_sub_channel")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.stitch_sub_channel = double([]);
            else
                S.stitch_sub_channel = str2double(value);
            end

        % z positions; If only stitching a cetrain z range from all the images    
        elseif (parameter == "stitch_sub_stack")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.stitch_sub_stack = double([]);
            else
                S.stitch_sub_stack = str2double(value);
            end

        % true, false. Subtrat background (similar to Fiji's rolling ball background subtraction)   
        elseif (parameter == "subtract_background")
            S.subtract_background = string(value);

        % integers; Update intensity adjustments only to certain channels    
        elseif (parameter == "update_intensity_channels")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.update_intensity_channels = double([]);
            else
                S.update_intensity_channels = str2double(value);
            end

        % true, false; Update z adjusment steps with new parameters. Otherwise pipeline will search for previously calculated parameters    
        elseif (parameter == "update_z_adjustment")
            S.update_z_adjustment = string(value);

         % 1xn_channels numeric; Upper intensity for rescaling    
        elseif (parameter == "upperThresh")
            value = split(value, ';');
            if (strlength(value) == 0)
                S.upperThresh = double([]);
            else
                S.upperThresh = str2double(value);
            end

        % "use_processed_images": false or name of sub-directory in output directory (i.e. aligned, stitched...); Load previously processed images in output directory as input images    
        elseif (parameter == "use_processed_images")
            S.use_processed_images = string(value);

        % 1xn_channels-1 interger; Predicted initial z displacement between reference channel and secondary channel (i.e.    
        elseif (parameter == "z_initial")
            value = split(value, ';');
            S.z_initial = str2double(value);

        % integer or numeric; Sampling positions along adjacent image stacks to determine z displacement. If <1, uses fraction of all images. Set to 0 for no adjustment, only if you're confident tiles are aligned along z dimension    
        elseif (parameter == "z_positions")
            S.z_positions = str2double(value);

        % integer; Search window for finding corresponding tiles (i.e. +/-n z positions)    
        else (parameter == "z_window")
            S.z_window = str2double(value)
        end
    end
    % Save the structure as a .mat file
    save(matFilePath, '-struct', 'S');
end

function fileList = getAllMatFiles(startPath)
    % Generate a path string that includes all subdirectories
    allSubDirs = genpath(startPath);
    
    % Split the path string into individual directory paths
    dirList = strsplit(allSubDirs, pathsep);
    
    % Initialize the list of file paths
    fileList = {};
    
    % Loop over all directories
    for d = 1:length(dirList)
        % Skip any empty directory paths (can happen due to the split)
        if isempty(dirList{d})
            continue;
        end
        
        % Get a list of all '.mat' files in the current directory
        matFiles = dir(fullfile(dirList{d}, '*.mat'));
        
        % Add their full paths to the list
        for k = 1:length(matFiles)
            fileList{end+1, 1} = fullfile(matFiles(k).folder, matFiles(k).name);
        end
    end
end

% function that parses commanline parameters
function params = parseInputArgs(varargin)
    % Create an input parser
    p = inputParser;

    % Define the required parameters
    addParameter(p, 'input_dir', '', @ischar);
    addParameter(p, 'output_dir', '', @ischar);
    addParameter(p, 'parameter_file', '', @ischar);
    addParameter(p, 'sample_name', '', @ischar);
    addParameter(p, 'stage', '', @ischar);

    % Parse the inputs
    parse(p, varargin{:});

    % Get the values of the parameters
    params.input_dir = p.Results.input_dir;
    params.output_dir = p.Results.output_dir;
    params.parameter_file = p.Results.parameter_file;
    params.sample_name = p.Results.sample_name;
    params.stage = p.Results.stage;

    % Check if required arguments exsit 
    if isempty (params.input_dir)
        error('Input directory is required');
    end
    if isempty (params.output_dir)
        error('Output directory is required');
    end
    if isempty (params.parameter_file)
        error('Parameter file is required');
    end
    if isempty (params.sample_name)
        error('Sample name is required');
    end
    if isempty (params.stage)
        error('Stage is required');
    end

end