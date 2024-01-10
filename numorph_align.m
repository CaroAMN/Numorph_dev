function align = numorph_align(img_directory, output_directory, config_file, sample, use_processed_images)
    NM_setup;

    % read in the config file which is a csv file
    % convert config file to a structure like in original numorph
    % Step 1: Read the CSV file into a table
    %process_template = readtable(config_file, "FileType","text",'Delimiter', '\t');
    % Step 2: Convert the table to a structure
    tsv_to_mat_file(config_file, img_directory, output_directory, sample, use_processed_images);

    % Additional processing (if necessary)
    config = NM_config('process', sample);
    NM_process(config, "align", sample);

end

function tsv_to_mat_file(tsv_file_path, img_directory, output_directory, sample, use_processed_images)
    % Define the output directory and file name
    output_dir = fullfile(pwd, 'data', 'tmp');
    mat_file_name = 'NM_variables.mat';
    mat_file_path = fullfile(output_dir, mat_file_name);

    % Create the output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Read the TSV file into a table
    tsv_data = readtable(tsv_file_path, 'FileType', 'text', 'Delimiter', '\t');

    % Create a structure to store the variables
    S = struct();

    % Loop through each row of the table and convert the data
    for i = 1:height(tsv_data)
        parameter = tsv_data.Parameter{i};
        value = tsv_data.Value{i};
        % Handle different data types for value
        S.(parameter) = convert_to_proper_type(value);    
    end

    % Add img_directory, output_directory, and sample to the structure
    % Ensure they are stored as strings with double quotes
    S.img_directory = string(img_directory);
    S.output_directory = string(output_directory);
    S.sample_id = string(sample);
    S.use_processed_images = string(use_processed_images);

    % Save the structure as a .mat file
    save(mat_file_path, '-struct', 'S');
    fprintf('The .mat file has been created at: %s\n', mat_file_path);
end

function out = convert_to_double_quoted_string(val)
    if ischar(val) || isstring(val)
        fprintf(val);
        % If the value is already a string with double quotes, do not add extra quotes
        if startsWith(val, '"') && endsWith(val, '"')
            out = val;  % Keep as is
        else
            % Add double quotes around the string
            out = ['"' val '"'];
        end
    else
        out = val;  % Non-string values are returned as is
    end
end






function out = convert_to_proper_type(val)
    % Check if the value looks like an array
    if startsWith(val, '[') && endsWith(val, ']')
        % Remove the brackets
        val = regexprep(val, '^\[|\]$', '');

        % Split by commas
        elements = strsplit(val, ',');
        elements = strtrim(elements); % Remove extra whitespace

        % Determine if the elements are numeric, regex patterns, or other strings
        numericElements = cellfun(@str2double, elements);
        if all(~isnan(numericElements))
            out = numericElements; % Array of numbers
        else
            % Special handling for 'position_exp' to save as a string array
            if contains(val, '\')
                % Convert to string array for regex patterns
                out = string(elements);
                % Remove double quotes from regex patterns, if present
                out = regexprep(out, '^""|""$', '');
            else
                % For other non-numeric arrays, use cell array of strings
                out = elements;
            end
        end
    else
        % Convert single values
        numericVal = str2double(val);
        if ~isnan(numericVal)
            out = numericVal;
        elseif strcmp(val, 'true') || strcmp(val, 'false')
            out = strcmp(val, 'true');
        else
            out = val; % Keep as string, remove extra quotes if present
            if isstring(out)
                out = regexprep(out, '^""|""$', '');
            end
        end
    end
end


function cellArray = regex_patterns_to_cell(str)
    str = strtrim(regexprep(str, '^\[|\]$', '')); % Remove square brackets
    cellArray = strsplit(str, ',');
    cellArray = cellfun(@strtrim, cellArray, 'UniformOutput', false);
    % Remove double quotes from each element, if present
    cellArray = cellfun(@(x) regexprep(x, '^"|"$', ''), cellArray, 'UniformOutput', false);
    % Convert elements to char if they are strings
    cellArray = cellfun(@char, cellArray, 'UniformOutput', false);
end






    
    %{
    process_template = convert_table_types(process_template);
  
   % Loop through the table rows and save each parameter-value pair to the .mat file
    for i = 1:height(process_template)
        % Convert the parameter name to a valid field name if necessary
        fieldName = matlab.lang.makeValidName(process_template{i, 1}{1});
        fieldValue = process_template{i, 2}{1};
        save_to_mat_file(fieldName, fieldValue, fullfile('data', 'tmp', 'NM_variables.mat'));
    end

    % Add additional fields to the .mat file
    save_to_mat_file('img_directory', img_directory, fullfile('data', 'tmp', 'NM_variables.mat'));
    save_to_mat_file('output_directory', output_directory, fullfile('data', 'tmp', 'NM_variables.mat'));
    save_to_mat_file('sample', sample, fullfile('data', 'tmp', 'NM_variables.mat'));

    % Additional processing (if necessary)
    config = NM_config('process', sample);
    NM_process(config, "align", sample);
end

function save_to_mat_file(varName, varValue, filename)
    % This function appends a single variable to the .mat file
    S.(varName) = varValue;
    if exist(filename, 'file') == 2 % If the file exists, append the variable
        save(filename, '-struct', 'S', varName, '-append');
    else % If the file does not exist, create it
        save(filename, '-struct', 'S', varName);
    end
end

function outTable = convert_table_types(inTable)
    outTable = inTable;
    for i = 1:width(outTable)
        for j = 1:height(outTable)
            val = outTable{j, i}{1};
            if isempty(val)
                continue; % Skip empty values
            end
            numericVal = str2double(val);
            if ~isnan(numericVal)
                outTable{j, i} = numericVal; % Assign numeric values directly
            elseif startsWith(val, '[') && endsWith(val, ']')
                % Check for numeric array or string array
                numericArray = str2num(val); %#ok<ST2NM>
                if ~isempty(numericArray)
                    outTable{j, i} = numericArray; % Assign numeric arrays directly
                else
                    strArray = str2cell(val);
                    if ~isempty(strArray)
                        % Convert strArray to a MATLAB string array
                        outTable{j, i} = string(strArray); % Convert to MATLAB string array
                    end
                end
            elseif strcmpi(val, 'true') || strcmpi(val, 'false')
                outTable{j, i} = strcmpi(val, 'true'); % Assign logical values directly
            else
                outTable{j, i} = string(val); % Assign strings directly
            end
        end
    end
end

function cellArray = str2cell(str)
    str = strtrim(str);
    if isempty(str) || str == "[]"
        cellArray = {}; % Return empty cell array for empty input
    else
        str = regexprep(str, '^\[|\]$', ''); % Remove square brackets
        cellArray = regexp(str, ',', 'split');
        cellArray = strtrim(cellArray); % Remove whitespace
    end
end
%}