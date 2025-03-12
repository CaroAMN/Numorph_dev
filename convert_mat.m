
function convertMatFile(matFilePath)
     
    outputPath = '';

    try
        % Load the .mat file
        loadedData = load(matFilePath);
        
        % Get field names
        dataFieldNames = fieldnames(loadedData);
        
        % Process single field
        if length(dataFieldNames) == 1
            data = loadedData.(dataFieldNames{1});
            
            % Handle table data
            if istable(data)
                [filePath, fileName, ~] = fileparts(matFilePath);
                outputPath = fullfile(filePath, [fileName, '.csv']);
                writetable(data, outputPath);

            % Handle structure data
            else
                [filePath, fileName, ~] = fileparts(matFilePath);
                outputPath = fullfile(filePath, [fileName, '.json']);
                
                % Convert to JSON
                if isstruct(data)
                    jsonStr = jsonencode(data);
                else
                    jsonStr = jsonencode(data);
                end
                
                % Write JSON file
                writeJSON(outputPath, jsonStr);
            end
            
        % Process multiple fields
        elseif length(dataFieldNames) > 1
            % Combine all fields into one structure
            allData = struct();
            for idx = 1:length(dataFieldNames)
                fieldName = dataFieldNames{idx};
                allData.(fieldName) = loadedData.(fieldName);
            end
            
            % Save as JSON
            [filePath, fileName, ~] = fileparts(matFilePath);
            outputPath = fullfile(filePath, [fileName, '.json']);
            writeJSON(outputPath, jsonencode(allData));
            
        else
            warning('File contains no variables: %s', matFilePath);
            return;
        end
        
        % Set success status
        
        fprintf('Converted %s to %s\n', matFilePath, outputPath);
        
    catch ME
        warning('Failed to convert file %s: %s', matFilePath, ME.message);
        return;
    end
end

function writeJSON(filepath, jsonStr)
    fid = fopen(filepath, 'w');
    if fid == -1
        error('Cannot create JSON file: %s', filepath);
    end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
end
