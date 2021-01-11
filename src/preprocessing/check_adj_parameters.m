function [adj_params_full, config, lowerThresh, upperThresh, signalThresh] = check_adj_parameters(adj_params, config, nrows, ncols)
%--------------------------------------------------------------------------
% Check for user defined adjustment parameters that will override measured
% values.
%--------------------------------------------------------------------------
% 
% Inputs:
% adj_params - structure containing adjustment parameters.
%
% config - config structure from NM_process.
%
% Ouputs:
%--------------------------------------------------------------------------
        
% Check for user-defined intensity threshold values
lowerThresh = config.lowerThresh;
signalThresh = config.signalThresh;
upperThresh = config.upperThresh;

% If given is 16-bit integer, rescale to unit interval
if lowerThresh>1
    lowerThresh = lowerThresh/65535;
end
if upperThresh>1
    upperThresh = upperThresh/65535;
end
if signalThresh>1
    signalThresh = signalThresh/65535;
end

% If no adjustment parameters specified, check lengths and return
if isempty(adj_params)
    if ~isempty([lowerThresh upperThresh signalThresh])
        n_markers = length(config.markers);
        msg = "User thresholds must be specified for all markers.";
        assert(length(lowerThresh) == n_markers,msg)
        assert(length(upperThresh) == n_markers,msg)
        assert(length(signalThresh) == n_markers,msg)
        fprintf('%s\t Using user-defined thresholds. \n',datetime('now'));  
        config.lowerThresh = lowerThresh;
        config.upperThresh = upperThresh;
        config.signalThresh = signalThresh;
    end
    adj_params_full = [];
    return
end

% Get measured thresholds from adj_params structure
adj_params_full = adj_params;
for i = 1:length(config.markers)
    if isfield(adj_params_full,config.markers(i))
        adj_params = adj_params_full.(config.markers(i));
    else
        error("Adjustment parameters for marker %s is not found in parameter structure. "+...
            "Update adjustment parameters.",config.markers(i))
    end
    lowerThresh_measured = adj_params.lowerThresh;
    upperThresh_measured = adj_params.upperThresh;
    signalThresh_measured = adj_params.signalThresh;

    % If empty, set to predicted threshold values
    if isempty(lowerThresh)
        adj_params.lowerThresh = lowerThresh_measured;
        config.lowerThresh(i) = lowerThresh_measured;
    else
        fprintf('%s\t Using user-defined lowerThresh of %s \n',datetime('now'),num2str(round(lowerThresh*65535)));    
        adj_params.lowerThresh = lowerThresh(i);
    end

    if isempty(upperThresh)
        adj_params.upperThresh = upperThresh_measured;
        config.upperThresh(i) = upperThresh_measured;
    else
        fprintf('%s\t Using user-defined upperThresh of %s \n',datetime('now'),num2str(round(upperThresh*65535)));    
        adj_params.upperThresh = upperThresh(i);
    end

    if isempty(signalThresh)
        adj_params.signalThresh = signalThresh_measured;
        config.signalThresh(i) = signalThresh_measured;
    else
        fprintf('%s\t Using user-defined signalThresh of %s \n',datetime('now'),num2str(round(signalThresh*65535)));    
        adj_params.signalThresh = signalThresh(i);
    end

    % Set Gamma equal to 1 if undefined
    if isempty(config.Gamma)
        adj_params.Gamma = 1;
    else
        adj_params.Gamma = config.Gamma(i); 
    end

    % Check img directories
    if ~any(adj_params.img_directory == config.img_directory)
        warning("Intensity adjustment parameters were calculated from %s "+...
            "and not the input image directory. Some adjustments may have been "+...
            "already applied. Consider updating adjustment parameters.",...
            adj_params.img_directory)
        pause(5)
    end
    
    % Update which adjustment to apply based on current configs
    adj_params.adjust_tile_shading = config.adjust_tile_shading(i);
    adj_params.adjust_tile_position = config.adjust_tile_position(i);
    
    % Check tile position
    if isequal(adj_params.adjust_tile_position,'true')
        assert(~isempty(adj_params.t_adj), "No tile position adjustments found for "+...
            "marker %s",config.markers(i))
        assert(size(adj_params.t_adj,1) == nrows(i), "Incorrect numner of rows "+...
               "in the tile adjustment matrix for marker %s",config.markers(i))
        assert(size(adj_params.t_adj,2) == ncols(i), "Incorrect numner of columns "+...
               "in the tile adjustment matrix for marker %s",config.markers(i))
    end
    
    % Save into full parameter structure
    adj_params_full.(config.markers(i)) = adj_params;
end

% Set default nucleus diameter (15um)
if isempty(config.nuc_radius)
    config.nuc_radius = round(15/config.resolution(1));
end

end