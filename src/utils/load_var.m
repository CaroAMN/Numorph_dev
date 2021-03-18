function varargout = load_var(config,variable,var_name)
%--------------------------------------------------------------------------
% Load specific variable from output directory.
%--------------------------------------------------------------------------
% Usage:
% varargout = load_var(config,variable,var_name)
%
%--------------------------------------------------------------------------
% Inputs:
% config: Configuration structure containg output directory location.
%
% variable: (string,char) Name of mat file containing variable. Leave
% '.mat' ending off. 
% 
% var_name: (string,char) Name of variable. (default: same as variable)
%--------------------------------------------------------------------------
% Outputs:
% varargout: Collected variable.
%--------------------------------------------------------------------------

if ~ischar(variable) && ~isstring(variable)
   error("2nd input should be the name of the variable to collect")
end

if nargin<3
    var_name = variable;
end

var_path = fullfile(config.output_directory,'variables',strcat(variable,'.mat'));

if isfile(var_path)
    out = load(var_path,var_name);
    f = fields(out);
    varargout = cell(1,length(f));
    for i = 1:length(f)
        varargout{i} = out.(f{i});
    end
else
    error("Could not locate specified variable in the output directory")
end

end