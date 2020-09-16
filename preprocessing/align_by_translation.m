function coreg_table = align_by_translation(config,path_table,z_displacement)
%--------------------------------------------------------------------------
% Align image channels by translation using MATLAB's imregister
%--------------------------------------------------------------------------

% Unpack important variables from config structure
output_directory = config.output_directory;
lowerThresh = config.lowerThresh;
align_stepsize = config.align_stepsize;
save_aligned_images = config.save_aligned_images;
markers = config.markers;
col = unique(path_table.x);
row = unique(path_table.y);

% Default displacement threshold for comparing to previous translation.
% Multiplied by the average dimension of the image.
distance_threshold = 0.005;

assert(length(row) == 1 & length(col) ==1 , "Image path table should contain only 1 unique tile position")

% Check lower threshold
if any(lowerThresh<1)
    lowerThresh = lowerThresh*65535;
end

% Set precision and number of cc peaks
usfac = 10;
peaks = 5;

% Create matrix with with z displacements for each channel
z_displacement = cat(2,0,z_displacement(:)');
list = zeros(length(unique(path_table.z)),length(markers));
for i = 1:length(markers)  
    list(:,i) = path_table.z(path_table.markers == markers(i)) + z_displacement(i);
end

% Indicate which z positions are not present in all channels
list(list<1) = NaN;
list(list>min(max(list)))=NaN;
list(any(isnan(list), 2),:)=[];

% Creat new table
path_new = cell2table(cell(0,7),'VariableNames',path_table.Properties.VariableNames);

% Add updated z positions for each channel to new table
for i = 1:length(markers)
    table1 = path_table(path_table.markers == markers(i),:);
    table1 = table1(ismember(table1.z,list(:,i)),:);
    table1.z = (1:numel(table1.z))';
    path_new = cat(1,path_new,table1);
end
path_full = path_new;

% Take every n images along z if coreg_stepsize > 1 
min_z = min(path_new.z);
max_z = max(path_new.z);
if align_stepsize > 1
    z = min_z:align_stepsize:max_z;
    path_new = path_new(any(path_new.z == z,2),:);
end

% Create a matrix for recording information. This info will get saved in
% structure array later
order_m = unique(path_new.z)';
order_m = cat(1,order_m,zeros(1+5*(numel(markers)-1),length(order_m)));

% Read each image in stack and measure number of bright pixels
for i = 2:length(markers)
    table2 = path_new(path_new.markers==markers(i),:);
    for j = 1:height(table2)
        mov_img = imread(table2.file{j});
        order_m(i,j) = numel(mov_img(mov_img>lowerThresh(i)))/numel(mov_img);
    end
end

% Sort by average of bright pixels in all channels
sort_row = 1+numel(markers);
order_m(sort_row,:) = mean(order_m(2:length(markers),:),1);

% Adjust distance threshold based on size of image
distance_threshold = distance_threshold * mean(size(mov_img));

% Sort the matrix to register brightest positions first
order_m = sortrows(order_m',sort_row,'descend');
reg_img = cell(1,length(markers)-1);

% Define intensity-based registration settings
% Use optimizer and metric setting
metric = registration.metric.MattesMutualInformation;
optimizer = registration.optimizer.RegularStepGradientDescent;
       
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
       
optimizer.GradientMagnitudeTolerance = 1.00000e-04;
optimizer.MinimumStepLength = 1.00000e-05;
optimizer.MaximumStepLength = 1.00000e-01;
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.5;

% Run channel coregistration
for i = 1:size(order_m,1)
    % Get z position 
    z_idx = order_m(i,1);
   
    % Read reference image
    path_ref = path_new(path_new.z==z_idx & path_new.markers==markers(1),:);
    ref_img = imread(path_ref.file{:});
   
    % For each non-reference channel, perform registration using translation
    for j = 2:numel(markers)
        if order_m(i,j) > 1E-3
        % Index within order matrix
        idx = 2+length(markers)+4*(j-2);

        % Read moving image
        path = path_new(path_new.z==z_idx & path_new.markers==markers(j),:);
        mov_img = imread(path.file{:});
       
       % Do simple crop in case image sizes don't match up
       if any(size(ref_img) ~= size(mov_img))
           mov_img = crop_to_ref(ref_img, mov_img);
       end
       
       % Register using phase correlation. Record 1,cross correlation,error
       [~,~,tform,~,~] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac);
       
       % If tform is empty, use nearest translation
       if isempty(tform)
           m_subset = order_m(order_m(:,idx)~=0,:);
           [~,nearest_idx] = min(abs(m_subset(:,1)-z_idx));
           tform = affine2d([1 0 0; 0 1 0; 0 0 1]);
           tform.T(3) = m_subset(nearest_idx,idx+1);
           tform.T(6) = m_subset(nearest_idx,idx+2);
           order_m(i,idx) = 1; % indicates poor initial phase correlation    
       end
   
       % Check for big translations in phase correlation. It's very important to use
       % a good initial translation before intensity-based registration
       if i>1
           % Take median of translations calculated so far            
           x_med = median(order_m(1:i-1,idx+1));
           y_med = median(order_m(1:i-1,idx+2));
      
           % Calculate distance from median value
           d = ((tform.T(3)-x_med).^2 + (tform.T(6)-y_med).^2).^0.5;
      
           % If translation is too far from the median, take nearest
           % translation
           if d > distance_threshold
               m_subset = order_m(order_m(:,idx)~=0,:);
               [~,nearest_idx] = min(abs(m_subset(:,1)-z_idx));
               tform.T(3) = m_subset(nearest_idx,idx+1);
               tform.T(6) = m_subset(nearest_idx,idx+2);
               order_m(i,idx)=1; %note poor initial phase correlation          
           end
       end
    
       % Use MATLAB's intensity-based registration to refine registration
       % Calculate transform
       tform = imregtform(mov_img,ref_img,'translation',optimizer,metric,...
           'PyramidLevels',2,'InitialTransformation',tform);
   
       % Apply the transform
       reg_img{j-1} = imwarp(mov_img,tform,'OutputView',imref2d(size(ref_img)));

       order_m(i,idx)=order_m(i,idx)+1; %note registration
           
       % Save x,y translations
       order_m(i,idx+1)=tform.T(3);
       order_m(i,idx+2)=tform.T(6);
       
       % Calculate cross correlation
       ref_img(reg_img{j-1} == 0) = 0;
       reg_img{j-1}(ref_img == 0) = 0;
       order_m(i,idx+3) = corr2(reg_img{j-1},ref_img);
       end
    end
    
   % Update status
   fprintf("%s\t %d out %d images completed\n", datetime('now'),i,size(order_m,1));
end

% If stepsize > 1, interpolate values to get missing translations
[~,sort_idx] = sort(order_m(:,1));
order_m = order_m(sort_idx,:);

% Check for outliers and fill in 'no signal' images
for j = 2:length(markers)
    % Index within order matrix
    idx = 2+length(markers)+4*(j-2);

    % Subset translations for this channel
    subset0 = order_m(:,idx:idx+2);
    
    % Remove rows not registered
    nonreg_idx = subset0(:,1)==0;
    reg_idx = subset0(:,1)>0;
    subset1 = subset0(reg_idx,:);
    
    % Detect and replace outlier using moving median sliding window.
    % Use 10 image sliding window
    window = min(10,length(subset1));
    subset1(:,2:3) = filloutliers(subset1(:,2:3),'linear','movmedian',window,1);
    
    % Replace outlier values in order_m
    order_m(reg_idx,idx:idx+2) = subset1;
    
    % Use nearest translation for non-registered images
    nonreg_z = order_m(nonreg_idx,1);
    reg_z = order_m(reg_idx,1);
    for k = 1:length(nonreg_z)
        [~,nearest_idx] = min(reg_z-nonreg_z(k));
        order_m(order_m(:,1) == nonreg_z(k),idx) = 3;
        order_m(order_m(:,1) == nonreg_z(k),idx+1) = order_m(order_m(:,1) == reg_z(nearest_idx),idx+1);
        order_m(order_m(:,1) == nonreg_z(k),idx+2) = order_m(order_m(:,1) == reg_z(nearest_idx),idx+2);
    end
end

% If not every image was registered, interpolate translations for images in
% between 
if align_stepsize >1
   z_range = min_z:max_z;
   z_interp = setdiff(z_range,z);
   interp_m = zeros(length(z_interp),size(order_m,2));
   interp_m(:,1) = z_interp;
   order_m = vertcat(order_m,interp_m);
   for j = 2:length(markers)
       % Index within order matrix
       idx = 2+length(markers)+4*(j-2);
       
       order_m(length(z)+1:end,idx) = 4;
       order_m(length(z)+1:end,idx+1) = interp1(z,order_m(1:length(z),idx+1),z_interp,'linear','extrap');
       order_m(length(z)+1:end,idx+2) = interp1(z,order_m(1:length(z),idx+2),z_interp,'linear','extrap');
   end
   [~,sort_idx] = sort(order_m(:,1));
    order_m = order_m(sort_idx,:);
end

% Save translations
% Remove combined channel column
order_m(:,length(markers)+1) = [];

% Make variable names
vars = cell(1,size(order_m,2));
vars{1} = 'Reference_Z';
for j = 2:length(markers)
   idx = length(markers) + 1 + 4*(j-2);
   c = char(markers(j));
   vars{j} =  strcat('Pixels_',c);
   vars{idx} = strcat('TransformType_',c);
   vars{idx+1} = strcat('X_Shift_',c);
   vars{idx+2} = strcat('Y_Shift_',c);
   vars{idx+3} = strcat('CC_',c);
end

% Convert to table
coreg_table = array2table(order_m,'VariableNames',vars);

% Replace transform type index with specified registration procedure
ttvars = {'Registered','PC_Outlier','Low_Signal','Interpolated'};
idxs = 1 + length(markers) + 4*((2:length(markers))-2);

% Add transform type to coreg_table
tt_table = cell2table(ttvars(order_m(:,idxs))','VariableNames',coreg_table(:,idxs).Properties.VariableNames);
coreg_table(:,idxs) = [];
coreg_table = [coreg_table tt_table];

% Add image filenames
for i = fliplr(1:length(markers))
   file_paths = path_full(path_full.markers == markers(i),1);
   file_paths.Properties.VariableNames = sprintf("file_%d",i);
   coreg_table = cat(2,file_paths,coreg_table);
end

% Save images if needed
if isequal(save_aligned_images,"true")   
   fprintf("%s\t Writing aligned images \n", datetime('now'));
   
    % Create directory to store images
    if exist(fullfile(output_directory,'aligned'),'dir') ~= 7
        mkdir(fullfile(output_directory,'aligned'));
    end
    
    % Get column index for translations for each marker
    t_idx = zeros(1,length(markers));
    for i = 2:length(markers)
        t_idx(i) = find(contains(coreg_table.Properties.VariableNames,'X') & contains(coreg_table.Properties.VariableNames,markers(i)));
    end
    
    for i = 1:height(coreg_table)
        path_sub = coreg_table(coreg_table.Reference_Z == i,:);        
        for j = 1:length(markers)
            mov_img = imread(string(table2cell(coreg_table(i,j))));
            % Apply intensity adjustment
            if isequal(config.shading_correction,'true')
                if i == 1
                   fprintf(strcat(char(datetime('now')),'\t Applying flatfield correction for marker: %s\n'),markers(j));
                end
                   mov_img = apply_intensity_adjustment(mov_img,'flatfield', config.adj_params.flatfield{j},'darkfield',config.adj_params.darkfield{j});
            elseif isequal(config.adjust_ls_width,'true')
                if i == 1
                   fprintf(strcat(char(datetime('now')),'\t Adjusting for light sheet width for marker: %s\n'),markers(j));
                end
                   mov_img = apply_intensity_adjustment(mov_img,'y_adj', config.adj_params.y_adj{j});
            end
            % Apply translations to non-reference image
            if j > 1
                mov_img = imtranslate(mov_img,[path_sub{1,t_idx(j)} path_sub{1,t_idx(j)+1}]);
            end
            % Write aligned images
            img_name = sprintf('%s_%s_C%d_%s_0%d_0%d_aligned.tif',config.sample_name,num2str(coreg_table.Reference_Z(i),'%04.f'),j,markers(j),col,row);
            img_path = fullfile(output_directory,'aligned',img_name);
            imwrite(mov_img,img_path)
        end
    end
end

end