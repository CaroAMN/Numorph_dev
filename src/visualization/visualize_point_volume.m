function visualize_point_volume(df, BW, I)
%%% Constants
rf = 1;                     % Rescale factor
res1 = [0.75,0.75,2.5];     % High resolution (centroids)
res2 = [25, 25, 25];        % Low resolution (mask)

%I = imresize3(I,round((res2./res_view).*size(I)));

df2 = df;
%df2(:,1:3) =df2(:,1:3).*(res1./res_view);
df2(:,1:3) =round(df2(:,1:3)*rf).*(res1./res2);

%I = permute(I,[3,1,2]);

% For padding 
%pad_amount = round(size(I)*0.2);
%I = padarray(I,pad_amount,'pre');
%pad_amount(1) = 0;
%I = padarray(I,pad_amount,'post');
%pad2 = pad_amount([2 1 3]);
%df2(:,1:3) = df2(:,1:3) + pad2;

% Create binary mask if provided image
if nargin>2 || isempty(BW)
    [nrows, ncols, nslices] = size(I);
    r1 = round(0.45*[nrows, ncols, nslices]);
    r2 = round(0.55*[nrows, ncols, nslices]);

    chunk = single(I(r1(1):r2(1), r1(2):r2(2), r1(3):r2(3)));
    lowerThresh = (prctile(single(I(:)),2) + prctile(chunk(:),2))/(2*65535);

    BW = imbinarize(I,lowerThresh);
    se = strel('disk',6);
    for i = 1:size(BW,3)
        BW(:,:,i) = imopen(BW(:,:,i),se);
        BW(:,:,i) = imclose(BW(:,:,i),se);
        BW(:,:,i) = imfill(BW(:,:,i),'holes');
    end
    
    % Cleanup volume and binarize
    cc = bwconncomp(BW,6);
    cc_vol = cellfun(@(s) numel(s),cc.PixelIdxList);
    max_cc = find(max(cc_vol));

    BW(:) = 0;
    BW(cc.PixelIdxList{max_cc}) = 1;
end


% Create surface, reduce, and smooth
s = isosurface(BW,0);
%s = reducepatch(s,0.1);
s = smoothpatch(s,1,5,1);

% Sort indexes 
[~,idx] = sort(df2(:,1));
df2 = df2(idx,:);

% Plotting
[scat] = make_plot(df2,s,2,1);


end

function keepAlpha(src,eventData,cEdgeColor)  
    scat = src.MarkerHandle;
    scat.EdgeColorData = EdgeColor;
    scat.FaceColorData = FaceColor;   
end

