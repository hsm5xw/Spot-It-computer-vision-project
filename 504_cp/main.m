% 504 Foundations of Computer Vision
% Challenge Project
% Hong Moon 
% Kevin Ong 

clear all; close all; clc % debug: clear stuffs

% set debug mode (for debugging)
debug_mode   = 0;
CARDS  = 1; IMAGES = 2; % NEVER CHANGE THIS !!!!

% set feature vector length
fvector_size = 128;

%% select images and type ********
img_type = CARDS; % images or cards?
img1_idx = 1;
img2_idx = 3;       

%% setting the img paths
img1_path = get_img_path( img1_idx, img_type);
img2_path = get_img_path( img2_idx, img_type);

% read the images
im1 = double(imread( img1_path ))/255;
im2 = double(imread( img2_path ))/255;


% preprocess to eliminate green background
if( norm([im1(1,1,1) im1(1,1,2) im1(1,1,3)] - [1 1 1]) > 0.2)
    im1 = bg_preprocess_f(im1);
    im2 = bg_preprocess_f(im2);
end

im1_gray = rgb2gray(im1);
im2_gray = rgb2gray(im2);


%% image 1: first compute the superpixels on the image we loaded
[S1,C1] = slic(im1,800);
cmap1 = rand(max(S1(:)),3);
mainfig = figure(99); 
subplot(2,3,1); imagesc(im1); title('input image')
subplot(2,3,2); imagesc(S1);   title('superpixel mask');
colormap(cmap1);

lambda=0.5;
subplot(2,3,3); imagesc( ind2rgb(S1,cmap1)*lambda + im1*(1-lambda));  title('superpixel overlay');

% next compute the feature reduction for the segmentation (histograms)
C1 = reduce(im1,C1,S1);

keyindex = 1; % top-left super-pixel is the key on which to base the foreground extraction
              % used to identify background in dataset (which is white
              % after preprocessing the images


B1 = graphcut(S1,C1,keyindex); 

subplot(2,3,4); imagesc(S1==keyindex); title('keyindex')
subplot(2,3,5); imagesc(B1); title('graphcut result mask')

% Visualize the output
sp4 = subplot(2,3,6);
hw7show(im1,C1,S1,B1,cmap1,sp4)
title('extracted boundaries and fg');
hold off;



%% image 2: first compute the superpixels on the image we loaded
[S2,C2] = slic(im2,800);
cmap2 = rand(max(S2(:)),3);
mainfig = figure(100); 
subplot(2,3,1); imagesc(im2); title('input image')
subplot(2,3,2); imagesc(S2);   title('superpixel mask');
colormap(cmap2);

lambda=0.5;
subplot(2,3,3); imagesc( ind2rgb(S2,cmap2)*lambda + im2*(1-lambda));  title('superpixel overlay');




% next compute the feature reduction for the segmentation (histograms)
C2 = reduce(im2,C2,S2);

keyindex = 1; % top-left super-pixel is the key on which to base the foreground extraction
              % used to identify background in dataset (which is white
              % after preprocessing the images


B2 = graphcut(S2,C2,keyindex); 

subplot(2,3,4); imagesc(S2==keyindex); title('keyindex')
subplot(2,3,5); imagesc(B2); title('graphcut result mask')

% Visualize the output
sp4 = subplot(2,3,6);
hw7show(im2,C2,S2,B2,cmap2,sp4)
title('extracted boundaries and fg');
hold off;
hold off;

%% Image 1 play

B1_flip = ~B1;
CC1 = bwconncomp( B1_flip); % find connected components of 0's in B1ind_x1

%% Image 2 play
B2_flip = ~B2;
CC2 = bwconncomp( B2_flip); % find connected components of 0's in B2


%% Image 1: find k largest connected components, and put boxes around them

list1 = CC1.PixelIdxList; % list of connected components

sizes1 = zeros( 1, length(list1) );

k = min(10, length(list1));

for i=1:length(list1);
    sizes1(i) = length( list1{i} );
end

[sort_sizes1, sort_idx1] = sort( sizes1, 'descend');

select_list1 = list1( sort_idx1(1:k) ); % find the k largest islands where the elemnents are 0 in B1 matrix

hold off;
hold off;
hold off;
hold off;
hold off;

figure(999);
hold off;
% imshow(B1_flip); 
imshow(im1_gray);
hold on;

cell_1 = cell(k,1);

for i=1:k
    [ind_y1, ind_x1] = ind2sub( CC1.ImageSize , select_list1{i} );
    cell_1{i} = [ind_y1, ind_x1];
end

X_IND = 2;
Y_IND = 1;

min_x = zeros(k,1);
min_y = zeros(k,1);

max_x = zeros(k,1);
max_y = zeros(k,1);

box_coord = zeros(k,4); % first 2 stores coord of bottom left corner of box, last 2 stores coord of top right corner of box

hold off;

for i=1:k
     tmp_msg  = ['\leftarrow ' num2str(i)];
     tmp_text = text( double( cell_1{i}(1, X_IND) ), double( cell_1{i}(1,Y_IND) ), tmp_msg);
     tmp_text.FontSize = 15;
     tmp_text.Color = 'red';     

     min_x(i) = min( cell_1{i}(:, X_IND) );
     min_y(i) = min( cell_1{i}(:, Y_IND) );
     
     max_x(i) = max( cell_1{i}(:, X_IND) );
     max_y(i) = max( cell_1{i}(:, Y_IND) );
     
     width   = max_x(i) - min_x(i);
     height  = max_y(i) - min_y(i);

     w_epsilon = floor(width*0.15);
     h_epsilon = floor(height*0.15);
     
     a = max(1,min_x(i)-w_epsilon);  % min_x ++
     b = max(1,min_y(i)-h_epsilon ); % max_x ++
     
     c = min(width + 2*w_epsilon, CC1.ImageSize(1) );  % width ++
     d = min(height + 2*h_epsilon, CC1.ImageSize(2) ); % height ++
     
     box_coord(i, 1) = b;
     box_coord(i, 2) = min(b+d, CC1.ImageSize(2));
     box_coord(i, 3) = a;
     box_coord(i, 4) = min(a+c, CC1.ImageSize(1));
     
     rectangle( 'Position', [a,b,c,d], 'EdgeColor', 'r' ); 
     hold off;
end

hold off;
hold off;

cell_box = cell(k, 1);


for i = 1:k
    
      cell_box{i} = im1( box_coord(i,1):box_coord(i,2), box_coord(i,3):box_coord(i,4), :);

end

hold off;


cell_keypoints = cell(k,3);


%% Image 2: find k largest connected components, and put boxes around them

list2 = CC2.PixelIdxList;

% k2 = length(list2);
k2 = min(10, length(list2));

sizes2 = zeros( 1, length(list2) );

for i=1:length(list2);
    sizes2(i) = length( list2{i} );
end

[sort_sizes2, sort_idx2] = sort( sizes2, 'descend');

select_list2 = list2( sort_idx2(1:k2) ); % find the k largest islands where the elemnents are 0 in B1 matrix

hold off;

figure(2999);
% imshow(B2_flip); 
imshow(im2_gray);
hold on;

cell_2 = cell(k2,1);

for i=1:k2
    [ind_y2, ind_x2] = ind2sub( CC2.ImageSize , select_list2{i} );
    cell_2{i} = [ind_y2, ind_x2];
end

X_IND = 2;
Y_IND = 1;

min_x2 = zeros(k2,1);
min_y2 = zeros(k2,1);

max_x2 = zeros(k2,1);
max_y2 = zeros(k2,1);

box_coord2 = zeros(k2,4); % first 2 stores coord of bottom left corner of box, last 2 stores coord of top right corner of box

hold off;

for i=1:k2
     tmp_msg  = ['\leftarrow ' num2str(i)];
     tmp_text = text( double( cell_2{i}(1, X_IND) ), double( cell_2{i}(1,Y_IND) ), tmp_msg);
     tmp_text.FontSize = 15;
     tmp_text.Color = 'red';     

     min_x2(i) = min( cell_2{i}(:, X_IND) );
     min_y2(i) = min( cell_2{i}(:, Y_IND) );
     
     max_x2(i) = max( cell_2{i}(:, X_IND) );
     max_y2(i) = max( cell_2{i}(:, Y_IND) );
     
     width2   = max_x2(i) - min_x2(i);
     height2  = max_y2(i) - min_y2(i);

     w_epsilon2 = floor(width2*0.15);
     h_epsilon2 = floor(height2*0.15);
     
     a2 = max(1,min_x2(i)-w_epsilon2 );  % min_x ++
     b2 = max(1,min_y2(i)-h_epsilon2 );  % max_x ++
     
     c2 = min(width2 + 2*w_epsilon2, CC2.ImageSize(1) );  % min_y ++
     d2 = min(height2 + 2*h_epsilon2, CC2.ImageSize(2) ); % max_y ++
     
     box_coord2(i, 1) = b2;
     box_coord2(i, 2) = min(b2+d2, CC2.ImageSize(2));
     box_coord2(i, 3) = a2;
     box_coord2(i, 4) = min(a2+c2, CC2.ImageSize(1));
     
     rectangle( 'Position', [a2,b2,c2,d2], 'EdgeColor', 'r' ); 
     hold on;
end

hold off;


cell_box2 = cell(k2, 1);

% for debugging

for i = 1:k2
     cell_box2{i} = im2( box_coord2(i,1):box_coord2(i,2), box_coord2(i,3):box_coord2(i,4), :);

     %figure(2123+i)
     %imshow( cell_box2{i} )  
end

hold off;



%% SURF

num_match_R = zeros(k,k2);
num_match_G = zeros(k,k2);
num_match_B = zeros(k,k2);

metric_R = zeros(k,k2);
metric_G = zeros(k,k2);
metric_B = zeros(k,k2);

for i = 1:k
    % Detect SURF keypoints in R, G, B channels for each object in image 1
    
    target_R = cell_box{i}(:,:,1);
    points_target_R = detectSURFFeatures(target_R);
    [f_target_R,vpts_target_R] = extractFeatures(target_R,points_target_R);
    
    target_G = cell_box{i}(:,:,2);
    points_target_G = detectSURFFeatures(target_G);
    [f_target_G,vpts_target_G] = extractFeatures(target_G,points_target_G);
    
    target_B = cell_box{i}(:,:,3);
    points_target_B = detectSURFFeatures(target_B);
    [f_target_B,vpts_target_B] = extractFeatures(target_B,points_target_B);
    
    
    for j = 1:k2
        % Detect SURF keypoints in R, G, B channels for each object in image 2
        % Match features from each object in image 1 to each object in image 2
        
        match_R = cell_box2{j}(:,:,1);
        points_match_R = detectSURFFeatures(match_R);
        [f_match_R, vpts_match_R] = extractFeatures(match_R, points_match_R);
        [index_pairs_R, matchmetric_R] = matchFeatures(f_target_R, f_match_R);
        
        match_G = cell_box2{j}(:,:,2);
        points_match_G = detectSURFFeatures(match_G);
        [f_match_G, vpts_match_G] = extractFeatures(match_G, points_match_G);
        [index_pairs_G, matchmetric_G] = matchFeatures(f_target_G, f_match_G);
        
        match_B = cell_box2{j}(:,:,3);
        points_match_B = detectSURFFeatures(match_B);
        [f_match_B, vpts_match_B] = extractFeatures(match_B, points_match_B);
        [index_pairs_B, matchmetric_B] = matchFeatures(f_target_B, f_match_B);
        
        
        % Calculate score based on number of matches and metrics from SURF
        % features
        
        num_match_R(i,j) = size(index_pairs_R,1);
        num_match_G(i,j) = size(index_pairs_G,1);
        num_match_B(i,j) = size(index_pairs_B,1);
        
        % If there are no matched features, change metric to 0
        if(numel(matchmetric_R) == 0)
            matchmetric_R = 0;
        end
        
        if(numel(matchmetric_G) == 0)
            matchmetric_G = 0;
        end
        
        if(numel(matchmetric_B) == 0)
            matchmetric_B = 0;
        end
        
        metric_R(i,j) = sum(matchmetric_R);
        metric_G(i,j) = sum(matchmetric_G);
        metric_B(i,j) = sum(matchmetric_B);
        
    end
    
end

num_match = num_match_R + num_match_G + num_match_B;
metric = metric_R + metric_G + metric_B;

c_m = 1;
b_m = 1.2;
exp_metric = b_m*exp(metric/c_m);

[~, sorted_idx] = sort(exp_metric(:), 'descend');
exp_metric_vec = exp_metric(:);
exp_metric_vec(sorted_idx(11:end)) = -1;
exp_metric_vec(exp_metric_vec == b_m) = -1; % make it so that all with same value is -1
exp_metric = reshape(exp_metric_vec, k, k2);


%% find color histogram similarity

sim_hist = zeros(k, k2);
num_bins = 10;

for i=1:k
    
    h1 = new_histvec( cell_box{i}, num_bins);
    
    for j = 1:k2
        h2 = new_histvec( cell_box2{j}, num_bins);
        d =  pdist2( h1', h2');
        
        sim_hist(i,j) = pdist2( h1', h2');
    end
    
end

c_h= 0.05;
b_h = 10;
exp_sim_hist = b_h*exp(-sim_hist/c_h);

[~, sorted_idx] = sort(exp_sim_hist(:), 'descend');
exp_sim_hist_vec = exp_sim_hist(:);
exp_sim_hist_vec(sorted_idx(11:end)) = -1;
exp_sim_hist = reshape(exp_sim_hist_vec, k, k2);


new_metric = exp_metric + exp_sim_hist;

[match_found, match_idx] = max(new_metric(:));

[match_box1, match_box2] = ind2sub([k k2], match_idx);

figure(1234)
subplot(1,2,1)
imshow(cell_box{match_box1});
subplot(1,2,2)
imshow(cell_box2{match_box2});


img_XXX_YYY = zeros(size(B1));

out_x = floor( ( (box_coord(match_box1,1) + box_coord(match_box1,2) )/2 ));
out_y = floor( ( (box_coord(match_box1,3) + box_coord(match_box1,4) )/2 ));


img_XXX_YYY( out_x, out_y) = 1;

%% setting the output img path
img_output_path = get_output_img_path( img1_idx, img2_idx, img_type);


%% store the result image into the output folder
%imwrite(img_XXX_YYY, 'images_results/001_002.png')
imwrite( img_XXX_YYY, img_output_path);




