function v = histvec(image,mask,b)

% function v = histvec(image,mask,b)
%
%     Foundation of Computer Vision;
%     Jason Corso
%
%  For each channel in the image, compute b-bin histogram (uniformly space
%  bins in the range 0:1) of the pixels in image where the mask is true. 
%  Then, concatenate the vectors together into one column vector (first
%  channel at top).
%
%  mask is a matrix of booleans the same size as image.
% 
%  normalize the histogram of each channel so that it sums to 1.
%
%  You CAN use the hist function.
%  You MAY loop over the channels.

chan = size(image,3);

c = 1/b;       % bin offset
x = c/2:c:1;   % bin centers


%%%%% IMPLEMENT below this line into a 3*b by 1 vector  v
%%  3*b because we have a color image and you have a separate 
%%  histogram per color channel

%edges = linspace(0,1,b+1);

im1 = image(:,:,1);
im2 = image(:,:,2);
im3 = image(:,:,3);

% mask -> (a matrix with either 0 or 1)
centers = x;
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end) ];

% find indices for region of interest (where mask is true)
idx = find(mask ==1); 

im1_roi = im1( idx);
im2_roi = im2( idx);
im3_roi = im3( idx);

% Get 3 histograms, one for each channel
[N1, edges1] = histcounts( im1_roi, edges, 'Normalization','probability' );
[N2, edges2] = histcounts( im2_roi, edges, 'Normalization','probability' );
[N3, edges3] = histcounts( im3_roi, edges, 'Normalization','probability' );

v = [N1';N2';N3'];







