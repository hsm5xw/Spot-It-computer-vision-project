function v = new_histvec(image, b)

% function v = histvec(image,mask,b)
%
%     Foundation of Computer Vision;
%     Jason Corso
%
%  For each channel in the image, compute a b-bin histogram (uniformly space
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

c = 1/b;
x = c/2:c:1;

w = ones(b,chan);
for i=1:chan
    
    ii = image(:,:,i);
    ii = reshape(ii, numel(ii), 1);
    w(:,i) = hist(ii,x);
    w(:,i) = w(:,i) / sum(w(:,i));
   
end

v = w(:);