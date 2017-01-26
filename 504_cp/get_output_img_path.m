
function [ img_output_path ] = get_output_img_path( img1_idx, img2_idx, img_type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CARDS  = 1;
IMAGES = 2;

%% setting output img path
img_appendix  = '.png';

if( img_type == CARDS)
    output_base_path = 'cards_results/';
elseif( img_type == IMAGES)
    output_base_path = 'images_results/';
else
    warning('img_type must be either CARDS or IMAGES');
end

img_output_path = [output_base_path, '00', num2str(img1_idx), '_', '00', num2str(img2_idx), img_appendix];  

end
