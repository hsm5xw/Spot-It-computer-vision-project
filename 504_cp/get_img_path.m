function [ img_path ] = get_img_path( img_idx, img_type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CARDS  = 1;
IMAGES = 2;

%% setting the img paths
img_base_path = '';
img_appendix  = '.png';

if( img_type == CARDS)
    img_base_path = 'cp_data/cards/00';
elseif( img_type == IMAGES)
    img_base_path = 'cp_data/images/00';
else
    warning('img_type must be either CARDS or IMAGES');
end

img_path = [img_base_path, num2str(img_idx), img_appendix];

end

