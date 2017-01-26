function [ im ] = bg_preprocess_f( im )


[width, length, ~] = size(im);

green = [im(1,1,1) im(1, 1, 2) im(1, 1, 3)];
white = [0.9 0.9 0.9];
edge1 = [0.4 0.5 0.4];
edge2 = [0.5 0.6 0.5];
edge3 = [0.7 0.8 0.7];

for i = 1:width
    for j = 1:length
        color = [im(i, j, 1) im(i, j, 2) im(i, j, 3)];
        if ( norm(green - color) < 0.3 || norm(white - color) < 0.2 || norm(edge1 - color) < 0.2 || norm(edge2 - color) < 0.2 || norm(edge3 - color) < 0.1)
            im(i, j, :) = [1 1 1];
        end
        
%         color_sum = zeros(1,3);
%         for k = 0:2
%             for l = 0:2
%                 color_sum = color_sum + [im(i+k, j+l, 1) im(i+k, j+l, 2) im(i+k, j+l, 3)];
%             end
%         end
%         color = color_sum/9;
%         
%         if ( norm(green - color) < 0.3 || norm(white - color) < 0.2 || norm(edge1 - color) < 0.2 || norm(edge2 - color) < 0.2 || norm(edge3 - color) < 0.1)
%             for k = 0:2
%                 for l = 0:2
%                     im(i+k, j+l, :) = [1 1 1];
%                 end
%             end
%         end
    end
end

% background_R = imopen(im(:,:,1),strel('disk',15));
% background_G = imopen(im(:,:,2),strel('disk',15));
% background_B = imopen(im(:,:,3),strel('disk',15));
% 
% im(:,:,1) = im(:,:,1) - background_R;
% im(:,:,2) = im(:,:,2) - background_G;
% im(:,:,3) = im(:,:,3) - background_B;
% 
% im(:,:,1) = imadjust(im(:,:,1));
% im(:,:,2) = imadjust(im(:,:,2));
% im(:,:,3) = imadjust(im(:,:,3));

% se = offsetstrel('ball',5,5);

% se = strel('square',3);
% im(:,:,1) = imdilate(im(:,:,1),se);
% im(:,:,2) = imdilate(im(:,:,2),se);
% im(:,:,3) = imdilate(im(:,:,3),se);


% im = imadjust(im);

% h = fspecial('gaussian', [5 5], 2);
% im = imfilter(im, h);

end

