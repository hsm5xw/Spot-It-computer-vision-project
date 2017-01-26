function Bmap = segNeighbors(svMap)

%%% function Bmap = segNeighbors(svMap)
%  Foundations of Computer Vision
%
%  Implement the code to compute the adjacency matrix for the superpixel graph
%  captured by svMap
%
%  INPUT:  svMap is an integer image with the value of each pixel being
%           the id of the superpixel with which it is associated
%  OUTPUT: Bmap is a binary adjacency matrix NxN (N being the number of superpixels
%           in svMap).  Bmap has a 1 in cell i,j if superpixel i and j are neighbors.
%           Otherwise, it has a 0.  Superpixels are neighbors if any of their
%           pixels are neighbors.

segmentList = unique(svMap);
segmentNum = length(segmentList);

%%%% IMPLEMENT the calculation of the adjacency

f_Bmap = zeros(segmentNum, segmentNum);

[m,n] = size(svMap);

% within box
for i=2:m-1
    for j=2:n-1
        % above
        if( svMap(i,j) ~= svMap(i-1,j) )
            f_Bmap( svMap(i,j),   svMap(i-1,j) ) = 1;
            f_Bmap( svMap(i-1,j), svMap(i,j)   ) = 1;
        end 

        % below
         if( svMap(i,j) ~= svMap(i+1,j) )
            f_Bmap( svMap(i,j),   svMap(i+1,j) ) = 1;
            f_Bmap( svMap(i+1,j), svMap(i,j)   ) = 1;
         end  

        % left 
        if( svMap(i,j) ~= svMap(i,j-1) )
            f_Bmap( svMap(i,j),   svMap(i,j-1) ) = 1;
            f_Bmap( svMap(i,j-1), svMap(i,j)   ) = 1;
        end     

        % right
        if( svMap(i,j) ~= svMap(i,j+1) )
            f_Bmap( svMap(i,j),   svMap(i,j+1) ) = 1;
            f_Bmap( svMap(i,j+1), svMap(i,j)   ) = 1;
        end
        
        % top-left diagonal
        if( svMap(i,j) ~= svMap(i-1,j-1) )
            f_Bmap( svMap(i,j),   svMap(i-1,j-1) ) = 1;
            f_Bmap( svMap(i-1,j-1), svMap(i,j)   ) = 1;
        end 
        
        % top-right diagonal
        if( svMap(i,j) ~= svMap(i-1,j+1) )
            f_Bmap( svMap(i,j),   svMap(i-1,j+1) ) = 1;
            f_Bmap( svMap(i-1,j+1), svMap(i,j)   ) = 1;
        end        
        
        % bottom-left diagonal
        if( svMap(i,j) ~= svMap(i+1,j-1) )
            f_Bmap( svMap(i,j),   svMap(i+1,j-1) ) = 1;
            f_Bmap( svMap(i+1,j-1), svMap(i,j)   ) = 1;
        end        
        
        % bottom-right diagonal
         if( svMap(i,j) ~= svMap(i+1,j+1) )
            f_Bmap( svMap(i,j),   svMap(i+1,j+1) ) = 1;
            f_Bmap( svMap(i+1,j+1), svMap(i,j)   ) = 1;
        end       
        
    end
end

% create a spare matrix
Bmap = sparse( f_Bmap);
