function [B] = graphcut(segmentimage,segments,keyindex)

% function [B] = graphcut(segmentimage,segments,keyindex
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Function to take a superpixel set and a keyindex and convert to a 
%  foreground/background segmentation.
%
% keyindex is the index to the superpixel we wish to use as foreground and
% find its relevant neighbors that would be in the same macro-segment
%
% Similarity is computed based on segments(i).fv (which is a color histogram)
%  and spatial proximity.
%
% segmentimage and segments are returned by the superpixel function
%  segmentimage is called S and segments is called Copt  
%
% OUTPUT:  B is a binary image with 1's for those  connected to the
%          source node and hence in the same segment as the keyindex.
%          B has 0's for those nodes connected to the sink.

% compute basic adjacency information of superpixels
%%%%  Note that segNeighbors is code you need to implement
adjacency = segNeighbors(segmentimage);
%debug
%figure; imagesc(adjacency); title('adjacency');

% normalization for distance calculation based on the image size
% for points (x1,y1) and (x2,y2), distance is
% exp(-||(x1,y1)-(x2,y2)||^2/dnorm)
dnorm = 2*prod(size(segmentimage)/2)^2;

k = length(segments);

capacity = zeros(k+2,k+2);  % initialize the zero-valued capacity matrix
source = k+1;  % set the index of the source node
sink = k+2;    % set the index of the sink node

% this is a single planar graph with an extra source and sink
%
% Capacity of a present edge in the graph (adjacency) is to be defined as the product of
%  1:  the histogram similarity between the two color histogram feature vectors.
%      use the provided histintersect function below to compute this similarity 
%  2:  the spatial proximity between the two superpixels connected by the edge.
%      use exp(-D(a,b)/dnorm) where D is the euclidean distance between superpixels a and b,
%      dnorm is given above.
%
% source gets connected to every node except sink
%  capacity is with respect to the keyindex superpixel
% sink gets connected to every node except source; 
%  capacity is opposite that of the corresponding source-connection (from each superpixel)
%  in our case, the max capacity on an edge is 3; so, 3 minus corresponding capacity
% 
% superpixels get connected to each other based on computed adjacency matrix
%  capacity defined as above. EXCEPT THAT YOU ALSO NEED TO MULTIPLY BY A SCALAR 0.25 for
%  adjacent superpixels.


%%% IMPLEMENT CODE TO fill in the capacity values using the description above.

f_adjacency = full(adjacency); % k x k

for i=1:k
    for j=1:k   
        % 1. histogram similarity    
        c1 = histintersect( segments(i).fv, segments(j).fv) ;

        % 2. spatial proximity
    
        % Euclidean distance between two superpixels
        D =  sqrt( (segments(i).x-segments(j).x)^2 + (segments(i).y-segments(j).y)^2 );
        c2 = exp( -2*D/dnorm ); % !!!!!! temporary change !!!!
 
        % if two super-pixels are adjacent
        if( f_adjacency(i,j) )
            capacity(i,j) = 0.25*c1*c2;
        else   
            capacity(i,j) = 0;
        end
    end
end

% source (W_s,i)
for i=1:k
   % 1. histogram similarity    
   c1 = histintersect( segments(i).fv, segments(keyindex).fv) ;
 
   % 2. spatial proximity  
   % Euclidean distance between two superpixels
   D =  sqrt( (segments(i).x-segments(keyindex).x)^2 + (segments(i).y-segments(keyindex).y)^2 );
   c2 = exp( -D/dnorm );
   
   capacity(source,i) = c1*c2; % W_s,i
   capacity(i,source) = 0; 
end

% sink (W_i,t)
for i=1:k
  capacity(i,sink) = 3 - capacity(source,i); % W_i,t = 3 - W_s,i
  capacity(sink,i) = 0;  
end

capacity(source:sink,source:sink) = 0;


%debug
%figure; imagesc(capacity); title('capacity');

% compute the cut (this code is provided to you)
[~,current_flow] = ff_max_flow(source,sink,capacity,k+2);

% extract the two-class segmentation.
%  the cut will separate all nodes into those connected to the
%   source and those connected to the sink.
%  The current_flow matrix contains the necessary information about
%   the max-flow through the graph.
%
% Populate the binary matrix B with 1's for those nodes that are connected
%  to the source (and hence are in the same big segment as our keyindex) in the
%  residual graph.
% 
% You need to compute the set of reachable nodes from the source.  Recall, from
%  lecture that these are any nodes that can be reached from any path from the
%  source in the graph with residual capacity (original capacity - current flow) 
%  being positive.

%%%  IMPLEMENT code to read the cut into B

residual_capacity = capacity - current_flow;

% find set of reachable nodes (from source)

    % implement BFS (breadth-first search)
    
    reachable_nodes = [];
    num_nodes = length( f_adjacency);

    % mark all vertices as not visited
    visited = zeros( k+2,1);

    % set the start node as visited;
    visited(source) = 1;

    q = [];            % queue
    q(end+1) = source; % push start node

    while( length(q) ~= 0) 

        % Dequeue a vertex from queue
        v = q(1);
        q = q(2:end);

        % add to reachable nodes
        reachable_nodes(end+1) = v;

        % inspect all adjacent vertices
        for i=1:k+2                   

            % adjacent vertices that are not visited
            if ( ~visited(i) && (residual_capacity(v,i)>0) )
                visited(i) = 1;
                q(end+1) = i; % enqueue the visited vertex
            end    
        end         
    end   

[m,n] = size(segmentimage);

B = zeros(m,n);

    for r=1:length(reachable_nodes)
        for i=1:m
            for j=1:n
                 % super pixel associated with reachable node
                r_superpixel = reachable_nodes(r);

                if( segmentimage(i,j) == r_superpixel)
                    B(i,j) = 1;
                end    
            end
        end
    end

    
end

function c = histintersect(a,b)
    c = sum(min(a,b));
end
