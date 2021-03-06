% use this to find weights
% ...

function[weights] = get_weights(mask_xell, PixelIdxList, Xell_nearest)

PixelBoundaryList = find_boundary(mask_xell, PixelIdxList);
weights = get_boundary_size(PixelBoundaryList, Xell_nearest);

end

% find boundary weights:
% for each cell returns a vector of length (number of nearest neighbours)
% where each entry is the number of shared pixels along the boundary
% between the two cells
% kieranrcampbell@gmail.com
function[RelBoundarySize] = get_boundary_size(PixelBoundaryList, Xell_nearest)
RelBoundarySize = cell(size(PixelBoundaryList));

for i = 1:length(PixelBoundaryList)
    nearest_cells = Xell_nearest{i}(:,1);
    for j = nearest_cells'

        cell_intersect = length(intersect(PixelBoundaryList{i}, PixelBoundaryList{j}));
        
        if(cell_intersect == 0)
            error('Nearest neighbour finding gone wrong')
        end
        RelBoundarySize{i} = [RelBoundarySize{i} cell_intersect];
    end
end

end


% given a mask and the PixelIdxList returned by
% bwconncomp, this returns the positions of the pixels
% on the boundary (the '1's) surrounding a given cell
% output used in get_boundary_size
% kieranrcampbell@gmail.com

function[PixelBoundaryList] = find_boundary(mask_xell, PixelIdxList, connectivity)

if nargin < 3
    connectivity = 8; % default connectivity in bwconncomp
end

PixelBoundaryList = cell(size(PixelIdxList));

for i = 1:length(PixelIdxList)
    mask_temp = zeros(size(mask_xell));
    pixels = PixelIdxList{i};
    mask_temp(pixels) = 1;
    %imshow(mask_temp, [])
    PixelBoundaryList{i} = find_single_boundary(mask_temp, connectivity);
end

end


function[pixelList] = find_single_boundary(BW, connectivity)

if nargin < 2
    connectivity = 8;
end

B = bwboundaries(BW);

B = cell2mat(B);

B = B(1:end-1,:);

boundary_size = size(B,1);

d = [];
if connectivity == 8
    d = [ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1  ]; 
elseif connectivity == 4
    d = [ 1 0 ; -1 0 ; 0 -1 ; 0 1 ];
end
    
drep = [];

for i = d'
    %size(i)
    drep = vertcat(drep, repmat(i', [boundary_size,1]));
end


boundary_all = drep+repmat(B, [connectivity,1]);

remove_rows = any(boundary_all < 1, 2);
[x_size, y_size] = size(BW);

remove_rows = or(remove_rows, any(boundary_all(:, 1) > x_size, 2));
remove_rows = or(remove_rows, any(boundary_all(:, 2) > y_size, 2));

boundary_all(remove_rows,:) = [];

img = BW;

for i = boundary_all'
    img(i(1), i(2)) = 1;
end

%disp(img)

boundary = xor(img, BW);

%disp(boundary)

pixelList = find(boundary);

end


