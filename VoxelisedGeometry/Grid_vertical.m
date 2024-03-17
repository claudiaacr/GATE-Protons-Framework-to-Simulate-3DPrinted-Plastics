clear
clc

prompt = "Infill Percentage? ";
x = input(prompt);

xz_slab_dim = [100, 10]; %mm
y_slab_dim = 100; %mm
res_xz = 0.1; %mm - enough to capture detail of the pattern but not too large for computation
res_y = 0.2; %mm - should match the size of the shell on that direction

grid_2D = ones(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz);
grid_temp = ones(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz);
SL = ones(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz, y_slab_dim/res_y);
pristine = zeros(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz, y_slab_dim/res_y);

if x == 20
    a = 5.429/2; %mm
elseif x == 50
    a = 1.992/2; %mm
elseif x == 80
    a = 1.132/2; %mm
elseif x == 95
    a = 0.902/2; %mm
end

gap = 1.2; %mm

shell_x = 3*0.2; %mm
shell_y = res_y*3; 
shell_z = 3*0.2; 

N_shell_x = shell_x/res_xz;
N_shell_z = shell_z/res_xz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 2D SLICE MID OF THE SLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTRAL AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the possible centroid point of the polygon

%after the centroid point - start at the central point, add the 2*radius + gap,
%stop at maximum of the slab dimension plus radius (larger distances won't
%overlap with the slab)
centroid_polygon_post_x = xz_slab_dim(1)/2 + (2*a+gap)/2 : 2*a+gap : xz_slab_dim(1) + a;
centroid_polygon_post_z = xz_slab_dim(2)/2 : 2*a+gap : xz_slab_dim(2) + a;

%before the centroid point
centroid_polygon_neg_x =  xz_slab_dim(1)/2 + (2*a+gap)/2 : -(2*a+gap) : -a;
centroid_polygon_neg_z =  xz_slab_dim(2)/2 : - (2*a+gap) : -a;

%concatenate (join or combine) and keep only unique values
centroid_polygon_x = [centroid_polygon_neg_x centroid_polygon_post_x];
centroid_polygon_x = unique(centroid_polygon_x);
centroid_polygon_z = [centroid_polygon_neg_z centroid_polygon_post_z];
centroid_polygon_z = unique(centroid_polygon_z);

%all combinations of centroid points
comb_centroids = combvec(centroid_polygon_x, centroid_polygon_z)';
all_centroids = comb_centroids;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFF AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the possible centroid point of the polygon

%after the centroid point - start at the central point, add the 2*radius + gap,
%stop at maximum of the slag dimension plus radius (larger distances won't
%overlap with the slab)
% centroid_polygon_post_x = xz_slab_dim(1)/2 + (2*a+gap)/2 : 2*a +gap : xz_slab_dim(1) + a;
centroid_polygon_post_z = xz_slab_dim(2)/2 + (2*a+gap)/2 : 2*a +gap : xz_slab_dim(2) + a;
centroid_polygon_post_x = xz_slab_dim(1)/2 : 2*a +gap : xz_slab_dim(1) + a;

%before the centroid point
% centroid_polygon_neg_x =  xz_slab_dim(1)/2 + (2*a+gap)/2 : -(2*a+gap) : -a;
centroid_polygon_neg_z =  xz_slab_dim(2)/2 + (2*a+gap)/2 : -(2*a+gap) : -a;
centroid_polygon_neg_x =  xz_slab_dim(1)/2 : -(2*a+gap) : -a;

%concatenate and keep only unique values
centroid_polygon_x = [centroid_polygon_neg_x centroid_polygon_post_x];
centroid_polygon_x = unique(centroid_polygon_x);
centroid_polygon_z = [centroid_polygon_neg_z centroid_polygon_post_z];
centroid_polygon_z = unique(centroid_polygon_z);

%all combinations of centroid points
comb_centroids = combvec(centroid_polygon_x, centroid_polygon_z)';
all_centroids = [all_centroids; comb_centroids];

%coordinate grid
[colGrid,rowGrid] = meshgrid(0:res_xz:xz_slab_dim(2)-res_xz, 0:res_xz:xz_slab_dim(1)-res_xz);

for i = 1:size(all_centroids,1)
    
    %define the vertices x and y coords
    x_coords(1) = all_centroids(i,1) + a; % right vertice
    x_coords(2) = all_centroids(i,1);     % central vertice (up or down)
    x_coords(3) = all_centroids(i,1) - a; % left vertice   
    x_coords(4) = all_centroids(i,1);     % central vertice (up or down)
     
	z_coords(1) = all_centroids(i,2);     % central vertice
    z_coords(2) = all_centroids(i,2) + a; % upper vertice
    z_coords(3) = all_centroids(i,2);     % central vertice
    z_coords(4) = all_centroids(i,2) - a; % lower vertice
    
    %create the polygin
    pgon = polyshape(z_coords, x_coords);
    % figure(1)
    % plot(pgon); axis square
    
    idx = isinterior(pgon,[colGrid(:),rowGrid(:)]); % true = inside of a polyshape or on a boundary of a polyshap
    idx = reshape(idx,size(grid_2D));

    % Set elements inside the polygon to 0
    grid_2D(idx) = 0;
    figure(1)
    imshow(grid_2D)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD SHELL AT THE END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_2D(1:N_shell_x, :) = 1; % first rows
grid_2D(:, 1:N_shell_z) = 1; % first columns

grid_2D(size(grid_2D,1)-N_shell_x+1 : size(grid_2D,1), :) = 1; % last rows
grid_2D(:, size(grid_2D,2)-N_shell_z+1 : size(grid_2D,2)) = 1; % last columns

figure(2)
imshow(grid_2D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 3D COMPLETE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_3D = repmat(grid_2D, 1, 1, y_slab_dim/res_y);
matrix_3D(:,:,1) = ones(size(grid_2D)); % bottom layer all 1's
matrix_3D(:,:,2) = ones(size(grid_2D)); % bottom layer all 1's
matrix_3D(:,:,3) = ones(size(grid_2D)); % bottom layer all 1's
matrix_3D(:,:, (y_slab_dim/res_y)-2) = ones(size(grid_2D)); % top layer all 1's
matrix_3D(:,:, (y_slab_dim/res_y)-1) = ones(size(grid_2D)); % top layer all 1's
matrix_3D(:,:, y_slab_dim/res_y) = ones(size(grid_2D)); % top layer all 1's

% figure(2)
% subplot(1,4,1)
% imshow(matrix_3D(:,:,3))
% subplot(1,4,2)
% imshow(matrix_3D(:,:,4))
% subplot(1,4,3)
% imshow(matrix_3D(:,:,y_slab_dim/res_y - 3))
% subplot(1,4,4)
% imshow(matrix_3D(:,:,y_slab_dim/res_y - 2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE NII IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
final_matrix_3D = permute(matrix_3D, [1 3 2]);
final_SL = permute(SL, [1 3 2]);
final_pristine = permute(pristine, [1 3 2]);

nii = make_nii(final_matrix_3D, [res_xz res_y res_xz], [0 0 0], 2);
nii_SL = make_nii(final_SL, [res_xz res_y res_xz], [0 0 0], 2);
nii_pristine = make_nii(final_pristine, [res_xz res_y res_xz], [0 0 0], 2);

% save_nii(nii_SL, 'vertical_SL_biggerRes.nii')
% save_nii(nii_pristine, 'vertical_pristine_biggerRes.nii')

if x == 20
    save_nii(nii, 'vertical_20.nii')
elseif x == 50
    save_nii(nii, 'vertical_50.nii')
elseif x == 80
    save_nii(nii, 'vertical_80.nii')
elseif x == 95
    save_nii(nii, 'vertical_95.nii')
end
