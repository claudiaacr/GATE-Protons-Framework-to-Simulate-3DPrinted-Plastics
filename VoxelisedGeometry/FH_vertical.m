clear
clc

x = input('Infill Percentage? ');

xz_slab_dim = [100, 10]; %mm
y_slab_dim = 100; %mm
res_xz = 0.1; %mm - enough to capture detail of the pattern but not too large for computation
res_y = 0.2; %mm - should match the size of the shell on that direction

FH_2D = ones(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz);
FH_temp = ones(xz_slab_dim(1)/res_xz, xz_slab_dim(2)/res_xz);

if x == 20
    radius = 2.391; %mm
elseif x == 50
    radius = 0.759; %mm
elseif x == 80
    radius = 0.351; %mm
end

% hex_margin = 0.065 cm
gap = 1.2; %mm

shell_x = 3*0.2; %mm
shell_y = res_y*3; 
shell_z = 3*0.2; 

N_shell_x = shell_x/res_xz;
N_shell_z = shell_z/res_xz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTRAL AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_corners = 1.154701*(radius*2); %mm 
hex_edge = 2*(radius/tan(pi/3)); % mm

centroid_hexagon_post_x = xz_slab_dim(2)/2 : dist_corners+hex_edge+(2*gap) : xz_slab_dim(2) + dist_corners/2;
centroid_hexagon_post_z = xz_slab_dim(1)/2 : 2*radius+gap : xz_slab_dim(1) + radius;

%before the centroid point
centroid_hexagon_neg_x =  xz_slab_dim(2)/2 : - (dist_corners+hex_edge+2*gap) : -(dist_corners/2);
centroid_hexagon_neg_z =  xz_slab_dim(1)/2 : - (2*radius+gap): -radius;

%concatenate (join or combine) and keep only unique values
centroid_hexagon_x = [centroid_hexagon_neg_x centroid_hexagon_post_x];
centroid_hexagon_x = unique(centroid_hexagon_x);
centroid_hexagon_z = [centroid_hexagon_neg_z centroid_hexagon_post_z];
centroid_hexagon_z = unique(centroid_hexagon_z);

%all combinations of centroid points
comb_centroids = combvec(centroid_hexagon_x, centroid_hexagon_z)';
all_centroids = comb_centroids;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFF AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centroid_hexagon_post_x = xz_slab_dim(2)/2 + ((dist_corners+hex_edge)/2)+gap : dist_corners+hex_edge+2*gap : xz_slab_dim(2) + dist_corners/2;
centroid_hexagon_post_z = xz_slab_dim(1)/2 + (2*radius+gap)/2 : 2*radius+gap : xz_slab_dim(1) + radius;

%before the centroid point
centroid_hexagon_neg_x =  xz_slab_dim(2)/2 - (((dist_corners+hex_edge)/2)+gap) : -(dist_corners+hex_edge+2*gap) : -(dist_corners/2);
centroid_hexagon_neg_z =  xz_slab_dim(1)/2 - (2*radius+gap)/2 : -(2*radius+gap) : -radius;

%concatenate and keep only unique values
centroid_hexagon_x = [centroid_hexagon_neg_x centroid_hexagon_post_x];
centroid_hexagon_x = unique(centroid_hexagon_x);
centroid_hexagon_z = [centroid_hexagon_neg_z centroid_hexagon_post_z];
centroid_hexagon_z = unique(centroid_hexagon_z);

%all combinations of centroid points
comb_centroids = combvec(centroid_hexagon_x, centroid_hexagon_z)';
all_centroids = [all_centroids; comb_centroids];

%coordinate grid
[colGrid,rowGrid] = meshgrid(0:res_xz:xz_slab_dim(2)-res_xz, 0:res_xz:xz_slab_dim(1)-res_xz);

for i = 1:size(all_centroids,1)
    
    x_coords(1) = all_centroids(i,1) + dist_corners/2;
    x_coords(2) = all_centroids(i,1) + hex_edge/2;    
    x_coords(3) = all_centroids(i,1) - hex_edge/2;    
    x_coords(4) = all_centroids(i,1) - dist_corners/2;
    x_coords(5) = all_centroids(i,1) - hex_edge/2;
    x_coords(6) = all_centroids(i,1) + hex_edge/2;
     
	z_coords(1) = all_centroids(i,2); 
    z_coords(2) = all_centroids(i,2) + radius;
    z_coords(3) = all_centroids(i,2) + radius; 
    z_coords(4) = all_centroids(i,2);
    z_coords(5) = all_centroids(i,2) - radius; 
    z_coords(6) = all_centroids(i,2) - radius;
    
    %create the polygin
    pgon = polyshape(x_coords, z_coords);
    
    idx = isinterior(pgon,[colGrid(:),rowGrid(:)]); % true = inside of a polyshape or on a boundary of a polyshap
    idx = reshape(idx,size(FH_2D));

    % Set elements inside the polygon to 0
    FH_2D(idx) = 0;
    
    figure(1)
    imshow(FH_2D)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD SHELL AT THE END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FH_2D(1:N_shell_x, :) = 1; % first rows
FH_2D(:, 1:N_shell_z) = 1; % first columns

FH_2D(size(FH_2D,1)-N_shell_x+1 : size(FH_2D,1), :) = 1; % last rows
FH_2D(:, size(FH_2D,2)-N_shell_z+1 : size(FH_2D,2)) = 1; % last columns

figure(2)
imshow(FH_2D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 3D COMPLETE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_3D = repmat(FH_2D, 1, 1, y_slab_dim/res_y);
matrix_3D(:,:,1) = ones(size(FH_2D)); % bottom layer all 1's
matrix_3D(:,:,2) = ones(size(FH_2D)); % bottom layer all 1's
matrix_3D(:,:,3) = ones(size(FH_2D)); % bottom layer all 1's
matrix_3D(:,:, (y_slab_dim/res_y)-2) = ones(size(FH_2D)); % top layer all 1's
matrix_3D(:,:, (y_slab_dim/res_y)-1) = ones(size(FH_2D)); % top layer all 1's
matrix_3D(:,:, y_slab_dim/res_y) = ones(size(FH_2D)); % top layer all 1's

figure(3)
subplot(1,4,1)
imshow(matrix_3D(:,:,3))
subplot(1,4,2)
imshow(matrix_3D(:,:,4))
subplot(1,4,3)
imshow(matrix_3D(:,:,y_slab_dim/res_y - 3))
subplot(1,4,4)
imshow(matrix_3D(:,:,y_slab_dim/res_y - 2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE NII IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

final_matrix_3D = permute(matrix_3D, [1 3 2]);
nii = make_nii(final_matrix_3D, [res_xz res_y res_xz], [0 0 0], 2);

% save_nii(nii_SL, 'vertical_SL_biggerRes.nii')
% save_nii(nii_pristine, 'vertical_pristine_biggerRes.nii')

if x == 20
    save_nii(nii, 'vertical_20_FH_gap12.nii')
elseif x == 50
    save_nii(nii, 'vertical_50_FH_gap12.nii')
elseif x == 80
    save_nii(nii, 'vertical_80_FH_gap12.nii')
end
