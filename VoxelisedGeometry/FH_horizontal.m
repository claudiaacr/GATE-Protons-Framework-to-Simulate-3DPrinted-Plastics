clear
clc

% prompt = "Infill Percentage? ";
% x = input(prompt);

x = input('Infill Percentage? ');

xy_slab_dim = [100, 100]; %mm
z_slab_dim = 10; %mm
res_xy = 0.1; %mm - enough to capture detail of the pattern but not too large for computation
res_z = 0.2; %mm - should match the size of the shell on that direction

FH_2D = ones(xy_slab_dim(1)/res_xy, xy_slab_dim(2)/res_xy);
FH_temp = ones(xy_slab_dim(1)/res_xy, xy_slab_dim(2)/res_xy);

if x == 20
    radius = 2.391; %mm
elseif x == 50
    radius = 1.992/2; %mm
elseif x == 80
    radius = 1.132/2; %mm
end

% hex_margin = 0.065 cm
gap = 0.4; %mm

shell_x = 3*0.2; %mm
shell_y = 3*0.2; %mm
shell_z = res_z*3; 

N_shell_x = shell_x/res_xy;
N_shell_y = shell_y/res_xy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTRAL AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist_corners = 1.154701*(radius*2); %mm 
hex_edge = 2*(radius/tan(pi/3)); % mm

centroid_hexagon_post_x = xy_slab_dim(1)/2 : dist_corners+hex_edge+(2*gap) : xy_slab_dim(1) + dist_corners/2;
centroid_hexagon_post_y = xy_slab_dim(2)/2 : 2*radius+gap : xy_slab_dim(2) + radius;

%before the centroid point
centroid_hexagon_neg_x =  xy_slab_dim(1)/2 : - (dist_corners+hex_edge+2*gap) : -(dist_corners/2);
centroid_hexagon_neg_y =  xy_slab_dim(2)/2 : - (2*radius+gap): -radius;

%concatenate and keep only unique values
centroid_hexagon_x = [centroid_hexagon_neg_x centroid_hexagon_post_x];
centroid_hexagon_x = unique(centroid_hexagon_x);
centroid_hexagon_y = [centroid_hexagon_neg_y centroid_hexagon_post_y];
centroid_hexagon_y = unique(centroid_hexagon_y);

%all combinations of centroid points
comb_centroids = combvec(centroid_hexagon_x, centroid_hexagon_y)';
all_centroids = comb_centroids;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFF AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centroid_hexagon_post_x = xy_slab_dim(1)/2 + ((dist_corners+hex_edge)/2)+gap : dist_corners+hex_edge+2*gap : xy_slab_dim(1) + dist_corners/2;
centroid_hexagon_post_y = xy_slab_dim(2)/2 + (2*radius+gap)/2 : 2*radius+gap : xy_slab_dim(2) + radius;

%before the centroid point
centroid_hexagon_neg_x =  xy_slab_dim(1)/2 - (((dist_corners+hex_edge)/2)+gap) : -(dist_corners+hex_edge+2*gap) : -(dist_corners/2);
centroid_hexagon_neg_y =  xy_slab_dim(2)/2 - (2*radius+gap)/2 : -(2*radius+gap) : -radius;

%concatenate and keep only unique values
centroid_hexagon_x = [centroid_hexagon_neg_x centroid_hexagon_post_x];
centroid_hexagon_x = unique(centroid_hexagon_x);
centroid_hexagon_y = [centroid_hexagon_neg_y centroid_hexagon_post_y];
centroid_hexagon_y = unique(centroid_hexagon_y);

%all combinations of centroid points
comb_centroids = combvec(centroid_hexagon_x, centroid_hexagon_y)';
all_centroids = [all_centroids; comb_centroids];

%coordinate grid
[colGrid,rowGrid] = meshgrid(0:res_xy:xy_slab_dim(2)-res_xy,0:res_xy:xy_slab_dim(1)-res_xy);

for i = 1:size(all_centroids,1)
    
    x_coords(1) = all_centroids(i,1) + dist_corners/2;
    x_coords(2) = all_centroids(i,1) + hex_edge/2;    
    x_coords(3) = all_centroids(i,1) - hex_edge/2;    
    x_coords(4) = all_centroids(i,1) - dist_corners/2;
    x_coords(5) = all_centroids(i,1) - hex_edge/2;
    x_coords(6) = all_centroids(i,1) + hex_edge/2;
     
	y_coords(1) = all_centroids(i,2); 
    y_coords(2) = all_centroids(i,2) + radius;
    y_coords(3) = all_centroids(i,2) + radius; 
    y_coords(4) = all_centroids(i,2);
    y_coords(5) = all_centroids(i,2) - radius; 
    y_coords(6) = all_centroids(i,2) - radius; 
    
    %create the polygin
    pgon = polyshape(x_coords, y_coords);
    
    idx = isinterior(pgon,[colGrid(:),rowGrid(:)]);
    idx = reshape(idx,size(FH_2D));

    % Set elements inside the polygon to 0
    FH_2D(idx) = 0;
   
    figure(1)
    imshow(FH_2D)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD SHELL AT THE END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FH_2D(1:N_shell_x, :) = 1; % 1st lines
FH_2D(:, 1:N_shell_y) = 1; % 1st columns

FH_2D(size(FH_2D,1)-N_shell_x+1:size(FH_2D,1), :) = 1; % last lines
FH_2D(:, size(FH_2D,2)-N_shell_y+1:size(FH_2D,2)) = 1; % last columns

figure(2)
imshow(FH_2D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 3D COMPLETE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix_3D = repmat(FH_2D , 1, 1, z_slab_dim/res_z);
matrix_3D(:,:,1) = ones(size(FH_2D));
matrix_3D(:,:,2) = ones(size(FH_2D));
matrix_3D(:,:,3) = ones(size(FH_2D));
matrix_3D(:,:, (z_slab_dim/res_z)-2) = ones(size(FH_2D));
matrix_3D(:,:, (z_slab_dim/res_z)-1) = ones(size(FH_2D));
matrix_3D(:,:, z_slab_dim/res_z) = ones(size(FH_2D));
%%
figure(3)
subplot(1,4,1)
imshow(matrix_3D(:,:,3))
subplot(1,4,2)
imshow(matrix_3D(:,:,4))
subplot(1,4,3)
imshow(matrix_3D(:,:,z_slab_dim/res_z - 3))
subplot(1,4,4)
imshow(matrix_3D(:,:,z_slab_dim/res_z - 2))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE NII IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nii = make_nii(matrix_3D, [res_xy res_xy res_z], [0 0 0], 2);

if x   == 20
    save_nii(nii, 'horizontal_FH_20_gap12.nii')
elseif x == 50
    save_nii(nii, 'horizontal_FH_50_gap12.nii')
elseif x == 80
    save_nii(nii, 'horizontal_FH_80_gap12.nii')
end
