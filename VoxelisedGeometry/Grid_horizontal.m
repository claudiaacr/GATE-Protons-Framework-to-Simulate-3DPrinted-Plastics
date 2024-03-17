clear
clc

prompt = "Infill Percentage? ";
x = input(prompt);

xy_slab_dim = [100, 100]; %mm
z_slab_dim = 10; %mm
res_xy = 0.1; %mm - enough to capture detail of the pattern but not too large for computation
res_z = 0.2; %mm - should match the size of the shell on that direction

grid_2D = ones(xy_slab_dim(1)/res_xy, xy_slab_dim(2)/res_xy);
grid_temp = ones(xy_slab_dim(1)/res_xy, xy_slab_dim(2)/res_xy);

if x == 20
    a = 5.429/2; %mm
elseif x == 50
    a = 1.992/2; %mm
elseif x == 80
    a = 1.132/2; %mm
end

gap = 0.4; %mm

shell_x = 3*0.2; %mm
shell_y = 3*0.2; %mm
% shell_x = 0.3; %mm
% shell_y = 0.3; %mm
shell_z = res_z*3; 

N_shell_x = shell_x/res_xy;
N_shell_y = shell_y/res_xy;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 2D SLICE MID OF THE SLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTRAL AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the possible centroid point of the polygon

%after the centroid point - start at the central point, add the 2*radius + gap,
%stop at maximum of the slag dimension plus radius (larger distances won't
%overlap with the slab)
centroid_polygon_post_x = xy_slab_dim(1)/2 + (2*a+gap)/2 : 2*a+gap : xy_slab_dim(1) + a;
centroid_polygon_post_y = xy_slab_dim(2)/2 : 2*a+gap : xy_slab_dim(2) + a;

%before the centroid point
centroid_polygon_neg_x =  xy_slab_dim(1)/2 + (2*a+gap)/2 : - (2*a+gap) : -a;
centroid_polygon_neg_y =  xy_slab_dim(2)/2 : - (2*a+gap) : -a;

%concatenate and keep only unique values
centroid_polygon_x = [centroid_polygon_neg_x centroid_polygon_post_x];
centroid_polygon_x = unique(centroid_polygon_x);
centroid_polygon_y = [centroid_polygon_neg_y centroid_polygon_post_y];
centroid_polygon_y = unique(centroid_polygon_y);

%%

%all combinations of centroid points
comb_centroids = combvec(centroid_polygon_x, centroid_polygon_y)';
all_centroids = comb_centroids;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFF AXIS POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the possible centroid point of the polygon

%after the centroid point - start at the central point, add the 2*radius + gap,
%stop at maximum of the slag dimension plus radius (larger distances won't
%overlap with the slab)
centroid_polygon_post_x = xy_slab_dim(1)/2 : 2*a+gap : xy_slab_dim(1) + a;
centroid_polygon_post_y = xy_slab_dim(2)/2 + (2*a+gap)/2 : 2*a+gap : xy_slab_dim(2) + a;

%before the centroid point
centroid_polygon_neg_x =  xy_slab_dim(1)/2 : - (2*a+gap) : -a;
centroid_polygon_neg_y =  xy_slab_dim(2)/2 + (2*a+gap)/2 : - (2*a+gap) : -a;

%concatenate and keep only unique values
centroid_polygon_x = [centroid_polygon_neg_x centroid_polygon_post_x];
centroid_polygon_x = unique(centroid_polygon_x);
centroid_polygon_y = [centroid_polygon_neg_y centroid_polygon_post_y];
centroid_polygon_y = unique(centroid_polygon_y);

%all combinations of centroid points
comb_centroids = combvec(centroid_polygon_x, centroid_polygon_y)';
all_centroids = [all_centroids; comb_centroids];

%coordinate grid
[colGrid,rowGrid] = meshgrid(0:res_xy:xy_slab_dim(2)-res_xy,0:res_xy:xy_slab_dim(1)-res_xy);

for i = 1:size(all_centroids,1)
    
    %define the vertices x and y coords
    x_coords(1) = all_centroids(i,1) + a;
    x_coords(2) = all_centroids(i,1);    
    x_coords(3) = all_centroids(i,1) - a;    
    x_coords(4) = all_centroids(i,1);
     
	y_coords(1) = all_centroids(i,2); 
    y_coords(2) = all_centroids(i,2) + a;
    y_coords(3) = all_centroids(i,2); 
    y_coords(4) = all_centroids(i,2) - a;  
    
    %create the polygin
    pgon = polyshape(y_coords, x_coords);
    % figure(1)
    % plot(pgon); axis square
    
    idx = isinterior(pgon,[colGrid(:),rowGrid(:)]);
    idx = reshape(idx,size(grid_2D));

    % Set elements inside the polygon to 0
    grid_2D(idx) = 0;
    
    figure(1)
    imshow(grid_2D)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD SHELL AT THE END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_2D(1:N_shell_x, :) = 1; % 1st lines
grid_2D(:, 1:N_shell_y) = 1; % 1st columns

grid_2D(size(grid_2D,1)-N_shell_x+1:size(grid_2D,1), :) = 1; % last lines
grid_2D(:, size(grid_2D,2)-N_shell_y+1:size(grid_2D,2)) = 1; % last columns

figure(2)
imshow(grid_2D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE 3D COMPLETE IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_3D = repmat(grid_2D , 1, 1, z_slab_dim/res_z);
matrix_3D(:,:,1) = ones(size(grid_2D));
matrix_3D(:,:,2) = ones(size(grid_2D));
matrix_3D(:,:,3) = ones(size(grid_2D));
matrix_3D(:,:, (z_slab_dim/res_z)-2) = ones(size(grid_2D));
matrix_3D(:,:, (z_slab_dim/res_z)-1) = ones(size(grid_2D));
matrix_3D(:,:, z_slab_dim/res_z) = ones(size(grid_2D));

figure(2)
subplot(1,4,1)
imshow(matrix_3D(:,:,3))
subplot(1,4,2)
imshow(matrix_3D(:,:,4))
subplot(1,4,3)
imshow(matrix_3D(:,:,z_slab_dim/res_z - 3))
subplot(1,4,4)
imshow(matrix_3D(:,:,z_slab_dim/res_z - 2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE NII IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nii = make_nii(matrix_3D, [res_xy res_xy res_z], [0 0 0], 2);

if x   == 20
    save_nii(nii, 'horizontal_Grid_20_gap12.nii')
elseif x == 50
    save_nii(nii, 'horizontal_Grid_50_gap12.nii')
elseif x == 80
    save_nii(nii, 'horizontal_Grid_80_gap12.nii')
end
