clc
clear

slab = load_untouch_nii('0_NOVO/images Grid/vertical_20_gap12.nii');
slab_im = slab.img;
n_pla = 100*length(slab_im(slab_im==1))/numel(slab_im);
n_air = 100*length(slab_im(slab_im==0))/numel(slab_im);