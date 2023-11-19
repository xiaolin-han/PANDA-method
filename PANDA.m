
function [ X_re_3d ] = PANDA(X_ori,Y_R,Y_L,R_pan,psf_space,Dksvd_H,A_0,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues)
%    inputs:
%    X_ori         : ground truth image 
%    Y_R           : multispectral image
%    Y_L           : panchromatic image
%    R_pan         : spectral response function
%    psf_space     : spatial degradation operator
%    Dksvd_H       : initialization of spectral dictionary 
%    A_02          : initialization of sparse representation
%    space_rate    : downsampling rate
%    L             : Image radiometric resolution; 
%    Q_blocks_size : block size of the Q-index locally applied;
%    flag_cut_bounds:cut the boundaries of the viewed Panchromatic image;
%    dim_cut       : define the dimension of the boundary cut;
%    th_values:    : flag. If th_values == 1, apply an hard threshold to the dynamic range.
% 
%   outputs:
%   X_re_3d        : fused multispectral image

    space_rate        =   4;
    X_ori_2d =  change3dto2d(X_ori);
    % Subpixel-shift decomposition
    Y_L_3d = reshape(Y_L',size(X_ori,1),size(X_ori,2));
    Y_L_1 = change3dto2d(Y_L_3d(1:2:end,1:2:end));
    Y_L_2 = change3dto2d(Y_L_3d(2:2:end,1:2:end));
    Y_L_3 = change3dto2d(Y_L_3d(1:2:end,2:2:end));
    Y_L_4 = change3dto2d(Y_L_3d(2:2:end,2:2:end));
    Y_L_hat     = [Y_L_1;Y_L_2;Y_L_3;Y_L_4];
    L_hat       = blkdiag(R_pan,R_pan,R_pan,R_pan);
    D_hat       = blkdiag(Dksvd_H,Dksvd_H,Dksvd_H,Dksvd_H);
    % Sparse Coefficient Optimization
    ADMM_lowrank
    % Spectral dictionary Optimization
    ADMM_D_update
    % Reconstruction
    X_re     = D_update*A_update;
    X_re_3d    = reshape(X_re',size(X_ori,1),size(X_ori,2),4);
    % Quality matric (reference [])
    [Q_avg_Our, SAM_Our, ERGAS_Our, SCC_GT_Our, Q_Our] = indexes_evaluation(X_re_3d,X_ori,space_rate,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
    [Q_Our,SAM_Our,ERGAS_Our,SCC_GT_Our]
end

