
D         =   L_hat*D_hat;  
Y_2       =   Y_L_hat;

par.lambda     =  0.00001; 
par.mu         =  0.00001; 
eta       =   0.1;

B         =   zeros( size(D,2), size(Y_2,2) );
C         =   zeros( size(D,2), size(Y_2,2) );
mu        =   par.mu;
D2        =   D'*D;
DTY       =   D'*Y_2;
Ek        =   eye(size(D2,2));
V1        =   zeros( size(D,2), size(Y_2,2) );
V2        =   zeros( size(D,2), size(Y_2,2) );
par.psf   =   psf_space;
par.scale =   space_rate;
sz        =   [size(X_ori,1) size(X_ori,2)];
[X par]   =   H_L(X_ori_2d,par,sz);
par.B     =   Set_blur_matrix(par); 
par.KK    =   (par.B)';
[hd wd]   =   size(Dksvd_H);
A_ori_3d  =   reshape(A_0',size(X_ori,1),size(X_ori,2),[]);
A(1:wd,:) =   change3dto2d(A_ori_3d(1:2:end,1:2:end,:));
A(wd+1:2*wd,:)   =  change3dto2d(A_ori_3d(2:2:end,1:2:end,:));
A(2*wd+1:3*wd,:) =  change3dto2d(A_ori_3d(1:2:end,2:2:end,:));
A(3*wd+1:4*wd,:) =  change3dto2d(A_ori_3d(2:2:end,2:2:end,:));

[u s v]      =  svds(A + V2/(2*mu),10);
s_s          =  soft(s, eta/(2*mu));
C            =  u *s_s *v';
B            =  soft(A - V1/(2*mu), par.lambda /(2*mu));
A            =  inv( D2 + 2*mu*Ek ) * (DTY + mu*B+V1/2 + mu*C+V2/2);

A_hat_3d = zeros(size(Y_L_3d,1),size(Y_L_3d,2),size(Dksvd_H,2));
A_hat_3d(1:2:end,1:2:end,:) = reshape(A(1:wd,:)',size(X_ori,1)/2,size(X_ori,2)/2,[]);
A_hat_3d(2:2:end,1:2:end,:) = reshape(A(wd+1:2*wd,:)',size(X_ori,1)/2,size(X_ori,2)/2,[]);
A_hat_3d(1:2:end,2:2:end,:) = reshape(A(2*wd+1:3*wd,:)',size(X_ori,1)/2,size(X_ori,2)/2,[]);
A_hat_3d(2:2:end,2:2:end,:) = reshape(A(3*wd+1:4*wd,:)',size(X_ori,1)/2,size(X_ori,2)/2,[]);
A_update = change3dto2d(A_hat_3d);

