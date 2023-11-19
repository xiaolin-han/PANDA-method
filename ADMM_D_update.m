

gamma = 100;
mu_1 = 0.1;
C1 = pinv(R_pan'*R_pan + mu_1 * eye(size(R_pan,2))) * ( R_pan'* Y_L_1 + mu_1* Dksvd_H * A(1:wd,:) );
C2 = pinv(R_pan'*R_pan + mu_1 * eye(size(R_pan,2))) * ( R_pan'* Y_L_2 + mu_1* Dksvd_H * A(wd+1:2*wd,:) );
C3 = pinv(R_pan'*R_pan + mu_1 * eye(size(R_pan,2))) * ( R_pan'* Y_L_3 + mu_1* Dksvd_H * A(2*wd+1:3*wd,:) );
C4 = pinv(R_pan'*R_pan + mu_1 * eye(size(R_pan,2))) * ( R_pan'* Y_L_4 + mu_1* Dksvd_H * A(3*wd+1:4*wd,:) );
A_B = A_update * par.KK ;
D_update = ( gamma * Y_R * A_B' +  mu_1 * ( C1*A(1:wd,:)' + C2*A(wd+1:2*wd,:)' + C3*A(2*wd+1:3*wd,:)' + C4*A(3*wd+1:4*wd,:)' )) * pinv( gamma * A_B * A_B' +...
    mu_1 * ( A(1:wd,:)*A(1:wd,:)' +  A(wd+1:2*wd,:)* A(wd+1:2*wd,:)' + A(2*wd+1:3*wd,:)*A(2*wd+1:3*wd,:)' +  A(3*wd+1:4*wd,:)* A(3*wd+1:4*wd,:)'));