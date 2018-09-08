function  [z1, z2] = TestData

%Test
%Number of Observations
n = 200;


% one dimensional latent predictor variable
a = [1 1.3 1.8];
b = [1 1.1 1.7];
eta = randn(n,1);
z1 = randn(n,1);                                %norm
% z1 = rand(n,1);                                 %unif
% z2 = sin(2.*pi.*z1) + 0.5.*eta;                 %sin2pi
% z2 = 0.2.*exp(5.*z1)-25.*z1.^3 + 0.75.*eta;     %exp
z2 = (1./(5.*z1+1)) + sin(5*z1) + 0.75*eta;     %sin
X = zeros(3,n);
Y = zeros(3,n);
for i = 1:3
    X(i,:) = a(i).*z1 + 0.15.*randn(n,1);
    Y(i,:) = b(i).*z2 + 0.25*randn(n,1);
end
Matrix = [X;Y];
Z = [z1, z2]; 
dlmwrite('Data_1d_unif_N200_sin',Matrix)
dlmwrite('Data_1d_unif_N200_sin_Z',Z)



% % tow dimensional latent predictor variable
% a1 = [1 1.3 1.8];
% a2 = [1 1.5 2 1.3];
% b = [1 1.1 1.7];
% eta = randn(n,1);
% z1_1 = rand(n,1);   %unif
% z1_2 = rand(n,1);   %unif
% % z1_1 = randn(n,1);    %norm
% % z1_2 = randn(n,1);    %norm
% z2 = z1_1.*sin(z1_1.^2) - z1_2.*sin(z1_2.^2) + 0.15.*eta;   %sin2d
% % z2 = 4./(1+4*(z1_1-0.5).^2 + 4*(z1_2-0.5).^2) + 0.15.*eta;    %quad
% X_1 = zeros(3,n);
% X_2 = zeros(4,n);
% Y = zeros(3,n);
% for i = 1:3
%     X_1(i,:) = a1(i).*z1_1 + 0.15.*randn(n,1);
%     Y(i,:) = b(i).*z2 + 0.2.*randn(n,1);
% end
% for i = 1:4
%     X_2(i,:) = a2(i).*z1_2 + 0.15.*randn(n,1);
% end
% Matrix = [X_1;X_2;Y];
% Z = [z1_1; z1_2; z2]; 
% dlmwrite('Data_2d_unif_N200_sin2d',Matrix)
% dlmwrite('Data_2d_unif_N200_sin2d_Z',Z)
% 
% 
% 
% % three dimensional latent predictor variable
% a1 = [1 1.3 1.8];
% a2 = [1 1.5 2 1.3];
% a3 = [1 0.8 1.4 0.9 1.7];
% b = [1 1.1 1.7];
% eta = randn(n,1);
% z1_1 = rand(n,1);   %unif
% z1_2 = rand(n,1);   %unif
% z1_3 = rand(n,1);   %unif
% z2 = z1_1.*sin(z1_1.^2) - z1_2.*sin(z1_2.^2)...
%     + z1_3.*sin(z1_3.^2) + 0.15.*eta;                   %sin3d
% X_1 = zeros(3,n);
% X_2 = zeros(4,n);
% X_3 = zeros(5,n);
% Y = zeros(3,n);
% for i = 1:3
%     X_1(i,:) = a1(i).*z1_1 + 0.15.*randn(n,1);
%     Y(i,:) = b(i).*z2 + 0.2.*randn(n,1);
% end
% for i = 1:4
%     X_2(i,:) = a2(i).*z1_2 + 0.15.*randn(n,1);
% end
% for i = 1:5
%     X_3(i,:) = a3(i).*z1_3 + 0.15.*randn(n,1);
% end
% Matrix = [X_1;X_2;X_3;Y];
% Z = [z1_1; z1_2; z1_3; z2]; 
% dlmwrite('Data_3d_unif_N200_sin3d',Matrix)
% dlmwrite('Data_3d_unif_N200_sin3d_Z',Z)
