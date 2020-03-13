%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simple example for how to use the m_PSF    %
%                                              %
%   Sophie Schauman 2020                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Define k-space sampling pattern
%   We will use golden angle radial sampling within each encoding. 4
%   encodings with 10 spokes each.
%
%   Option 1: Same sampling pattern for each encoding
%   Option 2: Sample pattern rotated by 45 degrees between each encoding

k1 = create_radial_traj(4,10,'ga','same'); % option 1
k2 = create_radial_traj(4,10,'ga','uniform'); % option 2

% visualise sampling patterns
figure(1)
subplot(1,2,1)
for enc = 1:4
scatter(k1(:,1,enc,1), k1(:,1,enc,2),50/enc)
hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
title('same spokes every encoding')

subplot(1,2,2)
for enc = 1:4
scatter(k2(:,1,enc,1), k2(:,1,enc,2),50/enc)
hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
title('different spokes each encoding')
%% 2. Choose encoding matrix
%   We will use a 4x4 hadamard encoding scheme

H = [-1 -1 -1  1;...
    -1  1  1  1;...
    1 -1  1  1;...
    1  1 -1  1];

%% 3. Generate multidimensional PSFs
%   We use a 64x64 field of view, with 10 spokes per encoding and 
%   approximately uniform sampling R = 64*pi*0.5/10 = approx. 10. We assume
%   that every component has equal energy.

pause(1)
figure(2)
mpsf1 = m_psf(k1, H, [64 64]);
title('m-PSF for same spokes every encoding')

pause(1)
figure(3)
mpsf2 = m_psf(k2, H, [64 64]);
title('m-PSF for different spokes each encoding')
 