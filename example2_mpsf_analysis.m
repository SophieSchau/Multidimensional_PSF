%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simple example for how to analyse the m_PSF     %
%                                                   %
%   Sophie Schauman 2020                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Define k-space sampling pattern
%   We will use golden angle radial sampling within each encoding. 
%   4 encodings with 8 spokes each. We will then reorder the spokes across
%   encodings.
rng(1)
Nencs  = 4;
Nspokes = 8;

k1 = create_radial_traj(Nencs,Nspokes,'ga','uniform');
tmp = reshape(k1,[],Nencs*Nspokes,2);
k2 = reshape(tmp(:,randperm(Nencs*Nspokes),:),size(k1));
k3 = reshape(tmp(:,randperm(Nencs*Nspokes),:),size(k1));
k4 = reshape(tmp(:,randperm(Nencs*Nspokes),:),size(k1));

% visualise sampling patterns
figure(1)
subplot(2,2,1)
for enc = 1:4
    scatter(k1(:,1,enc,1), k1(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,2)
for enc = 1:4
    scatter(k2(:,1,enc,1), k2(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,3)
    for enc = 1:4
    scatter(k3(:,1,enc,1), k3(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,4)
for enc = 1:4
    scatter(k4(:,1,enc,1), k4(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow


set(gcf,'Position', [204 284 867 517])

%% 2. Choose encoding matrix
%   We will use a 4x4 hadamard encoding scheme for option 1 and 2, and a
%   modified 8x4 encoding for option 3.

H = [-1 -1 -1  1;...
    -1  1  1  1;...
    1 -1  1  1;...
    1  1 -1  1];

%% 3. Generate multidimensional PSFs
%   We use a 64x64 field of view, with 16 spokes per encoding and 
%   approximately uniform sampling R = 64*pi*0.5/10 = approx. 6, and for 8
%   spokes per encoding the undersampling factor is therefore approx. 12.

figure(2)
subplot(2,2,1)
mpsf1 = m_psf(k1, H, [64 64]);
imshow_mpsf(mpsf1,[-2,-0])
drawnow

subplot(2,2,2)
mpsf2 = m_psf(k2, H, [64 64]);
imshow_mpsf(mpsf2,[-2,-0])
drawnow

subplot(2,2,3)
mpsf3 = m_psf(k3, H, [64 64]);
imshow_mpsf(mpsf3,[-2,-0])
drawnow

subplot(2,2,4)
mpsf4 = m_psf(k4, H, [64 64]);
imshow_mpsf(mpsf4,[-2,-0])
drawnow

%% 4. Analyse m-PSF's with spectral flatness of off-diagonals
%   Wiener entropy, or spectral flatness, is a measurement of 'peakiness'.
%   (https://en.wikipedia.org/wiki/Spectral_flatness).
%   Reconstruction algorithms based on sparsity work best when aliasing
%   from undersampling is a similar to white noise as possible. Therefore,
%   'non-peaky' aliasing should predict reconstruction quality
%   (hypothesis).
%   
%   Since the diagonal of the m-PSF is independent of spoke ordering we
%   will only measure the wiener entropy of the off-diagonal elements. The
%   m-PSF is also symmetric, so we only need to consider one half of it.

off_diag_comps1 = [];
off_diag_comps2 = [];
off_diag_comps3 = [];
off_diag_comps4 = [];

for n = 1:size(H,2)
    for m = 1:size(H,2)
        if m==n
            % do nothing
        elseif m>n
            % do nothing
        else
            off_diag_comps1 = abs(cat(1,off_diag_comps1,reshape(mpsf1(:,:,m,n),[],1)));
            off_diag_comps2 = abs(cat(1,off_diag_comps2,reshape(mpsf2(:,:,m,n),[],1)));
            off_diag_comps3 = abs(cat(1,off_diag_comps3,reshape(mpsf3(:,:,m,n),[],1)));
            off_diag_comps4 = abs(cat(1,off_diag_comps4,reshape(mpsf4(:,:,m,n),[],1)));
        end
    end
end

spec_flat1 = exp(mean(log(off_diag_comps1)))/mean(off_diag_comps1);
spec_flat2 = exp(mean(log(off_diag_comps2)))/mean(off_diag_comps2);
spec_flat3 = exp(mean(log(off_diag_comps3)))/mean(off_diag_comps3);
spec_flat4 = exp(mean(log(off_diag_comps4)))/mean(off_diag_comps4);

%% 5. Simulate acquisition and reconstruction with the different trajectories




%% 6. Correlate spectral flatness and image reconstruction quality



