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
rng(10)
Nencs  = 4;
Nspokes = 8;
Npermutations = 10;

k{1} = create_radial_traj(Nencs,Nspokes,'ga','uniform');
tmp = reshape(k{1},[],Nencs*Nspokes,2);
for n = 2:Npermutations
    k{n} = reshape(tmp(:,randperm(Nencs*Nspokes),:),size(k{1}));
end

% visualise the first 4 sampling patterns
figure(1)
subplot(2,2,1)
for enc = 1:4
    scatter(k{1}(:,1,enc,1), k{1}(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,2)
for enc = 1:4
    scatter(k{2}(:,1,enc,1), k{2}(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,3)
for enc = 1:4
    scatter(k{3}(:,1,enc,1), k{3}(:,1,enc,2),50/enc)
    hold on
end
legend('encoding 1', 'encoding 2','encoding 3','encoding 4')
axis equal
drawnow

subplot(2,2,4)
for enc = 1:4
    scatter(k{4}(:,1,enc,1), k{4}(:,1,enc,2),50/enc)
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
%   We use a 64x64 field of view, with 8 spokes per encoding and 
%   approximately uniform sampling R = approx. 12. We will visualize the
%   first 4 m-PSFs

figure(2)
for n = 1:Npermutations
    mpsf{n} = m_psf(k{n}, H, [64 64]);
    if n < 5
        subplot(2,2,n)
        imshow_mpsf(mpsf{n},[-2,-0])
        drawnow
    end        
end


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
off_diag_comps = cell(Npermutations);

for c = 1:Npermutations
    for n = 1:size(H,2)
        for m = 1:size(H,2)
            if m==n
                % do nothing
            elseif m>n
                % do nothing
            else
                off_diag_comps{c} = abs(cat(1,off_diag_comps{c},reshape(mpsf{c}(:,:,m,n),[],1)));
            end
        end
    end
end

for n = 1:Npermutations
    spec_flat{n} = exp(mean(log(off_diag_comps{n})))/mean(off_diag_comps{n});
end


%% 5. Simulate acquisition and reconstruction with the different trajectories
%   We will use a vessel encoded ASL phantom with total backgriund
%   suppression. 10 augmentations (components slightly shifted, rotated and
%   intensity modulated) of a the phantom are simulated. Each dataset is 
%   transformed to k-space using the forward transform where 10 instances 
%   of complex gaussian noise (standard deviation = 10) is added before a 
%   simple compressed sensing reconstruction(with FISTA).

phantomname = 'phantom64x64';
load(['data/', phantomname])
N_augmentations = 10;
N_noise_instances = 10;
lambda = 5;
GT(:,:,:,:,4) = 0;
GT_aug = zeros(size(GT));

% Forward models
for n = 1:Npermutations
    E{n} = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H,2)],ones(64),[],k{n},H);
end

for aug = 1:N_augmentations
    for vc = 1: size(GT,5)
        rng_seed = 2^aug*3^vc;
        GT_aug(:,:,:,:,vc) = perturb_image(GT(:,:,:,:,vc),1,5,0.5 ,rng_seed);
    end
    for n = 1:N_noise_instances
        savefilename = ['data/recons/' phantomname '_recons_spokes' num2str(Nspokes*Nencs) '_lambda' num2str(lambda) '_pertubation' num2str(aug) '_noise' num2str(n) '_bgsupp' num2str(1) '_spoke_order_experiment.mat'];
        if ~exist(savefilename,'file')
            rng(n);
            noise = (10/sqrt(2))*(randn(prod(E{1}.dsize),1)+1i*randn(prod(E{1}.dsize),1));
            for c = 1:Npermutations
                kdata{c} = E{c}*GT_aug + reshape(noise, E{c}.dsize);
                image{c} = abs(fista_general(kdata{c}, E{c} , 1,...
                    lambda, 300, 0.01, 'showProgress', false));
            end

            save(savefilename,...
                'image')
        end
    end
end


%% 6. Correlate spectral flatness and image reconstruction quality

for aug = 1:N_augmentations
    for vc = 1: size(GT,5)
        rng_seed = 2^aug*3^vc;
        GT_aug(:,:,:,:,vc) = perturb_image(GT(:,:,:,:,vc),1,5,0.5 ,rng_seed);
    end
    for n = 1:N_noise_instances
        load(['data/recons/' phantomname '_recons_spokes' num2str(Nspokes*Nencs) '_lambda' num2str(lambda) '_pertubation' num2str(aug) '_noise' num2str(n) '_bgsupp' num2str(1) '_spoke_order_experiment.mat']);

        for c = 1:Npermutations
            [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image{c}(:,:,:,:,1:end-1),[],1),'poly1');
            rsquare{c}(aug,n) = gof.rsquare;
        end  
    end
end
figure(3)
for c = 1:Npermutations
    x(c) = spec_flat{c}(:);
    y(c) = mean(rsquare{c}(:));
    err(c) = std(rsquare{c}(:));
end

errorbar(x,y,err,'x');
xlabel('spectral flatness')
ylabel('R^2')
set(gca,'FontSize',20)
