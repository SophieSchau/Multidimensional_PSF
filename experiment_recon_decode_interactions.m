%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Experiments to test different decoding and reconstruction methods for 
%   simulated vessel encoded angiography with different encoding and
%   sampling schemes. Details presented in abstract 629 at ISMRM 2020
% 
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% 1. Setup for all experiments
phantomname = 'phantom64x64';
N_augmentations = 10;
N_noise_instances = 10;
N_total_spokes = 16;
full_bg_supp = false;
lambda = 5;

load(['data/' phantomname], 'GT');
if full_bg_supp
    GT(:,:,:,:,4) = 0;
end

%% 2a. Setup forward model - Same spokes each encoding + 4x4 encoding

% trajectory
k_same = create_radial_traj(4,N_total_spokes/4,'ga','same'); 

% encoding matrix
H_4 = [-1 -1 -1  1;...
    -1  1  1  1;...
    1 -1  1  1;...
    1  1 -1  1];

% acquisition operator
E_same = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_4,2)],ones(64),[],k_same,H_4);

%% 2b. Setup forward model - Varying spokes each encoding + 4x4 encoding

% trajectory
k_vary = create_radial_traj(4,N_total_spokes/4,'ga','uniform'); 

% encoding matrix = H_4

% acquisition operator
E_vary = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_4,2)],ones(64),[],k_vary,H_4);

%% 2c. Setup forward model - Varying spokes every other encoding + 8x4 encoding

% trajectory
k_alt = create_radial_traj(8,N_total_spokes/8,'ga','same_uniform_alt'); 

% encoding matrix
H_8 = [ -1 -1 -1  1;...
        1  1  1  1;...
        -1  1  1  1;...
        1 -1 -1  1;...
        1 -1  1  1;...
        -1  1 -1  1;...
        1  1 -1  1;...
        -1 -1  1  1];

% acquisition operator
E_alt = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_8,2)],ones(64),[],k_alt,H_8);

%% 3. Simulate k-space data and reconstruct

for aug = 1:N_augmentations
    for vc = 1: size(GT,5)
        rng_seed = 2^aug*3^vc;
        GT_aug(:,:,:,:,vc) = perturb_image(GT(:,:,:,:,vc),1,5,0.5 ,rng_seed);
    end
    for n = 1:N_noise_instances
        savefilename = ['data/recons/' phantomname '_recons_spokes' num2str(N_total_spokes) '_lambda' num2str(lambda) '_pertubation' num2str(aug) '_noise' num2str(n) '_bgsupp' num2str(full_bg_supp) '.mat'];
        if ~exist(savefilename,'file')
            rng(n);
            noise = (10/sqrt(2))*(randn(prod(E_same.dsize),1)+1i*randn(prod(E_same.dsize),1));

            kdata_same = E_same*GT_aug + reshape(noise, E_same.dsize);
            kdata_vary = E_vary*GT_aug + reshape(noise, E_vary.dsize);
            kdata_alt = E_alt*GT_aug + reshape(noise, E_alt.dsize);

            % Decode + recon (only possible for same spokes every encoding)
            kdata_same_decoded = reshape((H_4'*reshape(kdata_same, [],size(H_4,1))')',[size(kdata_same,1),size(kdata_same,2),size(H_4,2)]);
            H1 = eye(size(H_4,1));
            E1 = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H1,2)],ones(64),[],k_same,H1);

            image_same_decode_recon = abs(fista_general(kdata_same_decoded, E1 , 1,...
                lambda, 300, 0.01, 'showProgress', false));

            % Recon + decode
            H1 = eye(size(H_4,1));
            H2 = eye(size(H_8,1));
            E1_same = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H1,2)],ones(64),[],k_same,H1);
            E1_vary = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H1,2)],ones(64),[],k_vary,H1);
            E1_alt = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H2,2)],ones(64),[],k_alt,H2);

            image_same_recon = fista_general(kdata_same, E1_same , 1,...
                lambda, 300, 0.01, 'showProgress', false);
            image_vary_recon = fista_general(kdata_vary, E1_vary , 1,...
                lambda, 300, 0.01, 'showProgress', false);
            image_alt_recon = fista_general(kdata_alt, E1_alt , 1,...
                lambda, 300, 0.01, 'showProgress', false);


            image_same_recon_decode = abs(reshape((H_4'*reshape(image_same_recon, [],size(H_4,1))')',size(GT)));
            image_vary_recon_decode = abs(reshape((H_4'*reshape(image_vary_recon, [],size(H_4,1))')',size(GT)));
            image_alt_recon_decode = abs(reshape((H_8'*reshape(image_alt_recon, [],size(H_8,1))')',size(GT)));


            % joint recon and decoding
            image_same_joint = abs(fista_general(kdata_same, E_same , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_vary_joint = abs(fista_general(kdata_vary, E_vary , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_alt_joint = abs(fista_general(kdata_alt, E_alt , 1,...
                lambda, 300, 0.01, 'showProgress', false));

            save(savefilename,...
                'image_same_decode_recon','image_same_joint','image_same_recon_decode', 'image_vary_joint', 'image_vary_recon_decode', 'image_alt_joint', 'image_alt_recon_decode')
        end
    end

end

%% 4. Compare with ground truth

for aug = 1:N_augmentations
    for vc = 1: size(GT,5)
        rng_seed = 2^aug*3^vc;
        GT_aug(:,:,:,:,vc) = perturb_image(GT(:,:,:,:,vc),1,5,0.5 ,rng_seed);
    end
    for n = 1:N_noise_instances
        loadfilename = ['data/recons/' phantomname '_recons_spokes' num2str(N_total_spokes) '_lambda' num2str(lambda) '_pertubation' num2str(aug) '_noise' num2str(n) '_bgsupp' num2str(full_bg_supp) '.mat'];
        load(loadfilename);
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_same_joint(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_same_joint(aug,n) = gof.rsquare;
        rmserror_same_joint(aug,n) = gof.rmse;
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_same_recon_decode(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_same_recon_decode(aug,n) = gof.rsquare;
        rmserror_same_recon_decode(aug,n) = gof.rmse;
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_same_decode_recon(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_same_decode_recon(aug,n) = gof.rsquare;
        rmserror_same_decode_recon(aug,n) = gof.rmse;
        
        
        
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_vary_joint(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_vary_joint(aug,n) = gof.rsquare;
        rmserror_vary_joint(aug,n) = gof.rmse;
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_vary_recon_decode(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_vary_recon_decode(aug,n) = gof.rsquare;
        rmserror_vary_recon_decode(aug,n) = gof.rmse;

        
        
        
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_alt_joint(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_alt_joint(aug,n) = gof.rsquare;
        rmserror_alt_joint(aug,n) = gof.rmse;
        
        [~, gof] =fit(reshape(GT_aug(:,:,:,:,1:end-1),[],1),reshape(image_alt_recon_decode(:,:,:,:,1:end-1),[],1),'poly1');
        rsquare_alt_recon_decode(aug,n) = gof.rsquare;
        rmserror_alt_recon_decode(aug,n) = gof.rmse;
        
      
    end

end

%% 5. Vizualise result

rsquare_same_decode_recon_mean = mean(rsquare_same_decode_recon(:));
rsquare_same_recon_decode_mean = mean(rsquare_same_recon_decode(:));
rsquare_same_joint_mean = mean(rsquare_same_joint(:));

rsquare_vary_recon_decode_mean = mean(rsquare_vary_recon_decode(:));
rsquare_vary_joint_mean = mean(rsquare_vary_joint(:));

rsquare_alt_recon_decode_mean = mean(rsquare_alt_recon_decode(:));
rsquare_alt_joint_mean = mean(rsquare_alt_joint(:));

barh([rsquare_same_decode_recon_mean 0 0;...
      rsquare_same_recon_decode_mean, rsquare_vary_recon_decode_mean, rsquare_alt_recon_decode_mean;...
      rsquare_same_joint_mean, rsquare_vary_joint_mean, rsquare_alt_joint_mean]')
