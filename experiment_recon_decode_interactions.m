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
full_bg_supp = true;
lambda = 1;

load(['data/' phantomname], 'GT');
if full_bg_supp
    GT(:,:,:,:,4) = 0;
end

%% 2a. Setup forward model - Same spokes each encoding + 4x4 encoding

% trajectory
k_a = create_radial_traj(4,N_total_spokes/4,'ga','same'); 

% encoding matrix
H_a = [-1 -1 -1  1;...
    -1  1  1  1;...
    1 -1  1  1;...
    1  1 -1  1];

% acquisition operator
E_a = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_a,2)],ones(64),[],k_a,H_a);

%% 2b. Setup forward model - Varying spokes each encoding + 4x4 encoding

% trajectory
k_b = create_radial_traj(4,N_total_spokes/4,'ga','uniform'); 

% encoding matrix
H_b = [-1 -1 -1  1;...
    -1  1  1  1;...
    1 -1  1  1;...
    1  1 -1  1];

% acquisition operator
E_b = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_b,2)],ones(64),[],k_b,H_b);

%% 2c. Setup forward model - Varying spokes every other encoding + 8x4 encoding

% trajectory
k_c = create_radial_traj(8,N_total_spokes/8,'ga','same_uniform_alt'); 

% encoding matrix
H_c = [ -1 -1 -1  1;...
        1  1  1  1;...
        -1  1  1  1;...
        1 -1 -1  1;...
        1 -1  1  1;...
        -1  1 -1  1;...
        1  1 -1  1;...
        -1 -1  1  1];

% acquisition operator
E_c = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H_c,2)],ones(64),[],k_c,H_c);

%% 3. Simulate k-space data and reconstruct

for aug = 1:N_augmentations
    for vc = 1: size(GT,5)
        GT(:,:,:,:,vc) = perturb_image(GT(:,:,:,:,vc),1,5,0.5 ,aug+vc);
    end
    for n = 1:N_noise_instances
        savefilename = ['data/recons/' phantomname '_recons_spokes' num2str(N_total_spokes) '_lambda' num2str(lambda) '_pertubation' num2str(aug) '_noise' num2str(n)];
        if ~exist(savefilename,'file')
            rng(n);
            noise = (10/sqrt(2))*(randn(prod(E_a.dsize),1)+1i*randn(prod(E_a.dsize),1));

            kdata_a = E_a*GT + reshape(noise, E_a.dsize);
            kdata_b = E_b*GT + reshape(noise, E_b.dsize);
            kdata_c = E_c*GT + reshape(noise, E_c.dsize);

            % Decode + recon (only possible for same spokes every encoding)
            kdata_a_decoded = reshape((H_a'*reshape(kdata_a, [],size(H_a,1))')',[size(kdata_a,1),size(kdata_a,2),size(H_a,2)]);
            H1 = eye(size(H_a,1));
            E1 = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H1,2)],ones(64),[],k_a,H1);

            image_a_decode_recon = abs(fista_general(kdata_a_decoded, E1 , 1,...
                lambda, 300, 0.01, 'showProgress', false));

            % Recon + decode
            H1 = eye(size(H_a,2));
            H2 = eye(size(H_c,1));
            E1_ab = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H1,2)],ones(64),[],k_a,H1);
            E1_c = xfm_NUFFT_LinearEncoding([64 64 1 1 size(H2,2)],ones(64),[],k_c,H2);

            image_a_recon = abs(fista_general(kdata_a, E1_ab , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_b_recon = abs(fista_general(kdata_b, E1_ab , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_c_recon = abs(fista_general(kdata_c, E1_c , 1,...
                lambda, 300, 0.01, 'showProgress', false));


            image_a_recon_decode = reshape((H_a'*reshape(image_a_recon, [],size(H_a,1))')',[size(image_a_recon,1),size(image_a_recon,2),size(image_a_recon,3),size(image_a_recon,4),size(H_a,2)]);
            image_b_recon_decode = reshape((H_b'*reshape(image_b_recon, [],size(H_b,1))')',[size(image_b_recon,1),size(image_b_recon,2),size(image_b_recon,3),size(image_b_recon,4),size(H_b,2)]);
            image_c_recon_decode = reshape((H_c'*reshape(image_c_recon, [],size(H_c,1))')',[size(image_c_recon,1),size(image_c_recon,2),size(image_c_recon,3),size(image_c_recon,4),size(H_c,2)]);


            % joint recon and decoding
            image_a_joint = abs(fista_general(kdata_a, E_a , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_b_joint = abs(fista_general(kdata_b, E_b , 1,...
                lambda, 300, 0.01, 'showProgress', false));
            image_c_joint = abs(fista_general(kdata_c, E_c , 1,...
                lambda, 300, 0.01, 'showProgress', false));

            save(savefilename,...
                'image_a_decode_recon','image_a_joint','image_a_recon_decode', 'image_b_joint', 'image_b_recon_decode', 'image_c_joint', 'image_c_recon_decode')
        end
    end

end

