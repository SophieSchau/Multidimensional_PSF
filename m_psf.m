function mpsf = m_psf(kspace_points,encoding_mat, psf_size, intensity_norm )
%M_PSF creates a map of the multidimensional point spread function
%
% INPUT:
%   kspace_points   -   matrix with k-space sampling locations 
%                       [#samples, 1, #encodings, #dimensions]
%   encoding_mat    -   linear encoding matrix 
%                       [#encodings, #components]
%   psf_size        -   array with matrix size of psf to calculate
%                       [#voxels_in_x, #voxels_in_y, (#voxels_in_z)]
%   intensity_norm  -   array of intensity normalization constants 
%                       [#components] - defaults to ones
%
% OUTPUT:
%   mpsf            -   In 2D, mpsf(:,:,i,j) is the PSF of an impulse in i 
%                       onto component j. In 3D the equivalent is 
%                       mpsf(:,:,:,i,j)
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check inputs 
% kspace_points
if size(kspace_points,4)== 3
    error('3D not implemented yet')
    is3D = true;
elseif size(kspace_points,4)== 2
    is3D = false;
else
    error('kspace_points has the wrong dimensionality')
end

% encoding_mat
if size(size(encoding_mat)) ~= 2
    error('encoding_mat is the wrong size')
elseif size(encoding_mat,1)~=size(kspace_points,3)
    error('encoding_mat or kspace_points is the wrong size')
end
[Nenc, Ncomp] = size(encoding_mat);

%psf_size
if length(psf_size) == 2
    if is3D
        error('kspace_points 2D and psf 3D')
    end
    psf_size = [psf_size, 1];
elseif length(psf_size) == 3
    if ~is3D
        error('kspace_points 3D and psf 2D')
    end
else 
    error('psf_size has the wrong dimensionality')
end

% intensity_norm
if nargin < 4 % set default
    intensity_norm = ones(Ncomp,1);
end

if length(intensity_norm) ~= Ncomp
    error('intensity_norm or encoding_mat is the wrong size')
end

%% create encoding operator
E = xfm_NUFFT_LinearEncoding([psf_size 1, Ncomp,2],ones([psf_size]),[],kspace_points,encoding_mat, 'wi', ones(size(kspace_points,1),1,Nenc));

%% calculate m-psf
for i = 1:Ncomp
    x = zeros([psf_size 1, Ncomp]);
    x(round(end/2),round(end/2),1,i) = intensity_norm(i);    
    if is3D
        mpsf(:,:,:,:,i) = E'*(E*x); 
    else
        mpsf(:,:,:,i) = E'*(E*x); 
    end
end
mpsf = abs(mpsf)/max(abs(mpsf(:)));

%% display result (if 2D)
if ~is3D
    M       =   Ncomp*psf_size(1);
    N       =   Ncomp*psf_size(2);
    im   =   zeros(M,N);

    for c = 1:Ncomp
        for cc = 1:Ncomp
            im((c-1)*psf_size(1)+1:c*psf_size(1),(cc-1)*psf_size(2)+1:cc*psf_size(2)) = mpsf(:,:,c,cc);
        end
    end

    imshow(log10(im),[],'colormap',jet);
    caxis([-2,0])
end

end
