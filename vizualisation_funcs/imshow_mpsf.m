function [] = imshow_mpsf(mpsf,intensity_window)
%IMSHOW_MPFS Displays a multidimensional PSF
%   This function displays a grid of images that describes a
%   multidimensional point spread function (m-PSF). The images are log10
%   transformed, and the user can set a window. Only works for 2D m-PSFs so
%   far.
%
%   INPUTS:
%           mpsf - 4D multidimensional point spread function 
%                  [Nx, Ny, Ncomp, Ncomp]
%
%           intensity_window  - works in same way as imshow(), array of two
%                               numbers for min and max. Empty brackets
%                               scales the image to [min(log10(mpsf(:))),
%                               max(log10(mpsf(:)))]. Default is empty
%                               brackets.
%
%   SOPHIE SCHAUMAN 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        intensity_window = [];
    end

    M       =   size(mpsf,4)*size(mpsf,1);
    N       =   size(mpsf,4)*size(mpsf,2);
    im   =   zeros(M,N);

    for c = 1:size(mpsf,4)
        for cc = 1:size(mpsf,4)
            im((c-1)*size(mpsf,1)+1:c*size(mpsf,1),(cc-1)*size(mpsf,2)+1:cc*size(mpsf,2)) = mpsf(:,:,c,cc);
        end
    end

    imshow(log10(abs(im)),intensity_window,'colormap',parula);

end

