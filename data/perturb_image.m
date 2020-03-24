function [new_image] = perturb_image(image, displacement_sd, rotation_sd, intensity_sd, rng_seed)
%PERTURB_IMAGE Takes a 2D image and adds a random shift, rotation and
%scaling
%   INPUTS: image           - 2D matrix (image)
%           displacement_sd - standard deviation of gaussian that
%                             determines how much to shift the image
%                             (pixels)
%           rotation_sd     - standard deviation of gaussian that
%                             determines how much to rotate the image
%                             (degrees)
%           intensity_sd    - standard deviation of gaussian that
%                             determines the intensity modulation (mean =
%                             1).
%           rng_seed        - seed to use in random number generator
%                             (optional)
%
%   OUTPUT: new_image       - perturbed image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('rng_seed','var')
    rng(rng_seed)
end

new_image = image;

%change intensity
new_image = new_image*(intensity_sd*randn(1)+1);
%translate x
new_image = circshift(new_image, round(displacement_sd*randn(1)),1);
%translate y
new_image = circshift(new_image, round(displacement_sd*randn(1)),2);
%rotate
new_image = imrotate(new_image, rotation_sd*randn(1),'crop');


end

