function [k] = create_radial_traj(Nenc,Nspokes_per_enc,method_within_enc,...
    method_across_enc, random_seed)
%CREATE_RADIAL_TRAJ creates a matrix with radial k-space sampling locations 
%   
%   INPUTS:
%       Nenc                 -   Number of encodings (int)
%       Nspokes_per_enc      -   Number of spokes per encoding (int)
%       method_within_enc    -   Sampling method within encoding (str)
%                                  'uniform' (radially uniform sampling)
%                                  'GA' (golden angle sampling)
%       method_across_enc    -   Sampling method across encodings (str)
%                                  'same' (same spokes each encoding)
%                                  'uniform' (uniform sampling across encs)
%                                  'GA' (golden angle rotaion between encs)
%                                  'same_uniform_alt' (alternating between
%                                                     same and uniform 
%                                                     between encs)
%                                  'same_GA_alt' (alternating between
%                                                same and golden angle 
%                                                between encs)
%                                  'random' (random rotation across encs)
%       random_seed (optional) -  Seed for RNG (int)
%
%   OUTPUT:
%       k                      - Matrix with k-space sampling locations. 
%                                Size: [#samples, 1, #encodings, 2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = (1+sqrt(5))/2; % Golden ratio

% setup first encoding angles
if strcmpi(method_within_enc, 'ga')
    angles_enc1 = (180/g)*(0:Nspokes_per_enc-1);
elseif strcmpi(method_within_enc, 'uniform')
    angles_enc1 = (180/Nspokes_per_enc)*(0:Nspokes_per_enc-1);
else
    error(['method_within_enc named ' upper(method_within_enc) ' is not implemented'])
end
    

k = zeros(Nspokes_per_enc*256, 1,Nenc,2);
% figure out how much to rotate between encodings
switch lower(method_across_enc)
    case 'same'
        incr = zeros(Nenc-1,1);
    case 'uniform'
        incr = 180/Nenc*ones(Nenc-1,1);
    case 'ga'
        incr = 180/g*ones(Nenc-1,1);
    case 'same_uniform_alt'
        incr = 360/Nenc*ones(Nenc-1,1);
        incr(1:2:end) = 0;
    case 'same_ga_alt'
        incr = 180/g*ones(Nenc-1,1);
        incr(1:2:end) = 0;
    case 'random'
        if exist('random_seed', 'var')
            rng(random_seed)
        end
        incr = rand(Nenc-1,1)*180;
    otherwise
        error(['method_across_enc named ' upper(method_across_enc) ' is not implemented'])
end

angles = angles_enc1;
incr = cat(1,0,incr);
for enc = 1:Nenc
    angles = angles + incr(enc);
    points = reshape((2*pi)*linspace(-0.5,0.5,256)'.*(exp(+1j*pi.*(angles)./180)),[],1);
    k(:,1,enc,:) = cat(3, real(points),imag(points));

end



end

