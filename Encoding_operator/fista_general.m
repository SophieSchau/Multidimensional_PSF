function [recon, cost] = fista_general(kdata, E, S, lambda_s, niter, step, varargin)

%   Based on fista.m by Mark Chiew  
%   May 2017
%
%   Implementation of the FISTA iteration scheme for L1-Regularised problems
%   Adapted from Beck and Teboulle, SIAM J Im Sci 2009
%
%   This will solve a problem of the form:
%
%   min{x} |E*x - d|_2 + lambda_s*|S*u|_1 + other l2 terms
%
%   where x (recon) is the estimate,  S is an invertible  transform (often Wavelet), 
%   or identity transform
%   E is the measurement operator, lambda_s is the L1-weighting
%   and d is the measured raw k-space data (kdata)

%   To make things easier, we work with our estimate x in the main loop
%   in the sparse domain, so that our measurement operator is actually 
%   effectively E*S_inverse
%   Then we simply transform by S_inverse as a final step
%   
%   Another way of viewing this is that we perform a change of variables
%   such that x_new = S*x
%
%   Constant stepsize variant
%
% Edited by Sophie Schauman August 2019

%   Initialise
d2 = E'*kdata;

x   =   S*(d2);
v   =   x;
t   =   1;
x_p  =   x;

cost = zeros(niter,1,'single');


l2_constraint_gradients = zeros(size(d2));

Sinv = S';
clear('S')

T = calcToeplitzEmbedding(E);

z = mtimes_Toeplitz(E, T, Sinv*v);


%   Parse remaining inputs
p   =   inputParser;

%   Input options
addParameter(p, 'lambdas',       []);
addParameter(p, 'l2_operators', {});
addParameter(p, 'function_inputs', {});
addParameter(p,'showProgress',0);

parse(p, varargin{:});
p   =   p.Results;

lambdas = p.lambdas;
l2_operators = p.l2_operators;
function_inputs = p.function_inputs;
showProgress = p.showProgress;

if length(lambdas) ~= length(l2_operators)
error('need to the same number of lambdas as additional l2 terms')
end





fprintf(1, '%-5s %-16s %-16s %-16s %-16s\n', 'Iter','DataCon','L1','L2','Cost');

%   Main loop
for ii = 1:niter
    
    % additional l2 constraints
    l2_constraint_gradients = zeros(size(d2));
    for n = 1:length(l2_operators)  
        l2_constraint_gradients = l2_constraint_gradients+...
            lambdas(n) * l2_operators{n}(Sinv*v, function_inputs{n}{:}); 
    end
    
    %   Data consistency
    x   =   v + step * (Sinv' * ( d2       -   z    - l2_constraint_gradients));
   %x   =   v + step * (Sinv' * ( E'*kdata - E'*E*v - ...)
       
    %   Solve proximal sub-problem
    x   =   shrink(x, lambda_s*step);

    %   Compute momentum parameter
    t2  =   (1+sqrt(1+4*t^2))/2;

    %   Compute momentum update
    v   =   x + ((t-1)/t2)*(x-x_p);
%    v = u; %ista

    
    %   Update variables
    t   =   t2;
    x_p  =   x;
    z   =   mtimes_Toeplitz(E, T, Sinv*v);
   if showProgress 
      if mod(ii,10) == 0  

          im = Sinv*x;
        %   Error terms and cost function
        err1(ii)  =   reshape(im,[],1)'*(z(:)-2*d2(:))+kdata(:)'*kdata(:);
        err2(ii)  =   lambda_s*norm(x(~isnan(x)),1);
        err3(ii)  =   sum(abs(reshape(l2_constraint_gradients,[],1)).^2);
        cost(ii)  =   0.5*err1(ii) + err2(ii) + 0.5*err3(ii);

        %   Display iteration summary data
        fprintf(1, '%-5d %-16G %-16G %-16G %-16G\n', ii, err1(ii), err2(ii), err3(ii), cost(ii));
    %     imagesc(squeeze(max(abs(mean(u(:,:,:,:,1),4)),[],2)))
        figure(1)
        imagesc(max(squeeze(abs(im(:,:,:,1,1))),[],3))
        colorbar
        drawnow
%         figure(2)
%         imagesc(max(squeeze(abs(im(:,:,:,1,4))),[],3))
%         colorbar
%         drawnow
    %     
    %     imshow(abs(u(:,:,1,1,1)),[]);
    %     drawnow

      end
   else
      fprintf(['iteration: ' num2str(ii) '\n']);
   end
end

recon   =   Sinv*x;
end


function y = shrink(x, thresh)
    y = exp(1j*angle(x)).*max(abs(x)-thresh,0);
end
