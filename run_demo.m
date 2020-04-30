
function run_demo
    gamma = 150;       % Increase this to increase the number of points
    accelerateX = 1;   % Increase this to reduce the points in the x direction
    accelerateY = 2;   % Increase this to reduce the points in the y direction

    addpath( 'dworkLib' );

    Ns = [ 512 512 ];
    Delta = 0.15;
    r = @(x) ( norms( x, 2, 2 ) + Delta ) / gamma;
    min_r = r( [0 0] );

    bounds = [ [ -0.5, 0.5 ] / accelerateX;     % vertical bounds
               [ -0.5, 0.5 ] / accelerateY; ];  % horizontal bounds
    ks = poissonDisc2( r, 'min_r', min_r, 'bounds', bounds );
    ks(:,1) = ks(:,1) * accelerateX;
    ks(:,2) = ks(:,2) * accelerateY;
    dks = 1 ./ (Ns-1);
    [~,samples] = movePointsToGrid( ks', [-0.5, -0.5], 0.5-dks, Ns );
    sampleMask = samples > 0;
    figure;  imshow( imresize( sampleMask, 5, 'nearest'), [] );
end



function out = norms( A, p, dim )
  % out = norms( A [, p, dim ] )
  %
  % Calculates the p norm of A along dimension dim
  %
  % Inputs:
  % A - a matrix of dimension <= dim
  %
  % Optional Inputs:
  % p - Compute the Lp norm (default is 2)
  % dim - (default is last dimension)
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 3, dim = ndims(A); end
  if nargin < 2, p = 2; end

  absA = abs( A );

  if mod( p, 1 ) == 0
    % p is an integer
    tmp = absA;
    for i=2:p
      tmp = tmp .* absA;
    end    
  else
    % p is not an integer
    tmp = absA.^p;
  end

  tmp = sum( tmp, dim );
  out = tmp .^ (1/p);
end
