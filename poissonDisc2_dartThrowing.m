
function pts = poissonDisc2_dartThrowing( r, N, varargin )
  % Written by Nicholas Dwork, Copyright 2019
  %
  % Inputs:
  % r - a scalar or a function
  % N - the number of points to find
  %
  % Optional Inputs:
  % bounds - the bounds on the region to find points
  % maxIter - the maximum number of iterations to conduct
  %   ( default is 100 * N )
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addRequired( 'r' );
  p.addRequired( 'N', @ispositive );
  p.addParameter( 'bounds', [], @isnumeric );
  p.addParameter( 'maxIter', [], @ispositive );
  p.parse( r, N, varargin{:} );
  r = p.Results.r;
  N = p.Results.N;
  bounds = p.Results.bounds;
  maxIter = p.Results.maxIter;

  nd = 2;

  if numel( bounds ) == 0
    bounds = [ [ -0.5, 0.5 ];     % vertical bounds
               [ -0.5, 0.5 ]; ];  % horizontal bounds
  end

  if numel( maxIter ) == 0, maxIter = 100 * N; end

  rIsFunction = true;
  if ispositive( r ) && numel( r ) == 1, rIsFunction = 0; end
  if rIsFunction == false, rPt = r; end

  figure; axis( [ bounds(1,1) bounds(1,2) bounds(2,1) bounds(2,2) ] );

  pts = zeros( N, nd );
  pts_r = zeros( N, 1 );
  nPts = 0;
  tic
  for i = 1 : maxIter
    if mod( i, 100 ) == 0
      timeTaken = toc;
      disp([ 'Working on iteration ', num2str(i), ' of ', num2str(maxIter), '.  ', ...
        'Num points: ', num2str(nPts), ' of ', num2str(N), '  time taken: ', ...
        num2str(timeTaken/60), ' minutes' ]);
    end

    ptUnscaled = rand(nd,1);
    pt = ptUnscaled .* ( bounds(:,2) - bounds(:,1) ) + bounds(:,1);
    pt = pt';

    if rIsFunction, rPt = r( pt ); end

    ptIsGood = true;
    if nPts > 0
      % Compute the distance between the current point and all other points
      dists = norms( pts(1:nPts,:) - repmat( pt, [nPts 1] ), 2, 2 );

      if min( dists ) < rPt || max( dists < pts_r(1:nPts) ) == 1
        ptIsGood = false;
      end
    end

    if ptIsGood == true
      nPts = nPts + 1;
      pts(nPts,:) = pt';
      pts_r(nPts) = rPt;

      hold on;  scatter( pts(nPts,1), pts(nPts,2), '.k' ); drawnow;
      if nPts >= N, break; end
    end
  end

  if i >= maxIter, disp('Max number of iterations reached'); end
  pts = pts(1:nPts,:);
end

