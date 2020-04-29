
function pts = poissonDisc2_tulleken( r, max_r, varargin )
  % Variable density poisson disc sampling
  %
  % Written according to "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson
  % and "Poisson Disk Sampling" written by Herman Tulleken located at
  % http://devmag.org.za/2009/05/03/poisson-disk-sampling/
  %
  % Written by Nicholas Dwork

  p = inputParser;
  p.addRequired( 'r' );
  p.addRequired( 'max_r', @ispositive );
  p.addParameter( 'bounds', [], @isnumeric );
  p.addParameter( 'bIncSize', 50, @ispositive );
  p.addParameter( 'nK', 10, @ispositive );
  p.addParameter( 'incSize', 5000, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( r, max_r, varargin{:} );
  r = p.Results.r;
  max_r = p.Results.max_r;
  bounds = p.Results.bounds;
  bIncSize = p.Results.bIncSize;
  nK = p.Results.nK;
  incSize = p.Results.incSize;
  verbose = p.Results.verbose;

  if numel( bounds ) == 0
    bounds = [ [ -0.5, 0.5 ];     % vertical bounds
               [ -0.5, 0.5 ]; ];  % horizontal bounds
  end


  nd = 2;  % Number of dimensions
  
  if ispositive( r )
    rImg = r;
    sImg = size( rImg );

    if numel( bounds ) > 0
      oldBounds = bounds;
      oldBDists = oldBounds(:,2) - oldBounds(:,1);
    end
    bounds = [ [ 1, sImg(2) ]; ...   % horizontal bounds
               [ 1, sImg(1) ]; ];    % vertical bounds

    rScaling = mean( sImg(:) ./ oldBDists(:) );
    rImg = rImg * rScaling;
    max_r = max( rImg(:) );

    r = @(pts) rImg( round(pts(:,1)) + round( pts(:,2) - 1 ) * sImg(1) );
  end

  
  bDists = bounds(:,2) - bounds(:,1);

  pts = zeros( incSize, nd );
  nPts = 1;
  pts(nPts,:) = rand(1,nd) .* bDists' + bounds(:,1)';

  pts_r = zeros( incSize, 1 );
  pts_r(nPts) = r( pts(1,:) );

  cellSize = max_r / sqrt(nd);
  sBGrid = ceil( bDists' ./ cellSize );
  bGrid = cell( sBGrid );  % background grid
  tmp = zeros( bIncSize, 1 );
  bGrid(:) = {tmp};
  nBEls = bIncSize * ones( sBGrid );  % number of elements in each background grid cell
  nBPts = zeros( sBGrid );  % number of points in each background grid cell

  gc = getGridCoordinate( pts(1,:), bounds, cellSize );
  bGrid{ gc(1), gc(2) }(1) = 1;
  nBPts( gc(1), gc(2) ) = 1;

  pList = zeros( incSize, 1 );  % processing list
  nProcPts = 1;
  pList( nProcPts ) = 1;

  while nProcPts > 0

    %-- Choose a point at random from the processing list
    pIndx = ceil( rand() * nProcPts );  % index into the processing list
    aIndx = pList( pIndx );  % active point index
    aPt = pts( aIndx, : );
    active_r = pts_r( aIndx );

    %-- Make k random points in the annulus between r and 2r from the active point
    validPointFound = false;
    kDists = rand(nK,1) * active_r + active_r;
    kAngles = rand(nK,1) * 2*pi;
    kPts = [ aPt(1) + kDists .* cos( kAngles ),  ...
             aPt(2) + kDists .* sin( kAngles ) ];
    kPts = kPts( kPts(:,1) > bounds(1,1) & kPts(:,1) < bounds(1,2) & ...
                 kPts(:,2) > bounds(2,1) & kPts(:,2) < bounds(2,2), : );

    if numel( kPts ) > 0
      kgcs = getGridCoordinate( kPts, bounds, cellSize );
      krs = r( kPts );

      %-- Check to see if each candidate point is valid
      for kIndx = 1 : size( kPts, 1 )
        kPt = kPts( kIndx, : );
        kgc = kgcs( kIndx, : );
        kr = krs( kIndx );
        distSqThresh = kr * kr;

        % find the minimum distance to other nearby points
        nNearCells = ceil( kr / cellSize );
        nearbyGrid = bGrid( ...
          max( kgc(1) - nNearCells, 1 ) : min( kgc(1) + nNearCells, size(bGrid,1) ), ...
          max( kgc(2) - nNearCells, 1 ) : min( kgc(2) + nNearCells, size(bGrid,2) )  ...
        );
        nearbyGrid = cat( 1, nearbyGrid{:} );
          % The above line is much faster than nearbyGrid = cell2mat( nearbyGrid(:) );
        nearbyIndxs = nearbyGrid( nearbyGrid > 0 );
        nNearbyPts = numel( nearbyIndxs );
        if nNearbyPts > 0
          nearbyPts = pts( nearbyIndxs, : );
          tmp = nearbyPts - repmat( kPt, [ nNearbyPts 1 ] );
          dists2NearbyPtsSq = sum( tmp .* tmp, 2 );
          if min( dists2NearbyPtsSq ) < distSqThresh, continue; end
        end

        validPointFound = true;

        % add this point to the list of points
        nPts = nPts + 1;
        if size(pts,1) < nPts
          if verbose ~= false
            disp([ 'poissonDisc2: made ', num2str(nPts-1), ' points.' ]);
          end
          pts = [ pts; zeros(incSize,nd) ];   %#ok<AGROW>
        end
        pts( nPts, : ) = kPt;
        pts_r( nPts ) = kr;

        % add this point to the processing list
        nProcPts = nProcPts + 1;
        if numel( pList ) < nProcPts
          pList = [ pList; zeros(incSize,1); ];   %#ok<AGROW>
        end
        pList( nProcPts ) = nPts;

        % add this point into the background grid
        thisNBPts = nBPts( kgc(1), kgc(2) ) + 1;
        nBPts( kgc(1), kgc(2) ) = thisNBPts;
        if thisNBPts >= nBEls( kgc(1), kgc(2) )
          tmp = [ bGrid{ kgc(1), kgc(2) }; zeros( bIncSize, 1 ); ];
          tmp( thisNBPts ) = nPts;
          bGrid{ kgc(1), kgc(2) } = tmp;
          nBEls( kgc(1), kgc(2) ) = nBEls( kgc(1), kgc(2) ) + bIncSize;
        else
          bGrid{ kgc(1), kgc(2) }(thisNBPts) = nPts;
        end
      end
    end

    if validPointFound == false
      pList( pIndx ) = pList( nProcPts );
      nProcPts = nProcPts - 1;
    end
  end

  if exist( 'oldBounds', 'var' )
    pts(:,1) = ( pts(:,1) - bounds(1,1) ) / rScaling + oldBounds(1,1);
    pts(:,2) = ( pts(:,2) - bounds(2,1) ) / rScaling + oldBounds(2,1);
  end

  pts = pts( 1 : nPts, : );
end


function out = ispositive( x )
  out = isnumeric(x) && min( x(:) > 0 );
end


function gc = getGridCoordinate( pts, bounds, cellSize )
  gc = ceil( [ ( pts(:,1) - bounds(1,1) ) , ...
               ( pts(:,2) - bounds(2,1) ) ] / cellSize );
end

