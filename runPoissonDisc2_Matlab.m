
function runPoissonDisc2_Matlab
  close all;  clear;  rng(1);

  addpath( genpath( '.' ) );

  datacase = 1;
  pattern = 1;
  verbose = false;
  kDistThresh = 0.05;
  showScale = 1;
  nAvgs = 1;


  [r,min_r,max_r,img,lambda] = loadDatacase( datacase, pattern );


  avgTimes = zeros(nAvgs,1);
  avgTimesTulleken = zeros(nAvgs,1);
  for avgIndx = 1 : nAvgs
    tic
    if ispositive( r )
      if numel( r ) > 1
        bounds = [ [ -0.5, 0.5 ]; ...    % horizontal bounds
                   [ -0.5, 0.5 ]; ];     % vertical bounds

        ks = poissonDisc2( r, 'bounds', bounds );
      else
        ks = poissonDisc2( r );
      end
    else
      ks = poissonDisc2( r, 'min_r', min_r );
    end
    timeTaken = toc;
    avgTimes(avgIndx) = timeTaken;
    disp([ 'Time taken for sample generation: ', num2str(timeTaken), ' (s)' ]);

    if ~ispositive( r ) || numel( r ) > 0
      tic
      ksTulleken = poissonDisc2_tulleken( r, max_r );
      timeTakenTulleken = toc;
      avgTimesTulleken(avgIndx) = timeTakenTulleken;
    end
    disp([ 'Time taken for Tulleken sample generation: ', num2str(timeTakenTulleken), ' (s)' ]);

    ksAreas = voronoiAreas( ks );
    krs = sqrt( ks(:,1).^2 + ks(:,2).^2 );
    krs = krs( isfinite( ksAreas ) );
    ksAreas = ksAreas( isfinite( ksAreas ) );
    figure; s = scatternice( krs, ksAreas, 2, 'b' );
    s.MarkerFaceAlpha = 0.3;

    tAreas = voronoiAreas( ksTulleken );
    ksTrs = sqrt( ksTulleken(:,1).^2 + ksTulleken(:,2).^2 );
    ksTrs = ksTrs( isfinite( tAreas ) );
    tAreas = tAreas( isfinite( tAreas ) );
    hold all;  s = scatternice( ksTrs, tAreas, 2, 'r' );
    s.MarkerFaceAlpha = 0.3;
    ylim([ 0 0.00005 ]);
    legend( 'Fast Algorithm', 'Tulleken' );
    xlabel('Distance From Origin');  ylabel('Area of Voronoi Cell')
  end
  disp([ 'Avg time taken: ', num2str(mean(avgTimes)), ' (s)' ]);
  disp([ 'Tulleken Avg time taken: ', num2str(mean(avgTimesTulleken)), ' (s)' ]);

  kPts = ks';
  Ns = size( img );
  dks = 1 ./ (Ns-1);
  [~,samples] = movePointsToGrid( kPts, [-0.5, -0.5], 0.5-dks, Ns );
  sampleMask = samples > 0;

  if exist( 'ph', 'var' )
    figure; imshowscale( ph, showScale );
  end
  figure; imshowscale( sampleMask, showScale );  titlenice( 'Sample Mask' );
  disp([ 'Sample percentage: ', num2str( sum(sampleMask(:)) / numel( sampleMask) ) ]);

  kPtsTulleken = ksTulleken';
  [~,samplesTulleken] = movePointsToGrid( kPtsTulleken, [-0.5, -0.5], 0.5-dks, Ns );
  sampleMaskTulleken = samplesTulleken > 0;
  figure; imshowscale( sampleMaskTulleken, showScale );  titlenice( 'Tulleken Mask' );
  disp([ 'Tulleken sample percentage: ', num2str( sum(sampleMaskTulleken(:)) / ...
    numel( sampleMaskTulleken) ) ]);

  figure; showImageCube( cat( 3, sampleMask, sampleMaskTulleken ), ...
    'border', 5, 'borderValue', 'max' );

  coords = size2fftCoordinates( size(sampleMask) );
  [kxs,kys] = meshgrid( coords{1}, coords{2} );
  augMask = sampleMask;
  augMask( abs(kys) < kDistThresh & abs(kxs) < kDistThresh ) = 1;
  figure; imshowscale( augMask );  titlenice( 'Augmented Mask' );
  samples = augMask .* fftc( img );
  disp([ 'Aug Mask Sample Percentage: ', num2str( sum(augMask(:)) / numel( augMask) ) ]);

  augMaskTulleken = sampleMaskTulleken;
  augMaskTulleken( abs(kys) < kDistThresh & abs(kxs) < kDistThresh ) = 1;
  figure; imshowscale( augMaskTulleken );  titlenice( 'Augmented Mask Tulleken' );
  samplesTulleken = augMaskTulleken .* fftc( img );
  disp([ 'Aug Mask Tulleken Sample Percentage: ', ...
    num2str( sum(augMaskTulleken(:)) / numel( augMaskTulleken) ) ]);
  
  csRecon = csReconFISTA( samples, lambda, 'verbose', verbose );
  figure; imshowscale( abs(img), showScale );  titlenice( 'Original Image' );
  figure; imshowscale( abs(csRecon), showScale, 'range', abs(img) );  titlenice( 'CS Recon' );
  figure;
  showImageCube( cat( 3, img, csRecon ), showScale, 'border', 5, 'borderValue', 'max' );
  titlenice( 'New Sample scheme' );

  csReconTulleken = csReconFISTA( samplesTulleken, lambda, 'verbose', verbose );
  figure; imshowscale( abs(csReconTulleken), showScale, 'range', abs(img) );
  titlenice( 'CS Recon Tulleken' );
  figure;
  showImageCube( cat( 3, img, csReconTulleken ), showScale, 'border', 5, 'borderValue', 'max' );
  titlenice( 'Tulleken Sample scheme' );
end


function [r,min_r,max_r,img,lambda] = loadDatacase( datacase, pattern )

  if datacase == 1
    img = rgb2gray( single( imread( './bac.jpg' ) ) / 255 );
    img = imresize( img, [512 512], 'bilinear' );
    lambda = 1d-3;

  else
    load ../../data/brain.mat;
    img = im;  clear im;
    lambda = 1d-3;

  end

  Ns = size( img );

  switch pattern
    case 0
      r = 0.03;  min_r = r;  max_r = r;

    case 1
      Delta = 0.15;
      m = 150;
      %r = @(x) max( ( ( norms( x, 2, 2 ) + 0.2 ) / 100 ), min_r );
      r = @(x) ( norms( x, 2, 2 ) + Delta ) / m;
      min_r = Delta / m;
      max_r = r( [ 0.5 0.5 ] );

    case 2
      max_r = 0.1;
      min_r = 0.004;
      ph = imresize( phantom(), Ns, 'bilinear' );
      r = min_r + ph / max( ph(:) ) * ( max_r - min_r );

    case 3
      min_r = 0.004;
      r = @(x) abs( x(:,1) ) / 100 + min_r;
      max_r = r( [ 0.5 0.5 ] );
  end

end

