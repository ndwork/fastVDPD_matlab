
function run_makePaperImgs
  close all; clear; rng(1); addpath( 'ESPIRiT' );

  datacases = 2:-1:1;

  for datacase = datacases

    sub_kData = loadDatacase( datacase );

    maskDir = [ './samplingMasks/size_', indx2str(size(sub_kData,1),999), '_', ...
      indx2str(size(sub_kData,1),999) ];
    reconDir = [ './recons/', 'datacase_', indx2str(datacase,10) ];
    if ~exist( reconDir, 'dir' ), mkdir( reconDir ); end

    fftRecons = mri_fftRecon( squeeze(sub_kData(:,:,1,:)), 'multiSlice', true );
    recons2Show = abs( fftRecons ) / max( abs( fftRecons(:) ) );
    figure; showImageCube( recons2Show, 'nImgsPerRow', 4, 'range', [0 0.6] );
    saveas( gcf, [ reconDir, filesep(), 'fftRecons.png' ] );
    close all;

    ks = size2fftCoordinates( size( fftRecons(:,:,1) ) );
    [kxs,kys] = meshgrid( ks{1}, ks{2} );

    ssqRecon = mri_ssqRecon( sub_kData, 'multiSlice', true );
    figure; imshowscale( ssqRecon );
    saveas( gcf, [ reconDir, filesep(), 'ssqRecon.png' ] );
    close all;

    nCoils = size( fftRecons, 3 );

    maskFiles = dir( maskDir );
    maskFiles = maskFiles( 3 : end );

    parfor fileIndx = 1 : numel( maskFiles )
      dworkKData = zeros( size(fftRecons,1), size(fftRecons,2), 1, nCoils );
      tullekenKData = zeros( size(fftRecons,1), size(fftRecons,2), 1, nCoils );
      filename = maskFiles( fileIndx ).name;
      if regexp( filename, '^\.' ), continue; end
      if regexp( filename, '^dworkMask' ), continue; end
      dworkFilename = filename;
      filenameParts = split( filename, '_' );
      filenameParts{1} = 'tullekenMask';
      tullekenFilename = join( filenameParts(:), '_' );
      tullekenFilename = tullekenFilename{1};
      dworkFile2save = [ reconDir, filesep(), dworkFilename ];
      tullekenFile2save = [ reconDir, filesep(), tullekenFilename ];
      %if exist( file2save, 'file' ), continue; end

      dworkMask = double( imread( [ maskDir, filesep(), dworkFilename ] ) ) / 255.;
      %centerRegion = zeroOuterRegion( ones( size( mask ) ), 32 );
      %mask = mask | centerRegion;
      tmp = bsxfun( @times, squeeze( sub_kData ), dworkMask );
      dworkKData(:,:,1,:) = tmp;
      dworkSakeRecon = sakeRecon2D( dworkKData, 'type', 'espiritL1' );
      recon2Save = abs( dworkSakeRecon ) / max( abs( dworkSakeRecon(:) ) );
      imwrite( recon2Save, dworkFile2save );

      tullekenMask = double( imread( [ maskDir, filesep(), tullekenFilename ] ) ) / 255.;
      tmp = bsxfun( @times, squeeze( sub_kData ), tullekenMask );
      tullekenKData(:,:,1,:) = tmp;
      tullekenSakeRecon = sakeRecon2D( tullekenKData, 'type', 'espiritL1' );
      recon2Save = abs( tullekenSakeRecon ) / max( abs( tullekenSakeRecon(:) ) );
      imwrite( recon2Save, tullekenFile2save );

      dworkKxs = kxs( dworkMask == 1 );  tullekenKxs = kxs( tullekenMask == 1 );
      dworkKys = kys( dworkMask == 1 );  tullekenKys = kys( tullekenMask == 1 );
      dAreas = voronoiAreas( [ dworkKxs dworkKys ] );
      dRs = sqrt( dworkKxs.^2 + dworkKys.^2 );
      dAreas = dAreas( isfinite( dAreas ) );
      dRs = dRs( isfinite( dAreas ) );
      areaFig = figure;
      s = scatternice( dRs, dAreas, 2, 'b' );  s.MarkerFaceAlpha = 0.5;

      tAreas = voronoiAreas( [ tullekenKxs tullekenKys ] );
      tRs = sqrt( tullekenKxs.^2 + tullekenKys.^2 );
      tAreas = tAreas( isfinite( tAreas ) );
      tRs = tRs( isfinite( tAreas ) );

      hold all;  s = scatternice( tRs, tAreas, 2, 'r' );  s.MarkerFaceAlpha = 0.5;
      ylim([ 0 0.0002 ]);  legend( 'Fast Algorithm', 'Tulleken' );
      filenameParts{1} = 'areas';
      areasFilename = join( filenameParts(:), '_' );
      areasFilename = areasFilename{1};
      areasFile2Save = [ reconDir, filesep(), areasFilename ];
      saveas( areaFig, areasFile2Save ); close( areaFig );

    end
  end
end
