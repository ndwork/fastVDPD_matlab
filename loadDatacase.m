
function sub_kData = loadDatacase( datacase )

  dataDir = '/Volumes/NDWORK128GB/';
  if ~exist( dataDir, 'dir' )
    dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/';
    if ~exist( dataDir, 'dir' )
      dataDir = '/Users/ndwork/Documents/Data/';
    end
  end


  if datacase == 1
    dataMatFile = 'knee_sub_kData.mat';
  elseif datacase == 2
    dataMatFile = 'ankle_sub_kData.mat';
  end
  
  if exist( dataMatFile, 'file' )
    load( dataMatFile, 'sub_kData' );
    return
  end

  if datacase == 1
    data = readOldMriDataOrgData( [ dataDir, '/mriData.Org/P14/kspace' ] );
    data = ifftc( data, [], 3 );
    sub_kData = data( :, :, 5:10:end, : );
    sub_kData = sub_kData( :, :, 4, : );

  elseif datacase == 2
    [data,header] = read_MR_rawdata( [ dataDir, '/cartesian3d/body3d/P27648.7' ] );   %#ok<ASGLU>
    data = squeeze( data );
    data = ifftc( data, [], 3 );
    sub_kData = data(:,:,50,:);
    sub_kData = rot90( sub_kData, -1 );
    
  end

  save( dataMatFile, 'sub_kData' );
end
