#-------------------------------------------------------------------------------
# Name:        input_output.py
# Purpose:     Functions for controlling NIRAMS II input/output.
#
# Author:      James Sample
#
# Created:     17/01/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" Functions for reading arrays from HDF5 files and writing arrays to HDF5 and
    GeoTiff formats.
"""

def read_user_input(in_xls, sheet_name):
    """ Read the NIRAMS II input file.
    
    Args:
        in_xls:     Path to Excel template for NIRAMS II input.
        sheet_name: Sheet to read in Excel file. Either 'single_run' or
                    'param_combos'.
    
    Returns:
        Dict of user-specified model settings.
    """
    import pandas as pd
    
    # Read template
    df = pd.read_excel(in_xls, sheetname=sheet_name, skiprows=[0,1,2],
                       index_col=0)
    
    # Extract value of interest and convert to dict
    df = df['Value']
    param_dict = df.to_dict()    
    
    return param_dict

def merge_hdf5(hdf5_fold, out_path):
    """ Merge output HDF5 files from each run into a single file. All 
        attributes are also copied across.
        
    Args:
        hdf5_fold: Folder containing individual HDF5 files for each run. This
                   function will attempt to merge ALL HDF5 files in this 
                   folder, so make sure it only contains run output.
        out_path:  Path for merged HDF5 file.
    
    Returns:
        None. The datasets are merged and written to a single file.
    """
    import h5py, glob, os
       
    # Get a list of files to process
    search_path = os.path.join(hdf5_fold, '*.h5')
    file_list = glob.glob(search_path)

    # Open output HDF5
    out_h5 = h5py.File(out_path, 'w')
    
    # Loop over files
    for h5_path in file_list:
        # Get the run ID
        run_id = int(os.path.split(h5_path)[1].split('_')[1][:3])
        
        # Open file
        h5_file = h5py.File(h5_path, 'r')
        
        # Copy data
        h5_file.copy(r'/run_%03d' % run_id, out_h5)
        
        # Copy attributes from file to group
        grp = out_h5[r'/run_%03d' % run_id]
        for key, val in h5_file.attrs.iteritems():
            grp.attrs[key] = val 
    
        # Tidy up
        h5_file.close()
    out_h5.close()
    
def create_output_h5(params_dict):
    """ Create an output HDF5 file for a single NIRAMS II run. Store model
        parameter values as file attributes.
    
    Args:
        params_dict: Dict of user-specified params passed from main script.
    
    Returns:
        None. The specified file is created ready for adding data.
    """
    import h5py, os, datetime as dt
    
     # Create output H5
    hdf5_fold = params_dict['Output HDF5 folder']
    run_id = params_dict['Run ID']
    out_hdf5 = os.path.join(hdf5_fold, 'run_%03d.h5' % run_id)
    
    # Get current date and time
    now = dt.datetime.now()
    
    # Create the HDF5 file
    h5file = h5py.File(out_hdf5, mode='w',
                       title='NIRAMS II output. Created %s.'
                             % (now.strftime('%d/%m/%Y %H:%M')))
    
    # Add user attributes
    for key in params_dict.keys():
        h5file.attrs[key] = str(params_dict[key])
    
    h5file.close()                       
                             
def rebin(a, scale_factor):
    """ Resample an array to a new shape using the nearest neighbour technique. 
    
    Args:
        a:            Array to resample.
        scale_factor: Rescaling factor calculated as
                      orig_side_len / desired_side_len
    
    Returns:
        Resample array.
    """
    import numpy as np

    # Calculate new shape for array
    newshape = tuple([i*scale_factor for i in a.shape])

    # Perform resampling
    slices = [slice(0, old, float(old)/new) 
              for old, new in zip(a.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')

    return a[tuple(indices)]

def get_grid_indices(xmin, xmax, ymin, ymax):
    """ Take co-ordinates for bounding rectangle of interest, check they are 
        valid calculate array indices for desired area (for both a 1km and a 
        5km grid size).
    
    Args:
        xmin: Minimum Easting (in OSGB 1936 metres). Min 0; max 485000.
        xmax: Maximum Easting (in OSGB 1936 metres). Min 0; max 485000.
        ymin: Minimum Northing (in OSGB 1936 metres). Min 520000; max 1235000.
        ymax: Maximum Northing (in OSGB 1936 metres). Min 520000; max 1235000.
    
    Returns:
        Dict of grid indices.
    """
    # Check that co-ords lie in Scotland
    if (xmin < 0) or (xmax > 485000) or (ymin < 520000) or (ymax > 1235000):
        raise ValueError('The bounding box supplied does not lie within '
                         'Scotland.')

    # Check that the co-ordinates are an even multiple of 5000
    if ((xmin % 5000 != 0) or (xmax % 5000 != 0) or
        (ymin % 5000 != 0) or (ymax % 5000 != 0)):
        raise ValueError('Please enter co-ordinates that are an even multiple '
                         'of 5000m')

    # Calculate indices and add to dict
    indices_dict = {}
    indices_dict['xstart_1km'] = xmin / 1000
    indices_dict['xend_1km'] = xmax / 1000
    indices_dict['ystart_1km'] = (1235000-ymax) / 1000
    indices_dict['yend_1km'] = (1235000-ymin) / 1000
    indices_dict['xstart_5km'] = xmin / 5000
    indices_dict['xend_5km'] = xmax / 5000
    indices_dict['ystart_5km'] = (1235000-ymax) / 5000
    indices_dict['yend_5km'] = (1235000-ymin) / 5000

    return indices_dict
  
def read_static_grids(hdf5_file, node, wanted_list, indices_dict):
    """ Read the static, 1km grids from the HDF5 file.
    
    Args:
        hdf5_file:     Path to HDF5 file.
        node:          HDF5 node path.
        wanted_list:   List of grid names to read e.g. ['fc', 'sat'].
        indicies_dict: Dict from get_grid_indices().
        
    Returns:
        List of grids in same order specified in wanted_list.
    """
    import numpy.ma as ma, h5py

    # Open the HDF5 file
    h5file = h5py.File(hdf5_file, mode = "r")

    # Extract 1km indices from the indices dict
    xmin, xmax = indices_dict['xstart_1km'], indices_dict['xend_1km']
    ymin, ymax = indices_dict['ystart_1km'], indices_dict['yend_1km']

    # Empty list to store arrays
    grid_list = []

    for grid in wanted_list:
        # Get the array from the HDF5 file
        dset_path = node + grid
        dset = h5file[dset_path]
        array = dset[:].astype(float)
        
        # Read the array, mask no data
        array = ma.masked_equal(array, float(dset.attrs['No_Data_Value']))
        
        # Clip to BB
        array = 1.0*array[ymin:ymax, xmin:xmax]
        
        grid_list.append(array)

    h5file.close()

    return grid_list

def read_ts_table(hdf5_file, ts_node, day, month):
    """ Read a row from the time series table corresponding to a particular day
        and month.
    
    Args:
        hdf5_file: Path to HDF5 file.
        ts_node:   Node containing time series table.
        day:       Day of month (int)
        month:     Month (int)
        
    Returns:
        Structured array. Elements can be accessed as: str_arr['or_sp'] etc.
    """
    import h5py

    # Open the HDF5 file
    h5file = h5py.File(hdf5_file, mode = "r")

    # Get the table object
    dset_path = ts_node + 'time_series_table'
    table = h5file[dset_path]
    
    # Get the desired row
    ts_row = table[(table['day']==day) & (table['month']==month)]

    h5file.close()

    return ts_row

def read_annual_n_grids(hdf5_file, node_dict, year, indices_dict):
    """ Read the 1km organic, inorganic, crop uptake and N deposition grids 
        from the HDF5 file. 
    
    Args:
        hdf5_file:     Path to HDF5 file.
        node_dict:     Dict of node paths with keys ['in', 'or', 'up', 'n_dep']
        year:          Year (int)
        indicies_dict: Dict from get_grid_indices().        
    
    Returns:
        Dict of 13 grids: 4 each for organic, inorganic and crop uptake, and 
        one for deposition.
            dict[(N type, Land class)] = grid for org, inorg and uptake and
            dict['n_dep'] = grid for N deposition
    """
    import numpy.ma as ma, h5py

    # Open the HDF5 file
    h5file = h5py.File(hdf5_file, mode = "r")

    # Extract 1km indices from the indices dict
    xmin, xmax = indices_dict['xstart_1km'], indices_dict['xend_1km']
    ymin, ymax = indices_dict['ystart_1km'], indices_dict['yend_1km']

    # Empty dict to hold arrays
    annual_n_dict = {}

    # List of letter codes for broad land classes
    lnd_cls_list = ['gr', 'ot', 'sp', 'wi']

    # Process nodes
    for node in ['or', 'in', 'up', 'n_dep']:	
        for lnd_cls in lnd_cls_list:
            # n_dep requires slightly different processing
            if node != 'n_dep':
                # Get the array from the HDF5 file
                dset_path = node_dict[node] + '%s_%s_%s' % (node, lnd_cls, year)
                dset = h5file[dset_path]
                array = dset[:].astype(float)
				
                # Read the array, mask no data
                array = ma.masked_equal(array, 
                                        float(dset.attrs['No_Data_Value']))
				
                # Clip to BB
                array = array[ymin:ymax, xmin:xmax]

                # Add to dict
                annual_n_dict[(node, lnd_cls)] = array
            else:
                # For n_dep, the grid for 2010 isn't yet available, so re-use 
                # 2009 grid
                if year == 2010:
                    dset_path = node_dict[node] + 'n_dep_2009'
                else:
                    dset_path = node_dict[node] + 'n_dep_%s' % year
				
                # Get array
                dset = h5file[dset_path]
                array = dset[:].astype(float)
				
                # Read the array, mask no data
                array = ma.masked_equal(array, 
                                        float(dset.attrs['No_Data_Value']))
				
                # Clip to BB
                array = array[ymin:ymax, xmin:xmax]

                # Set masked values to default of 10
                array[array.mask] = 10

                # Add to dict
                annual_n_dict['n_dep'] = array

    h5file.close()

    return annual_n_dict

def read_met_data(hdf5_file, met_data_node, indices_dict, year, month, day,
                  days_in_month):
    """ Read the 5km meteorological data for the current time step. PET grids 
        are converted from monthly to daily values to match other variables.
        All grids are resampled to 1km resolution using the nearest neighbour
        technique.
    
    Args:
        hdf5_file:     Path to HDF5 file.
        met_data_node: Node in HDF5 file.
        indicies_dict: Dict from get_grid_indices().   
        year:          Year (int).
        month:         Month (int).        
        day:           Day of month (int).
        days_in_month: Number of days in month (int).
        
    Returns:
        List of grids: [rainfall, min_temp, max_temp, pet].
    """
    import numpy.ma as ma, h5py

    # Open the HDF5 file
    h5file = h5py.File(hdf5_file, mode = "r")

    # Extract 1km indices from the indices dict
    xmin, xmax = indices_dict['xstart_5km'], indices_dict['xend_5km']
    ymin, ymax = indices_dict['ystart_5km'], indices_dict['yend_5km']

    # List to store grids
    grid_list = []
    
    # List of variables to process
    var_list = ['rainfall', 'min_temp', 'max_temp', 'et']
    
    # Loop over variables
    for variable in var_list:
        # Get path to dataset
        if variable == 'rainfall':
            dset_path = (r'%sdaily/rainfall/rainfall_%s/rainfall_%s_%02d_%02d'
                         % (met_data_node, year, year, month, day))
        elif variable == 'min_temp':
            dset_path = (r'%sdaily/min_temp/min_temp_%s/min_temp_%s_%02d_%02d'
                         % (met_data_node, year, year, month, day))   
        elif variable == 'max_temp':
            dset_path = (r'%sdaily/max_temp/max_temp_%s/max_temp_%s_%02d_%02d'
                         % (met_data_node, year, year, month, day))
        elif (variable=='et') and ((year<1969) or (year>2005)):
            # Have to switch to Thornthwaite PET
            dset_path = (r'%smonthly/est_pm_pet/est_pm_pet_%s_%02d'
                         % (met_data_node, year, month))
        else:
            # P-M PET
            dset_path = (r'%smonthly/pm_pet/pm_pet_%s_%02d'
                         % (met_data_node, year, month))                         
                         
        # Get the array
        dset = h5file[dset_path]
        array = dset[:].astype(float)
    
        # Read the array and mask no data
        array = ma.masked_equal(array, float(dset.attrs['No_Data_Value']))                         

        # Clip to BB
        array = array[ymin:ymax, xmin:xmax]

        # Convert monthly to daily for PET
        if variable == 'et':
            array = array / days_in_month
            
        # Resample to 1km
        array = 1.0*rebin(array, 5)
               
        # Add to output
        grid_list.append(array)
        
    # Close the HDF5 file
    h5file.close()

    return grid_list

def write_array_to_h5(hdf5_file,
                      node,
                      array_name,
                      array,
                      units='Units',
                      srs_name='EPSG: 27700; OSGB 1936',
                      xmin=0,
                      xmax=485000,
                      ymin=520000,
                      ymax=1235000,
                      cell_size=1000,
                      no_data_value=-9999,
                      compress='gzip'):
    """ Write an array to an HDF5 file.
    
    Args:
        hdf5_file:     Path to HDF5 file.
        node:          HDF5 node to save to.
        array_name:    Name to use for svaed array.
        units:         String giving units for array.
        srs_name:      String describing co-ordinate system.
        xmin:          Minimum Easting (in OSGB 1936 metres). Min 0; max 485000.
        xmax:          Maximum Easting (in OSGB 1936 metres). Min 0; max 485000.
        ymin:          Minimum Northing (in OSGB 1936 metres). Min 520000; 
                       max 1235000.
        ymax:          Maximum Northing (in OSGB 1936 metres). Min 520000; 
                       max 1235000.
        cell_size:     Grid cell size in metres.
        no_data_value: Value to use to represent no data.
        compress:      Compression options. Either 'n' for no compression, or 
                       'lzf' or 'gzip'.
    """
    import h5py
	
    # Open HDF5 file
    h5 = h5py.File(hdf5_file, mode = "a")

	# Explicitly set the no data value to -9999
    array[array.mask] = no_data_value
	
	# Get dataset path
    dataset_path = r'%s/%s' % (node, array_name)
    
    if compress == 'n':    
        dset = h5.create_dataset(dataset_path, data=array)
        dset.attrs['Units'] = units
        dset.attrs['SRS_Name'] = srs_name
        dset.attrs['Bounding_Box'] = ("xmin, xmax, ymin, ymax = %s, %s, %s, %s"
                                      % (xmin,xmax,ymin,ymax))
        dset.attrs['Cell_Size'] = '%s m' % cell_size
        dset.attrs['No_Data_Value'] = str(no_data_value) 
        
    elif compress == 'lzf':           
        dset = h5.create_dataset(dataset_path, data=array,
                                 compression='lzf', fletcher32=True, 
                                 shuffle=True)
        dset.attrs['Units'] = units
        dset.attrs['SRS_Name'] = srs_name
        dset.attrs['Bounding_Box'] = ("xmin, xmax, ymin, ymax = %s, %s, %s, %s"
                                      % (xmin,xmax,ymin,ymax))
        dset.attrs['Cell_Size'] = '%s m' % cell_size
        dset.attrs['No_Data_Value'] = str(no_data_value)

    elif compress == 'gzip':           
        dset = h5.create_dataset(dataset_path, data=array, 
                                 compression='gzip', compression_opts=9,
                                 fletcher32=True, shuffle=True)
        dset.attrs['Units'] = units
        dset.attrs['SRS_Name'] = srs_name
        dset.attrs['Bounding_Box'] = ("xmin, xmax, ymin, ymax = %s, %s, %s, %s"
                                      % (xmin,xmax,ymin,ymax))
        dset.attrs['Cell_Size'] = '%s m' % cell_size
        dset.attrs['No_Data_Value'] = str(no_data_value)
        
    else:
        raise ValueError('Invalid compression option specified.')

    # Close
    h5.close()

def ma_to_gtiff(xmin, ymax, cell_size, out_path, data_array,
                no_data_value=-9999, srs_name='EPSG:27700'):
    """ Save numpy array as GeoTiff.
    
    Args:
        xmin:          Minimum Easting (in OSGB 1936 metres). Min 0; max 485000.
        ymax:          Maximum Northing (in OSGB 1936 metres). Min 520000; 
                       max 1235000.  
        cell_size:     Grid cell size in metres.
        out_path:      Path to GeoTiff.
        data:          Array to save.  
        no_data_value: Value to use to represent no data. 
        srs_name:      Valid GDAL co-ordinate system identifier.
        
    Returns:
        None. Array is saved to specified path.
    """
    # Import modules
    import gdal, gdalconst, osr

    # Explicitly set NDV
    data_array[data_array.mask] = no_data_value

    # Get array shape
    cols = data_array.shape[1]
    rows = data_array.shape[0]

    # Get driver
    driver = gdal.GetDriverByName('GTiff')

    # Create a new raster data source
    out_ds = driver.Create(out_path, cols, rows, 1, gdal.GDT_Float32)

    # Get spatial reference
    sr = osr.SpatialReference()
    sr.SetProjection (srs_name)
    sr_wkt = sr.ExportToWkt()

    # Write metadata
    # (xmin, cellsize, 0, ymax, 0, -cellsize)
    out_ds.SetGeoTransform((xmin, cell_size, 0.0, ymax, 0.0, -1*cell_size)) 
    out_ds.SetProjection(sr_wkt)
    out_band = out_ds.GetRasterBand(1)
    out_band.SetNoDataValue(-9999)
    out_band.WriteArray(data_array)

    # Tidy up
    del out_ds, out_band