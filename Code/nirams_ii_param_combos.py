#-------------------------------------------------------------------------------
# Name:        nirams_ii_param_combos.py
# Purpose:     Run NIRAMS II for multiple parameters sets.
#
# Author:      James Sample
#
# Created:     18/01/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" This script takes lists of parameter values to consider and calculates
    parameter sets for all possible parameter combination in the specified 
    lists. These parameter sets are then distributed to multiple cores using 
    the multiprocessing module.
"""

def main():
    """ User-specified options and parameters for NIRAMS II. Edit the code 
        below to set up the model. 
        
        Any of the "physical parameters" can be entered as lists. The script
        will calculate parameter combinations from these lists and run the 
        model once for each parameter set. All the other options are kept the same
        for each run.
        
        Runs are automatically assigned a unique ID.
    """
    import itertools as it, nirams_ii_main as nirams, multiprocessing as mp
    import time, pandas as pd, os, shutil
    
    # #########################################################################
    # Get user input
    # Number of processors available for multiprocessing
    n_proc = 4
    
    # Output path for CSV of parameter combinations
    param_csv = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Output\param_combos.csv'
    
    # Output path for HDF5 file
    out_path = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Output\nirams_output.h5'
    
    # HDF5 path containing input data
    hdf5_path = r'C:\Users\user\Documents\JHI\NIRAMS_2015\nirams_ii_input.h5'
    
    # Period of interest (years)
    st_yr, end_yr = 2001, 2001
    
    # Area of interest (co-ordinates should be integer multiples of 5000m)
    # National = 0, 485000, 520000, 1235000
    xmin, xmax, ymin, ymax = 0, 485000, 520000, 1235000
    
    # Choose the default PET to AET conversion factor grid
    pet_grid = 'lcms88_pet_fact'
    
    # Choose whether to use the iacs_pet_fact grids instead of the one  
    # above for the years from 2001 to 2010 inclusive
    use_iacs_pet_fact = True
    
    # Choose whether to write output to GeoTiffs as well as to HDF5 
    # (slower, but GeoTiffs can be added directly to ArcGIS)
    write_gtif = False
       
    # Output folder for GeoTiffs
    gtiff_fold = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Output\GeoTiffs'
    
    # Physical parameters
    # NOTE: must all be lists e.g. [1,] for a single value
    t_snow = [-2,]                  # Threshold temp. for snowfall in degrees
    t_melt = [1,]                   # Threshold temp. for snow melt in degrees
    ddf = [10,]                     # Degree-day factor in mm/C/day
    org_n_fact = [0.2,]             # Fraction of organic N immediately 
                                    # available forleaching
    calib_min = [0.15,]             # Parameter controlling mineralistion rate
    calib_denit = [0.01, 0.02]      # Parameter controlling denitrification rate
    calib_n_leach = [1.0, 1.1, 1.2] # Parameter controlling N leaching rate
    # #########################################################################

    st_time = time.time()

    # Create temp folder for intermediate output
    temp_path = os.path.split(out_path)[0]
    temp_fold = os.path.join(temp_path, 'temp')
    os.makedirs(temp_fold)
    
    # Add model params to dict
    input_dicts = {'Input HDF5 path'          : [hdf5_path,],
                   'Output HDF5 folder'       : [temp_fold,],
                   'Write GeoTiffs'           : [write_gtif,],
                   'Output GeoTiff folder'    : [gtiff_fold,],
                   'Start year'               : [st_yr,],
                   'End year'                 : [end_yr,],
                   'xmin'                     : [xmin,],
                   'xmax'                     : [xmax,],
                   'ymin'                     : [ymin,],
                   'ymax'                     : [ymax,],
                   'Default PET to AET grid'  : [pet_grid,],
                   'Use IACS'                 : [use_iacs_pet_fact,],
                   'T_snow'                   : t_snow,
                   'T_melt'                   : t_melt,
                   'Degree-day factor'        : ddf,
                   'Organic N factor'         : org_n_fact,
                   'Mineralisation parameter' : calib_min,
                   'Denitrification parameter': calib_denit,
                   'N leaching parameter'     : calib_n_leach}
    
    # Generate combinations. See
    # http://stackoverflow.com/questions/3873654/combinations-from-dictionary-with-list-values-using-python
    # for details
    param_dicts = sorted(input_dicts)
    param_combos = [dict(zip(param_dicts, prod)) 
                    for prod in it.product(*(input_dicts[param_dict] 
                    for param_dict in param_dicts))]
    
    # Check n_runs < 1000 (because my file naming is padded to 3 digits, so
    # more than 999 runs won't work. Could be easily extended if necessary)
    assert len(param_combos) < 1000, ('The maximum numbers of runs for this '
                                      'code is 999.')
    
    # Assign unique run IDs from 1 to n for n param combos
    for idx, param_dict in enumerate(param_combos):
        param_dict['Run ID'] = idx + 1
    
    # Write param combos to CSV
    df = pd.DataFrame(data=param_combos)
    df.index = df['Run ID']
    del df['Run ID']
    df.to_csv(param_csv, index_label='Run ID')
    
    # Setup multiprocessing pool
    pool = mp.Pool(n_proc)
    
    # Ditribute runs to processors
    pool.map(nirams.run_nirams_ii, param_combos)
    pool.close()

    # Merge outputs from each run to single HDF5
    merge_hdf5(temp_fold, out_path)
    
    # Remove temp folder
    shutil.rmtree(temp_fold)
    
    end_time = time.time()
    print 'Finished. Processing time: %.2f minutes.' % ((end_time-st_time)/60)

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
            
if __name__ == "__main__":
    main()