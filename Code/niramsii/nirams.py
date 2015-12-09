#-------------------------------------------------------------------------------
# Name:        nirams.py
# Purpose:     The main script for NIRAMS II.
#
# Author:      James Sample
#
# Created:     18/01/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" This is the main or parent script for NIRAMS II.
"""

def run_model(param_xls, run_type):
    """ Run the NIRAMS II model wtht he specified parameters and options.
    
    Args:
        param_xls: Path to complete Excel template for model setup.
        run_type:  The name of the worksheet to read data from in 'param_xls'.
                   Either 'single_run' for a single model run with the  
                   specified parameter set or 'param_combos' for multiple 
                   model runs.
    """
    import time, input_output as io
    
    st_time = time.time()

    # Check run_type is valid
    assert run_type in ('single_run', 'param_combos'), (
        "The run_type must be either 'single_run' or 'param_combos'.\n"
        "Please check and try again.")
             
    # Read user input from specified sheet
    user_dict = io.read_user_input(param_xls, run_type)
        
    # Determine which version of NIRAMS to run
    if run_type == 'single_run':
        print 'Running NIRAMS II for a single parameter set...'   

        # Run single model
        single_niramsii_run(user_dict)
    
    else:
        # run_type = 'param_combos'
        import multiprocessing as mp, pandas as pd, os, shutil
        import input_output as io

        print 'Running NIRAMS II for multiple parameter sets...'  

        # Create temp folder for intermediate output
        temp_path = os.path.split(user_dict['Output HDF5 path'])[0]
        temp_fold = os.path.join(temp_path, 'temp')
        os.makedirs(temp_fold)
        
        # Calculate parameter combinations based on user input
        param_combos = calculate_combinations(user_dict)
        
        # Write param combos to CSV
        df = pd.DataFrame(data=param_combos)
        df.index = df['Run ID']
        del df['Run ID']
        
        # Re-order columns
        col_order = ['Number of processors', 'Input HDF5 path', 'Parameter CSV',
                     'Output HDF5 path', 'Write GeoTiffs', 
                     'Output GeoTiff folder', 'Start year', 'End year',
                     'xmin', 'xmax', 'ymin', 'ymax', 'Default PET to AET grid',
                     'Use IACS', 'T_snow', 'T_melt', 'Degree-day factor',
                     'Organic N factor', 'Mineralisation parameter',
                     'Denitrification parameter', 'N leaching parameter']
        df[col_order].to_csv(user_dict['Parameter CSV'][0], 
                             index_label='Run ID')

        # Add the temp folder to each param_dict
        for idx, param_dict in enumerate(param_combos):
            param_dict['Output HDF5 folder'] = temp_fold
    
        # Setup multiprocessing pool
        pool = mp.Pool(user_dict['Number of processors'][0])
        
        # Distribute runs to processors
        pool.map(single_niramsii_run, param_combos)
        pool.close()
    
        # Merge outputs from each run to single HDF5
        io.merge_hdf5(temp_fold, user_dict['Output HDF5 path'][0])
        
        # Remove temp folder
        shutil.rmtree(temp_fold)
        
    end_time = time.time()
    print 'Finished. Processing time: %.2f minutes.' % ((end_time-st_time)/60)        
    
def single_niramsii_run(params_dict):
    """ Run the NIRAMS II model with the specified parameters.
    
    Args:
        params_dict: Dict containing all model parameters and user-specified
                     options.
    
    Returns:
        None. Model outptus are written to the specified HDF5 file (and 
        GeoTiffs if specified).
    """
    import input_output as io, snow as sn, drainage as dr
    import nitrate as ni, calendar, numpy.ma as ma, os
                   
    # Paths to static nodes in the input HDF5 file
    nodes_dict = {'land_props' : r'/one_km_grids/old_land_properties/',
                  'soil_props' : r'/one_km_grids/soil_properties/',
                  'met_data'   : r'/five_km_grids/meteorological_data/',
                  'iacs_pet'   : r'/one_km_grids/iacs_pet_facts/',
                  'or'         : r'/one_km_grids/organic_n/',
                  'in'         : r'/one_km_grids/inorganic_n/',
                  'up'         : r'/one_km_grids/n_uptake/',
                  'n_dep'      : r'/one_km_grids/n_deposition/',
                  'time_series': r'/time_series/'}
                      
    # Create output HDF5 file
    io.create_output_h5(params_dict)
    
    # Dicts storing number of days in each month (one for leap years; one for 
    # non-leap years)
    days_in_month_dict = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30,
                          10:31, 11:30, 12:31}
    days_in_month_lpyr_dict = {1:31, 2:29, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 
                               9:30, 10:31, 11:30, 12:31}
    
    # Extract the grid indices for the bounding box into a dict
    indices_dict = io.get_grid_indices(
                                      params_dict['xmin'], params_dict['xmax'],
                                      params_dict['ymin'], params_dict['ymax'])
    
    # Extract the static grids from the HDF5 file
    fc, sat, calibl, calibv = io.read_static_grids(
                                             params_dict['Input HDF5 path'], 
                                             nodes_dict['soil_props'], 
                                             ['fc', 'sat', 'calibl', 'calibv'],
                                             indices_dict)
    
    # Extract the PET to AET correction factor grid from the HDF5 file
    default_pet_fact = io.read_static_grids(
                                     params_dict['Input HDF5 path'], 
                                     nodes_dict['land_props'], 
                                     [params_dict['Default PET to AET grid'],],
                                     indices_dict)[0]
    
    # Set an initial water level halfway between field and saturation capacity
    wat_lev = (fc + sat)/2
    
    # Set an initial snow pack of zero
    rows = (params_dict['ymax']-params_dict['ymin'])/1000
    cols = (params_dict['xmax']-params_dict['xmin'])/1000
    snow_pk = ma.zeros((rows,cols))
    
    # Set the initial amount of available N using a simple annual balance for
    # 2001
    # Get the annual N grids for 2001 in a dict
    n_bud_dict = io.read_annual_n_grids(params_dict['Input HDF5 path'], 
                                        nodes_dict, 
                                        2001,
                                        indices_dict)
    avail_n = ni.initial_n_budget(n_bud_dict, params_dict['Organic N factor'])
    
    # Begin looping over time series data
    for year in range(params_dict['Start year'], params_dict['End year']+1):
        # Choose PET to AET conversion grids based on user input
        if (params_dict['Use IACS'] == True) and (year in range(2001, 2011)):
            # Get the iacs_pet_fact grid for this year
            pet_fact = io.read_static_grids(params_dict['Input HDF5 path'], 
                                            nodes_dict['iacs_pet'],
                                            ['pet_fact_%s' % year,], 
                                            indices_dict)[0]
        else:
            # Use the default pet_fact grid
            pet_fact = default_pet_fact
    
        # Read the annual N grids
        annual_n_dict = io.read_annual_n_grids(params_dict['Input HDF5 path'], 
                                               nodes_dict, 
                                               year,
                                               indices_dict)
    
        # Calculate daily n_dep rate for this year
        if calendar.isleap(year) == True:
            daily_n_dep = annual_n_dict['n_dep'] / 366.
        else:
            daily_n_dep = annual_n_dict['n_dep'] / 365.
    
        # Keep track of annual totals
        an_n_leach = ma.zeros((rows,cols))
        an_ssf = ma.zeros((rows,cols))
        an_gwf = ma.zeros((rows,cols))
        an_of = ma.zeros((rows,cols))
    
        # Loop over months
        for month in range(1,13):   
            # Allow for leap years
            if calendar.isleap(year) == True:
                days_in_month = days_in_month_lpyr_dict[month]
            else:
                days_in_month = days_in_month_dict[month]
    
            # Loop over days
            for day in range(1, days_in_month+1):
                # Get today's met data from the HDF5 file
                pptn, t_min, t_max, pet = io.read_met_data(
                                                params_dict['Input HDF5 path'],
                                                nodes_dict['met_data'],
                                                indices_dict,
                                                year,
                                                month,
                                                day,
                                                days_in_month)
    
                # Convert PET to AET using pet_fact
                aet = pet_fact*pet
    
                # Where the ground is already covered in snow, set AET to zero
                aet[snow_pk>0] = 0
    
                # Reduce the AET if the soil is dry i.e. if wat_lev < 0.7*fc
                aet = dr.reduce_aet_if_dry(aet, wat_lev, fc)
    
                # Split today's pptn into rain and snow components
                rain, snow = sn.estimate_snow_and_rain(pptn, t_min, t_max, 
                                                       params_dict['T_snow'])
    
                # Calculate today's snow melt
                melt = sn.estimate_snow_melt(snow_pk, t_min, t_max, 
                                             params_dict['T_melt'], 
                                             params_dict['Degree-day factor'])
    
                # Estimate temp and moisture factors
                t_fact = ni.est_temp_factor(t_min, t_max)
                moist_fact = ni.est_moisture_fact(wat_lev, fc)
    
                # Calculate today's mineralisation
                n_mineral = ni.est_mineralisation(
                                       params_dict['Mineralisation parameter'], 
                                       t_fact, 
                                       moist_fact)
    
                # Calculate today's denitrification
                n_denit = ni.est_denitrification(
                                       params_dict['Denitrification parameter'], 
                                       wat_lev, 
                                       t_fact, 
                                       moist_fact, 
                                       avail_n)
    
                # Estimate amount of N added today
                ts_row = io.read_ts_table(params_dict['Input HDF5 path'], 
                                          nodes_dict['time_series'],
                                          day, 
                                          month)
                                          
                n_added = ni.estimate_n_added(annual_n_dict, 
                                              daily_n_dep, 
                                              params_dict['Organic N factor'], 
                                              n_mineral, 
                                              n_denit, 
                                              ts_row)
    
                # Calculate today's drainage grids
                dr_list = dr.estimate_drainage(fc, sat, calibl, calibv, 
                                               wat_lev, snow_pk, rain, snow,
                                               melt, aet)
                                               
                snow_pk, wat_lev, surf_ro, lat_dr, vert_dr, tot_dr = dr_list
    
                # Calculate today's N leaching
                n_leach_list = ni.calculate_n_leaching(
                                           avail_n, 
                                           n_added, 
                                           dr_list, 
                                           fc, 
                                           params_dict['N leaching parameter'])
                                           
                leached_n, avail_n = n_leach_list
    
                # Increment annual totals
                an_n_leach += leached_n
                an_gwf += vert_dr
                an_ssf += lat_dr
                an_of += surf_ro
    
        # Calculate yearly drainage
        an_drain = an_ssf+an_gwf+an_of
        an_ss_drain = an_ssf+an_gwf
    
        # Get path to output HDF5
        hdf5_fold = params_dict['Output HDF5 folder']
        run_id = params_dict['Run ID']
        out_hdf5 = os.path.join(hdf5_fold, 'run_%03d.h5' % run_id)
    
        # Write to output file
        # Total drainage    
        io.write_array_to_h5(out_hdf5,
                             '/run_%03d' % run_id,
                             'total_drainage_%s' % year,
                             an_drain,
                             units='mm', 
                             xmin=params_dict['xmin'], 
                             xmax=params_dict['xmax'], 
                             ymin=params_dict['ymin'], 
                             ymax=params_dict['ymax'])
    
        # Sub-surface drainage
        io.write_array_to_h5(out_hdf5,
                             '/run_%03d' % run_id,
                             'sub-surface_drainage_%s' % year,
                             an_ss_drain,
                             units='mm', 
                             xmin=params_dict['xmin'], 
                             xmax=params_dict['xmax'], 
                             ymin=params_dict['ymin'], 
                             ymax=params_dict['ymax'])
    
        # N leached
        io.write_array_to_h5(out_hdf5,
                             '/run_%03d' % run_id,
                             'n_leached_%s' % year,
                             an_n_leach,
                             units='mm', 
                             xmin=params_dict['xmin'], 
                             xmax=params_dict['xmax'], 
                             ymin=params_dict['ymin'], 
                             ymax=params_dict['ymax'])
    
        # Write to GTiff
        if params_dict['Write GeoTiffs'] == True:
            # Total drainage
            tot_dr_path = os.path.join(params_dict['Output GeoTiff folder'], 
                                       'run_%03d_total_drainage_%s.tif' 
                                       % (run_id, year))
            io.ma_to_gtiff(params_dict['xmin'], params_dict['ymax'], 1000, 
                           tot_dr_path, an_drain)
    
            # Sub-surface drainage
            ss_dr_path = os.path.join(params_dict['Output GeoTiff folder'], 
                                      'run_%03d_sub-surface_drainage_%s.tif' 
                                      % (run_id, year))
            io.ma_to_gtiff(params_dict['xmin'], params_dict['ymax'], 1000, 
                           ss_dr_path, an_ss_drain)
            
            # N leached
            n_leach_path = os.path.join(params_dict['Output GeoTiff folder'], 
                                        'run_%03d_n_leached_%s.tif' 
                                        % (run_id, year))
            io.ma_to_gtiff(params_dict['xmin'], params_dict['ymax'], 1000, 
                           n_leach_path, an_n_leach)

def calculate_combinations(user_dict):
    """ Read user inut from 'param_combos' sheet in input Excel file and 
        calculate all combnations of the parameters entered.
    
    Args:
        user_dict: Dict of user-defined input from the 'param_combos' sheet in
                   the input Excel template.
                   
    Returns:
        List of dicts, where each dict is a valid input for 
        single_niramsii_run()
    """
    import itertools as it
    
    # List of params for which multiple parameter values can be entered
    phys_params = ['T_snow', 'T_melt', 'Degree-day factor', 
                   'Organic N factor', 'Mineralisation parameter', 
                   'Denitrification parameter', 'N leaching parameter']
    
    # To work with itertools, all dict elements must be in lists. Also need to
    # Parse any comma separated lists entered for physical params               
    for key in user_dict.keys():
        if key in phys_params:
            if isinstance(user_dict[key], (float, int)):
                # User has just entered a single value
                user_dict[key] = [user_dict[key],]
            else:
                # User has entered a comma-separated list
                user_dict[key] = [float(i) for i in 
                                  user_dict[key].split(',')]
        else:
            # Just add the param directly to a list
            user_dict[key] = [user_dict[key],]
       
    # Generate combinations. See
    # http://stackoverflow.com/questions/3873654/combinations-from-dictionary-with-list-values-using-python
    # for details
    param_dicts = sorted(user_dict)
    param_combos = [dict(zip(param_dicts, prod)) 
                    for prod in it.product(*(user_dict[param_dict] 
                    for param_dict in param_dicts))]
    
    # Check n_runs < 1000 (because my file naming is padded to 3 digits, so
    # more than 999 runs won't work. Could be easily extended if necessary)
    assert len(param_combos) < 1000, ('The maximum numbers of runs for this '
                                      'code is 999.')
    
    # Assign unique run IDs from 1 to n for n param combos
    for idx, param_dict in enumerate(param_combos):
        param_dict['Run ID'] = idx + 1
    
    return param_combos