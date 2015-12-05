#-------------------------------------------------------------------------------
# Name:        nirams_ii_main.py
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

def main():
    """ User-specified options and parameters for NIRAMS II. Edit the code 
        below to set up the model.
    """
    import time
    
    # #########################################################################
    # Get user input
    # HDF5 path containing inut data
    hdf5_path = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Input_Files\nirams_ii_input.h5'
    
    # Period of interest (years)
    st_yr, end_yr = 2001, 2001
    
    # Area of interest (co-ordinates should be integer multiples of 5000m)
    # National = 0, 485000, 520000, 1235000
    xmin, xmax, ymin, ymax = 0, 485000, 520000, 1235000
    
    # Choose the default PET to AET conversion factor grid
    pet_grid = 'lcms88_pet_fact'
    
    # Choose whether to use the iacs_pet_fact grids instead of the one above 
    # for the years from 2001 to 2010 inclusive
    use_iacs_pet_fact = True
    
    # Choose whether to write output to GeoTiffs as well as to HDF5 (slower, 
    # but GeoTiffs can be added directly to ArcGIS)
    write_gtif = True
    
    # Run ID. You must assign a unique integer to identify this model run
    run_id = 1
    
    # Output folders
    hdf5_fold = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Input_Files'
    gtiff_fold = r'C:\Users\user\Documents\JHI\NIRAMS_2015\NIRAMS_Input_Files\GeoTiffs'
    
    # Physical parameters
    t_snow = -2             # Threshold temp. for snowfall to begin in degrees
    t_melt = 1              # Threshold temp. for snow melt to begin in degrees
    ddf = 10                # Degree-day factor in mm/C/day
    org_n_fact = 0.2        # Fraction of organic N immediately available for
                            # leaching
    calib_min = 0.15        # Parameter controlling mineralistion rate
    calib_denit = 0.01      # Parameter controlling denitrification rate
    calib_n_leach = 1.0     # Parameter controlling N leaching rate
    # #########################################################################
    
    st_time = time.time()
    
    # Add model params to dict
    params_dict = {'Input HDF5 path'          : hdf5_path,
                   'Output HDF5 folder'       : hdf5_fold,
                   'Write GeoTiffs'           : write_gtif,
                   'Output GeoTiff folder'    : gtiff_fold,
                   'Run ID'                   : run_id,
                   'Start year'               : st_yr,
                   'End year'                 : end_yr,
                   'xmin'                     : xmin,
                   'xmax'                     : xmax,
                   'ymin'                     : ymin,
                   'ymax'                     : ymax,
                   'Default PET to AET grid'  : pet_grid,
                   'Use IACS'                 : use_iacs_pet_fact,
                   'T_snow'                   : t_snow,
                   'T_melt'                   : t_melt,
                   'Degree-day factor'        : ddf,
                   'Organic N factor'         : org_n_fact,
                   'Mineralisation parameter' : calib_min,
                   'Denitrification parameter': calib_denit,
                   'N leaching parameter'     : calib_n_leach}
    
    # Run model
    run_nirams_ii(params_dict)
    
    end_time = time.time()
    print 'Finished. Processing time: %.2f minutes.' % ((end_time-st_time)/60)

def run_nirams_ii(params_dict):
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

if __name__ == "__main__":
    main()