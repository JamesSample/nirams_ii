#-------------------------------------------------------------------------------
# Name:        nitrate.py
# Purpose:     Calculate nitrate leaching for NIRAMS II.
#
# Author:      James Sample
#
# Created:     06/08/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" Module containing functions to determine how much N is added at each time
    step and how much is leached away by the draining water. Many of these 
    equations are based on Sarah Dunn's STREAM-N model.
"""

def initial_n_budget(n_bud_dict, org_n_fact):
    """ Performs an rough annual N budget used to initialise the model.
    
    Args:
        n_bud_dict: Dict of N budgets for chosen initialisation year from
                    input_output.read_annual_n_grids()
        org_n_fact: User specified Organic N factor.
        
    Returns:
        Array of available N.
    """
    # avail_n = org_factor*org_n + inorg_n + n_dep + n_min - n_uptake - n_denit
    # but org, inorg and uptake all consist of 4 different grids (1 for each broad
    # land class). The pieces of this equation are built up below
    org_n = org_n_fact*(n_bud_dict[('or', 'gr')] + n_bud_dict[('or', 'sp')] +
                        n_bud_dict[('or', 'wi')] + n_bud_dict[('or', 'ot')])
                        
    inorg_n = (n_bud_dict[('in', 'gr')] + n_bud_dict[('in', 'sp')] +
               n_bud_dict[('in', 'wi')] + n_bud_dict[('in', 'ot')])
               
    n_uptake = (n_bud_dict[('up', 'gr')] + n_bud_dict[('up', 'sp')] +
                n_bud_dict[('up', 'wi')] + n_bud_dict[('up', 'ot')])
                
    n_dep = n_bud_dict['n_dep']
    
    n_min = 30.   # Rough starting values for mineralisation and 
    n_denit = 20. # denitrification based on Sarah's STREAM-N code. 
    
    # Perform the annual balance and tidy up
    avail_n = org_n + inorg_n + n_dep + n_min - n_uptake - n_denit
    avail_n[avail_n<0] = 0
    
    return avail_n

def est_temp_factor(t_min, t_max):
    """ Estimate temperature factor (linked to soil temperature) from air 
        temperature. Based on equations from the STREAM-N model. 
    
    Args:
        t_min:  Grid of minimum daily temperature.
        t_max:  Grid of maximum daily temperature.
    
    Returns:    
        Temperature factor (float).
    """
    t_soil = (0.57*(t_min + t_max)/2.0) + 3.3
    t_fact = 1.047**(t_soil - 20)

    return t_fact

def est_moisture_fact(wat_lev, fc):
    """ Estimate moisture factor from field capacity and water level. Based on 
        equations from the STREAM-N model. 
    
    Args:
        wat_lev: Water level grid
        fc:      Soil field capacity grid.
    
    Returns:    
        Moisture factor (float).
    """
    moist_ratio = wat_lev/fc
    moist_fact = moist_ratio.copy()
    moist_fact[moist_ratio>=0.95] = 1
    moist_fact[(0.3<=moist_ratio) & (moist_ratio<0.95)] = (
                moist_ratio[(0.3<=moist_ratio) & (moist_ratio<0.95)] - 0.3)/0.65
    moist_fact[moist_ratio<0.3] = 0
    moist_fact.mask = moist_ratio.mask

    return moist_fact

def est_mineralisation(calib_min, t_fact, moist_fact):
    """ Calculate mineralisation grid.
    
    Args:
        calib_min:  Calibration parameter determining rate of N mineralisation.
        t_fact:     Temperature factor from est_temp_factor().
        moist_fact: Moisture factor from est_moisture_factor().
        
    Returns:
        Grid of mineralisation rates.
    """
    return calib_min*t_fact*moist_fact

def est_denitrification(calib_denit, wat_lev, t_fact, moist_fact, avail_n):
    """ Calculate denitrification grid.
    
    Args:
        calib_denit: Calibration parameter determining rate of N denitrification.
        wat_lev:     Water level grid
        t_fact:      Temperature factor from est_temp_factor().
        moist_fact:  Moisture factor from est_moisture_factor().
        avail_n:     Grid of N available for leaching in each cell (kg/ha).
    
    Returns:
        Grid of denitrification rates.
    """
    # Calculate N conc. in cell from avail_n in kg/Ha and wat_lev in mm
    n_conc = (avail_n*1000.) / ((wat_lev/1000.)*100.*100.)  # in g/m3

    # Estimate denitrification rate
    return calib_denit*t_fact*moist_fact*n_conc  # in kg/Ha? Check units!

def estimate_n_added(annual_n_dict, daily_n_dep, org_n_fact, n_mineral,
                     n_denit, ts_row):
    """ Calculate amount of N added to each grid cell during this time step.
    
    Args:
        annual_n_dict: Dict returned by input_output.read_annual_n_grids().
        daily_n_dep:   Daily N deposition grid.
        org_n_fact:    Calibrating parameter determining the fraction of 
                       organic N that becomes immediately available for 
                       leaching.
        n_mineral:     Mineralisation grid from est_mineralisation().
        n_denit:       Denitrification grid from est_denitrification().
        ts_row:        Array object returned by input_output.read_ts_table().
        
    Returns:
        Array (grid of N added during this step).
    """
    # Calculate the amount of N added at this time step. The equation is
    # n_added = org_n + inorg_n + n_dep + n_min - n_uptake - n_denit
    # Each of these terms involves many other terms, so the equation is long, but
    # simple
    n_added = (org_n_fact*annual_n_dict[('or', 'gr')]*ts_row['or_gr'] +
               org_n_fact*annual_n_dict[('or', 'sp')]*ts_row['or_sp'] +
               org_n_fact*annual_n_dict[('or', 'wi')]*ts_row['or_wi'] +

               annual_n_dict[('in', 'gr')]*ts_row['in_gr'] +
               annual_n_dict[('in', 'sp')]*ts_row['in_sp'] +
               annual_n_dict[('in', 'wi')]*ts_row['in_wi'] -

               annual_n_dict[('up', 'gr')]*ts_row['up_gr'] -
               annual_n_dict[('up', 'sp')]*ts_row['up_sp'] -
               annual_n_dict[('up', 'wi')]*ts_row['up_wi'] -
               annual_n_dict[('up', 'ot')]*ts_row['up_ot'] +

               daily_n_dep +
               n_mineral -
               n_denit)

    return n_added

def calculate_n_conc(leached_n, tot_dr):
    """ Calculate N concentration in mg/l from leached N in kg/ha and total 
        drainage in mm.
    
    Args:
        leached_n: Grid of leached N in kg/ha.
        tot_dr:    Grid of total drainage in mm.
        
    Returns:
        Grid of N concentration in mg/l.
    """
    return (leached_n*1.0E6)/(100.*100.*tot_dr)

def calculate_n_leaching(avail_n, n_added, dr_list, fc, calib_n_leach):
    """ Estimate the amount of N leached by the draining water during this
        time step.
    
    Args:
        avail_n:       Grid of available N at end of previous step (kg/ha). 
        n_added:       N added during this step (kg/ha).  
        dr_list:       List of drainage grids from drainage.estimate_drainage().
        fc:            Soil field capacity grid.
        calib_n_leach: Calibrating parameter determining the leaching rate.
    
    Returns:
        List of grids: [leached_n, new_avail_n].
    """
    import numpy.ma as ma

    # Extract individual parameters from drainage_list
    snow_pk, wat_lev, surf_ro, lat_dr, vert_dr, tot_dr = dr_list

    # Add today's N to the existing amount in the soil
    intermed_n = avail_n + n_added
    intermed_n[intermed_n<0]=0

    # Calculate the fraction of the available N that leaches
    exponent = -1*calib_n_leach*(lat_dr+vert_dr)/fc
    leach_frac = 1 - ma.exp(exponent)

    # Calculate the amount of N that leaches
    leached_n = leach_frac*intermed_n

    # Calculate amount of N left
    new_avail_n = intermed_n - leached_n

    return [leached_n, new_avail_n]