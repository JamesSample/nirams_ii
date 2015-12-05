#-------------------------------------------------------------------------------
# Name:        drainage.py
# Purpose:     Calculate drainage for NIRAMS II.
#
# Author:      James Sample
#
# Created:     19/01/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" This module takes the pre-existing soil water level and snow pack grids,
    together with today's meteorological data, and estimates the drainage
    parameters.
"""

def estimate_drainage(fc, sat, calibl, caliv, wat_lev, snow_pk, rain, snow, 
                      melt, aet):
    """ Take the met data for the present step, plus the soils water level & 
        snowpack grids from the end of the previous time step and calculate 
        today's drainage.
    
    Args:
        fc:      Soil field capacity grid.
        sat:     Soil saturation capacity grid.
        calibl:  Lateral flow calibration grid.
        calibv:  Vertical flow calibration grid.
        wat_lev: Water level grid for end of previous time step.
        snow_pk: Snowpack thickness grid from end of previous time step.
        rain:    Rainfall grid for current time step.
        snow:    Snowfall grid for current time step.
        melt:    Snow melt grid for current time step.
        aet:     Actual evapotranspiration grid for current time step.
    
    Returns:
        List of arrays: [new_snow_pk, new_wat_lev, surface_runoff, 
                         lateral_drainage, vertical_drainage, total_drainage]
    """
    # Add today's water to the existing store
    intermed_wat_lev = wat_lev + rain + melt - aet
    intermed_wat_lev[intermed_wat_lev<0] = 0

    # Add today's snow to the existing snow pack
    new_snow_pk = snow_pk + snow - melt
    new_snow_pk[new_snow_pk<0] = 0

    # Surface runoff
    surf_ro = intermed_wat_lev - sat
    surf_ro[surf_ro<0] = 0

    # Available drainage
    avail_drain = intermed_wat_lev - surf_ro - fc
    avail_drain[avail_drain<0] = 0

    # Partition available drainage
    lat_drain = calibl*avail_drain
    vert_drain = caliv*avail_drain

    # Total drainage
    tot_drain = surf_ro + lat_drain + vert_drain

    # New water level
    new_wat_lev = intermed_wat_lev - tot_drain
    new_wat_lev[new_wat_lev<0] = 0

    return [new_snow_pk, new_wat_lev, surf_ro, lat_drain, vert_drain, tot_drain]

def reduce_aet_if_dry(aet, wat_lev, fc):
    """ Reduce actual evapotranspiration if the soil is dry. If the water level 
        in a cell is less than 0.7*fc, the rate of evapo-transpiration is 
        reduced by a factor. This factor is 1 when wat_lev = 0.7*fc and 
        decreases linearly to reach 0 when wat_lev = 0 i.e. where 
        wat_lev < 0.7*fc, apply a correction factor of
        wat_lev/(0.7*fc) to the aet grid.
    
    Args:
        aet:     "Raw" actual evapotranspiration grid.     
        wat_lev: Water level grid
        fc:      Soil field capacity grid.
        
    Returns:
        Array (modified AET grid with AET reduced where necessary).
    """
    # Get a boolean array showing which cells need correcting
    bool_array = wat_lev < (0.7*fc)

    # Calculate a correction factor for all cells, but subtract 1 from answer
    cor_facts_minus1 = (wat_lev / (0.7*fc)) - 1

    # Multiplying bool_array by cor_facts_minus1 gives a grid with values of
    # (cor_fact - 1) for cells which need correcting and zero otherwise. Add 1
    # to this to get a grid of cor_facts and ones
    cor_facts = (bool_array * cor_facts_minus1) + 1

    return aet*cor_facts