#-------------------------------------------------------------------------------
# Name:        snow.py
# Purpose:     Snowfall/melt module for NIRAMS II.
#
# Author:      James Sample
#
# Created:     18/01/2012
# Copyright:   (c) James Sample and JHI, 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
""" This module contains functions to determine the amount of snow and rain
    each day, based on the total precipitation observed and the max and min
    temperatures. Also estimates snow melt each day using a simple degree-day
    melting model.
"""

def estimate_snow_and_rain(pptn, t_min, t_max, t_snow):
    """ Split precipitation grid into snowfall and rainfall components
        according to temperature. Assumes that temperature increases linearly 
        from t_min to t_max and then decreases linearly again back to t_min in 
        each 24hr period (i.e. a symmetrical triangular temperature profile). 
        In this simplistic case, the proportion of the day above t_snow is:

            (t_max - t_snow)/(t_max - t_min)
    
    Args:
        pptn:   Grid of total daily precipitation.
        t_min:  Grid of minimum daily temperature.
        t_max:  Grid of maximum daily temperature.
        t_snow: Threshold temperature below which precipitation falls as snow
                (float).
    
    Returns:
        List of arrays: [rain, snow]
    """
    # Calculate fraction of rain. There are three cases to consider:
    #   1. t_max > t_snow and t_min < t_snow. Use equation above.
    #   2. t_max <= t_snow. fr_rn = 0
    #   3. t_min > t_snow. fr_rn = 1
    # Calculate case 1
    fr_rn = (t_max-t_snow)/(t_max-t_min)

    # Deal with case 2
    fr_rn[t_max<=t_snow] = 0

    # Deal with case 3
    fr_rn[t_min>t_snow] = 1

    # Calculate rain and snow
    rain = pptn*fr_rn
    snow = pptn*(1-fr_rn)

    return [rain, snow]

def estimate_snow_melt(snow_pk, t_min, t_max, t_melt, ddf):
    """ Calculte the amount of melt from the snowpack according to temperature.
        As above, assumes temperature increases linearly from t_min to t_max
        and then decreases linearly back to t_min in each 24hr period. The 
        proportion of the day above t_melt is:

            (t_max - t_melt)/(t_max - t_min)
            
        The melting model assumes that         
            
            melt =  k(t_av - t_melt)      for t_av > t_melt
            melt = 0                      for t_av < t_melt

        where melt in mm, k is the degree-day factor (in mm/C/d) and t_melt
        is the threshold temperature for melting. t_av is the average of 
        temperatures above t_melt i.e:

            t_av = (t_max + t_melt)/2          for t_max > t_melt
    
    Args:
        snow_pk: Grid of snow pack thickness.
        t_min:   Grid of minimum daily temperature.
        t_max:   Grid of maximum daily temperature.
        t_melt:  Threshold temperature above which melting begins (float).
        ddf:     Degree-day factor (float)

    Returns:
        Array of melt.
    """
    # There are three cases to consider:
    #   1. t_min >= t_melt:
    #         m = k(((t_max+t_min)/2) - t_melt)
    #   2. t_max >= t_melt and t_min < t_melt:
    #         m = kh(((t_max+t_melt)/2) - t_melt)
    #      where h = (t_max - t_melt)/(t_max - t_min)
    #   3. t_max < t_melt:
    #         melt = 0
    # Deal with case 1
    melt1 = ddf*(((t_max+t_min)/2) - t_melt)
    melt1[t_min<t_melt] = 0     # Case 1 doesn't apply where t_min<t_melt

    # Deal with case 2
    fr_above_t_melt = (t_max - t_melt)/(t_max - t_min)
    melt2 = ddf*fr_above_t_melt*(((t_max+t_melt)/2) - t_melt)
    melt2[t_min>=t_melt] = 0    # Case 2 doesn't apply where t_min>=t_melt

    # Deal with case 3
    # The first of these shouldn't be necessary, but the Met office grids have
    # some cells where t_min > t_max!
    melt1[t_max<t_melt] = 0     # Case 1 doesn't apply where t_max<t_melt
    melt2[t_max<t_melt] = 0     # Case 2 doesn't apply where t_max<t_melt

    # Calculate melt
    melt = melt1+melt2

    # The total amount of melt can not exceed the amount of snow available.
    # Set melt values greater than snow available back to max possible
    melt[melt>snow_pk] = snow_pk[melt>snow_pk]

    return melt