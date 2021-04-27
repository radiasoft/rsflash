#Functions for data extraction from FLASH yt files to support plotting
#No specific plotting features implemented here.
#May be scrapped or replaced!

# Stephen Coleman and Nathan Cook
# Most recent updates - 03/01/2021

import yt
import numpy as np
from itertools import groupby #used for interpolation

def get_lineouts(files, field, axis, Nmax, interpolate=True, interpolate_max=2000, coordinates=[(0,0,0),(0,0,0)], xyz=False, cubic=False):
    '''
    Leverage extract_line on multiple files to grab from multiple times

    Parameters: Required
    ----------

      files: list
         A list of filenames referencing FLASH datasets
      Nmax: integer
         The number of points to sample
      axis: string
         The axis ('r','z') in polar coordinates that you're sampling along
      field_name: string
         Which field from the file the user wishes to plot


    Options: Keyword Arguments
     interpolate: bool (Default[True])
          If the user wants a smooth rather than staircased sequence of values
     step_threshold: float (Default[1e-6])
          The maximum fractional difference between two
          step values below, which we assume the differences are
          just float rounding error
     interpolate_max: double (Default[2000])
          The maximum distance along the sampled axiswhere we try to use
          interpolation
     min_levels: int (Default[3])
          The minimum number of steps we'll accept for an interpolated plot.
          If the number of steps is below this, we return the un-interpolated plot
          even if interpolate=True
     coordinates: [(x1,y1,z1),(x2,y2,z2)] ( Default is [(0,0,0),(0,0,0)] )
          Sample along an arbitrary line for a lineout from sets of coorinates 1 to 2
    xyz: bool (Default[False])
          Return the sampled points in the flash space as ordered lists of
          coordinates. This changes the nature of the output.
          Be careful with outputs if mixing this with interpolation.
    cubic: bool (Deafult[False])
         Use Scipy UnivariateSpline routine to further smooth the plot output
         over 500 interpolated points. 


    '''

    times = []
    xlineouts = []
    ylineouts = []

    xcvals = []
    ycvals = []
    zcvals = []

    for file in files:

        ds = yt.load(file)
        time = ds.parameters['time']
        
        if xyz:
            xc, yc, zc, linex, liney = extract_line(ds, field, axis, Nmax, interp=interpolate, interpolate_max=interpolate_max, coordinates=coordinates, xyz=xyz, cubic=cubic)
            
            xcvals.append(xc)
            ycvals.append(yc)
            zcvals.append(zc)
            xlineouts.append(linex)
            ylineouts.append(liney)
            times.append(time)
            
            return xcvals, ycvals, zcvals, xlineouts, ylineouts, times
            
        else:
            xvals, yvals = extract_line(ds, field, axis, Nmax, interp=interpolate, interpolate_max=interpolate_max, coordinates=coordinates, xyz=xyz, cubic=cubic)
            
            xlineouts.append(xvals)
            ylineouts.append(yvals)
            times.append(time)
            
            return xlineouts, ylineouts, times

def interpolate_lineout(xarr,yarr,**kwargs):
    '''
    Convert sampled points of a field on a mesh into smooth lines for plotting

    Parameters: Required
    ----------

      xarr: list
         The independent variable values
      yarr: list
         The dependent variable values


    Options: Keyword Arguments
     step_threshold: float (Default[1e-6])
          The maximum fractional difference between two
          step values below, which we assume the differences are
          just float rounding error
     interpolate_max: double (Default[2000])
          The maximum distance along the sampled axiswhere we try to use
          interpolation
     min_levels: int (Default[3])
          The minimum number of steps we'll accept for an interpolated plot.
          If the number of steps is below this, we return the un-interpolated plot
          even if interpolate=True
     cubic: bool (Deafult[False])
          Use Scipy UnivariateSpline routine to further smooth the plot output
          over 500 interpolated points.
     smoothing: float (Default[1.0])
          Tune the s value for Scipy Univariate spline function. See
          https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
    '''

    plateaus = [list(g) for k,g in groupby(yarr)] #list of lists of identical values in a row
    new_yvals = [i[0] for i in plateaus]          #List of all unique values, ordered
    new_xvals = []                                #Where we'll store the new x values
    xarr_consumable = xarr.copy()                 #We'll mutate this list

    #Count the plateaus - if there aren't enough levels
    # represented then the interpolation won't look great.
    # Check against a minimum number of levels, and if the
    # threshold isn't met then kick out the original lineout

    # Sometimes levels are separated by < 0.000001%. We're going
    # to call that float rounding error and lump them together.
    new_plateaus = [plateaus[0]]

    for i in range(1,len(plateaus)):
        diff = plateaus[i][0]/plateaus[i-1][0]
        if np.abs(1.0 - diff) < kwargs.pop("step_threshold",1e-6):
            new_plateaus[-1].append(plateaus[i])
        else:
            new_plateaus.append(plateaus[i])
    plateaus = new_plateaus

    if len(plateaus) < kwargs.pop("min_levels",3):
        return xarr,yarr

    interp_max = kwargs.pop("interpolate_max",2000.)

    #Find the central values of the plateaus
    for j,jj in enumerate(plateaus):                   #Working through the steps. jj is repeating list of yvals
        if xarr_consumable[0] < interp_max:            #Make sure we're below cutoff
            new_xvals.append(xarr_consumable[0:len(jj)].mean())     #Store center of the plateau 'bin'
            xarr_consumable = xarr_consumable[len(jj):]             #Trim off previous 'bin', prepare for the next 'bin'
        else:                                          #We're above threshold, stop
            #flatten the remaining list of lists into a list
            remaining_values = plateaus[j:]
            remaining_values = [item for sublist in remaining_values for item in sublist] #list comprehension magic
            #Append the remaining data point-by-point
            new_yvals = new_yvals[:j] + remaining_values     #Join the unique y values with those above threshold
            new_xvals = new_xvals + xarr_consumable.tolist() #same for x

            break #Exit the loop, no more plateau values to find

    #Regardless of how an arbitrary sampling line is defined, we need
    # the x values (and y values) sorted together based on ascending x values
    zipped_lists = zip(new_xvals, new_yvals)
    sorted_pairs = sorted(zipped_lists)
    new_xvals, new_yvals = [ list(tuple) for tuple in zip(*sorted_pairs)]

    #Use some splining. May need to allow for tuning of smoothing factor s
    if kwargs.get("cubic",False):
        from scipy.interpolate import UnivariateSpline
        print(new_xvals)
        f2 = UnivariateSpline(new_xvals, new_yvals, s = kwargs.get("s",1.0))
        xnew = np.linspace(new_xvals[0], new_xvals[-1], num=501, endpoint=True)
        new_xvals = xnew
        new_yvals = f2(xnew)

    return new_xvals,new_yvals

def extract_line(data_set, field_name, axis, Nmax, interp=True, **kwargs):
    """
    Sample points along an arbitrary 1-d line in a FLASH dataset using yt

    Parameters: Required
    ----------

      data_set: A FLASH Dataset
         The FLASH dataset, preloaded via yt.load
      Nmax: integer
         The number of points to sample
      axis: string
         The axis ('r','z') in polar coordinates that you're sampling along
      field_name: string
         Which field from the file the user wishes to plot


    Options: Keyword Arguments
     interpo: bool (Default[True])
          If the user wants a smooth rather than staircased plot
     step_threshold: float (Default[1e-6])
          The maximum fractional difference between two
          step values below, which we assume the differences are
          just float rounding error
     interpolate_max: double (Default[2000])
          The maximum distance along the sampled axiswhere we try to use
          interpolation
     min_levels: int (Default[3])
          The minimum number of steps we'll accept for an interpolated plot.
          If the number of steps is below this, we return the un-interpolated plot
          even if interpolate=True
     coordinates: [(x1,y1,z1),(x2,y2,z2)] ( Default is [(0,0,0),(0,0,0)] )
          Sample along an arbitrary line for a lineout from sets of coorinates 1 to 2
          For cylindrical coordinates, this is [(r1,z1,theta1),(r2,z2,theta2)]
     xyz: bool (Default[False])
          Return the sampled points in the flash space as ordered lists of
          coordinates. This changes the nature of the output.
          Be careful with outputs if mixing this with interpolation.

    Outputs:
    -----------

    Default:
        xarr,yarr  (both numpy arrays)
            xarr are the positions along the sampled line, beginning at 0
            yarr are the selected field values at those positions
    if xyz is True:
        xcoords, ycoords, zcoords, xarr, yarr
            (xyz)coords are the orginial sampled coordinates in 3d space
            xarr are the magnitudes of the positions along the sampled line, beginning at 0
            yarr are the selected field values at those positions
"""

    fields = data_set.field_list

    if ('flash',field_name) not in fields:
        print("Field name "+field_name+" not in file. Available fields:")
        print(fields)
        return

    yarr = np.zeros(Nmax)

    #Handle the arbitrary option by checking if the default
    # start and end coordinates are the same. If they are not,
    # then recalculate the limits and overwrite the 'axis'
    # variable to 'arb' to catch that case
    coords1, coords2 = kwargs.pop("coordinates",[(0,0,0),(0,0,0)])

    if data_set.geometry == "cylindrical":
        #The limits for a cylindrical simulation, based on the file
        xmin = 0.0; xmax = data_set.domain_width[0] + xmin # the max range of r
        zmin = data_set.index.grid_left_edge[0][1];        # the min value of z
        zmax = data_set.domain_width[1] + zmin             # the max value of z

        #Check if the arbitrarily defined plot exceeds these limits
        if coords1 != coords2:
            xmin_arb = coords1[0]; xmax_arb = coords2[0]
            ymin_arb = coords1[1]; ymax_arb = coords2[1]
            zmin_arb = coords1[2]; zmax_arb = coords2[2]
            axis='arb'
            print("Sampling between coordinate pair")

            if( (xmin_arb < xmin) | (xmax_arb < xmin) |
                (xmin_arb > xmax) | (xmax_arb > xmax) |
                (ymin_arb < zmin) | (ymax_arb < zmin) |
                (ymin_arb > zmax) | (ymax_arb > zmax) ):
                    print("Plot axis out of bounds")
                    return


        if str(data_set.dimensionality) == '3':
            print("Warning: All extracted lines at theta = pi")

        if axis.lower() == 'r':
            x_coords = [xmin + (0.5 + j) * (xmax - xmin) / Nmax for j in range(Nmax)]
            y_coords = [0.5 * (zmin + zmax) for j in range(Nmax)]
            z_coords = [np.pi for j in range(Nmax)]
            xarr = np.array(x_coords)

        elif axis.lower() == 'z':
            x_coords = [xmin + (0.5) * (xmax - xmin) for j in range(Nmax)]
            y_coords = [zmin + (0.5 + j) * (zmax - zmin) / Nmax for j in range(Nmax)]
            z_coords = [np.pi for j in range(Nmax)]
            xarr = np.array(y_coords)

        elif axis.lower() == "arb":
            x_coords = [xmin_arb + (0.5 + j) * (xmax_arb - xmin_arb) / Nmax for j in range(Nmax)]
            y_coords = [zmin_arb + (0.5 + j) * (zmax_arb - zmin_arb) / Nmax for j in range(Nmax)]
            z_coords = [np.pi for j in range(Nmax)]
            xarr = np.sqrt(np.add(np.power(x_coords,2),np.power(y_coords,2)))

        else:
            print("Select 'r' or 'z' for cylindrical coordinate system")

    if data_set.geometry == "cartesian":
        #The limits for a cartesian simulation
        xmin = data_set.index.grid_left_edge[0][0]; xmax = data_set.domain_width[0] + xmin
        ymin = data_set.index.grid_left_edge[0][1]; ymax = data_set.domain_width[1] + ymin
        zmin = data_set.index.grid_left_edge[0][2]; zmax = data_set.domain_width[2] + zmin

        #Check if the arbitrarily defined plot exceeds these limits
        # BTW it's not really 'xmin' or 'xmax', more like x1 and x2
        if coords1 != coords2:
            xmin_arb = coords1[0]; xmax_arb = coords2[0]
            ymin_arb = coords1[1]; ymax_arb = coords2[1]
            zmin_arb = coords1[2]; zmax_arb = coords2[2]
            axis='arb'
            print("Sampling between coordinate pair")

            if( (xmin_arb < xmin) | (xmax_arb < xmin) |
                (xmin_arb > xmax) | (xmax_arb > xmax) |
                (ymin_arb < ymin) | (ymax_arb < ymin) |
                (ymin_arb > ymax) | (ymax_arb > ymax) |
                (zmin_arb < zmin) | (zmax_arb < zmin) |
                (zmin_arb > zmax) | (zmax_arb > zmax) ):
                    print("Plot axis out of bounds")
                    return

        if axis.lower() == 'x':
            x_coords = [xmin + (0.5 + j) * np.abs(xmax - xmin) / Nmax for j in range(Nmax)]
            y_coords = [ymin + 0.5 * (ymin + ymax) for j in range(Nmax)]
            z_coords = [zmin + 0.5 * (zmin + zmax) for j in range(Nmax)]
            xarr = np.array(x_coords)

        elif axis.lower() == 'y':
            x_coords = [xmin + (0.5) * (xmax + xmin) for j in range(Nmax)]
            y_coords = [ymin + (0.5 + j) * (ymax - ymin) / Nmax for j in range(Nmax)]
            z_coords = [zmin + (0.5) * (zmax + zmin) for j in range(Nmax)]
            xarr = np.array(y_coords)

        elif axis.lower() == 'z':
            x_coords = [xmin + (0.5) * (xmax + xmin) for j in range(Nmax)]
            y_coords = [ymin + (0.5) * (ymax + ymin) for j in range(Nmax)]
            z_coords = [zmin + (0.5 + j) * np.abs(zmax - zmin) / Nmax for j in range(Nmax)]
            xarr = np.array(y_coords)

        elif axis.lower() == "arb":
            x_coords = [xmin_arb + (0.5 + j) * np.abs(xmax_arb - xmin_arb) / Nmax for j in range(Nmax)]
            y_coords = [ymin_arb + (0.5 + j) * np.abs(ymax_arb - ymin_arb) / Nmax for j in range(Nmax)]
            z_coords = [zmin_arb + (0.5 + j) * np.abs(zmax_arb - zmin_arb) / Nmax for j in range(Nmax)]
            xarr = np.sqrt( np.power(np.array(x_coords)-xmin_arb,2) + 
                            np.power(np.array(y_coords)-ymin_arb,2) + 
                            np.power(np.array(z_coords)-zmin_arb,2) )

        else:
            print("Select 'x' 'y' or 'z' for cartesian coordinate system")


    #Sample from the dataset
    for j in range(Nmax):
        #convert kelvin to eV and drop units
        DataPoint = data_set.point([x_coords[j],y_coords[j],z_coords[j]])
        yarr[j] = DataPoint['flash', field_name].in_cgs().d[0]

    #Fix 'staircasing' by connecting dots between the
    # central values of each step.  This necessarily breaks the
    # relationship between the sampled field and the (x,y,z) coordinates
    # stored above, since we're determining values along a new
    # line r' with variable spacing.
    if interp:
        new_xvals, new_yvals = interpolate_lineout(xarr,yarr,**kwargs)

        if kwargs.pop("xyz",False):
            return np.asarray(x_coords), np.asarray(y_coords), np.asarray(z_coords), np.asarray(new_xvals), np.asarray(new_yvals)
        else:
            return np.asarray(new_xvals), np.asarray(new_yvals)


    else:  #If we didn't interpolate, return the raw values
        if kwargs.pop("xyz",False):
            return np.asarray(x_coords), np.asarray(y_coords), np.asarray(z_coords), np.asarray(xarr), np.asarray(yarr)
        else:
            return np.asarray(xarr),np.asarray(yarr)
