# 2D slice plot module for rsflash based upon yt implementation.
# Currently a work in progress. May be enhanced or scrapped!
#
# Nathan Cook and Stephen Coleman
# 01/07/2021
#

import yt


def slice_plot_2D(ds, field, axis, **kwargs):
    """
    Generate and save a 2D slice plot of a FLASH yt dataset.

    Parameters: Required
    ----------

        ds: a FLASHDataset
            The dataset, loaded via yt.load
        field: string
            The field being plotted
        axis: string
            The axis ('x', 'y', 'z' for Cartesian) about which the plot is made

    Options: Keyword Arguments
    --------
        title: string (Default[None])
            The plot title. Defaults to stating only the field.
        clabel: string (Default[None])
            The color bar label. Defaults to stating only the field.
        basename: string (Default[ds.basename])
            Prefix used for saving and other bookeeping
        zlim: two-tuple (Default[None])
            Specify color bar limit
        set_log: boolean (Default[False])
            Specify log-scale of plotted variable
        save_loc: string (Default['./'])
            Subdirectory in which to save the file.
        save_index: boolean (Default[True])
            If true, label file with last 4 digits of the filename.
            Otherwise, uses step number.
        dosave: boolean (Default[True])
            If true, save the file. Otherwise don't.
    """

    # read and update keyword arguments
    options = {
        "title": field,
        "clabel": field,
        "basename": ds.basename,
        "zlim": None,
        "set_log": False,
        "save_loc": "./",
        "dosave": True,
        "save_index": True,
    }

    options.update(kwargs)

    # define other metadata
    time = ds.parameters["time"]
    nstep = ds.parameters["nstep"]

    # create slice plot
    slc = yt.SlicePlot(ds, axis, field, origin="native", aspect=1)

    # adjust sizes and annotations
    slc.annotate_grids(edgecolors="white")
    slc.annotate_title(options["title"])
    slc.figure_size = 12
    slc.set_font_size(26)
    slc.set_cmap(field, cmap="viridis")

    # colorbar manipulations
    slc.set_colorbar_label(field, options["clabel"])
    if options["zlim"] is not None:
        slc.set_zlim(field, zmin=options["zlim"][0], zmax=options["zlim"][1])

    slc.set_log(field, options["set_log"])

    # save options
    if options["save_index"]:
        index = ds.basename[-4:]
        savename = "{}{}_{}_{}.png".format(
            options["save_loc"], options["basename"], field, index
        )
    else:
        savename = "{}_{}_{}.png".format(
            options["save_loc"], options["basename"], field, nstep
        )
    if options["dosave"]:
        slc.save(savename)

    return slc
