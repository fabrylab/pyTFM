'''
plotting and other accesory functions for the Monolayer stress Microscopy

'''
from contextlib import suppress

import matplotlib
from matplotlib import cm  # making list from colormap
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# matplotlib.use("TkAgg")
# matplotlib.use("Agg") # when running the script to avoid plots popping up
from pyTFM.TFM_functions import *
# from matplotlib.figure import Figure
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from scipy.ndimage import binary_erosion
from tqdm import tqdm


def make_discrete_colorbar():
    '''
    function to make a three pieced colormap, only used for represetatinof loads in labt atlk example.
    :param values:
    :return:
    '''

    cmap = plt.cm.jet  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey

    # create the new map
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize

    norm = matplotlib.colors.BoundaryNorm([-1, -0.5, 0.5, 1], cmap.N)
    return norm


def show_grid(ax):
    ax.grid(True, color="black", alpha=0.3)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))


def hide_ticks(ax, interval):
    for i, label in enumerate(ax.xaxis.get_ticklabels()[:]):
        if (i % interval) == 0:
            continue
        else:
            label.set_visible(False)
    for i, label in enumerate(ax.yaxis.get_ticklabels()[:]):
        if (i % interval) == 0:
            continue
        else:
            label.set_visible(False)


def plot_continuous_boundary_stresses(plot_values, mask_boundaries=None, plot_t_vecs=False, plot_n_arrows=False,
                                      figsize=(10, 7),
                                      scale_ratio=0.2, border_arrow_filter=1, cbar_str="line tension in N/m", vmin=None,
                                      vmax=None,
                                      cbar_width="2%", cbar_height="50%", cbar_axes_fraction=0.2,
                                      cbar_tick_label_size=20,
                                      background_color="white", cbar_borderpad=0.1, linewidth=4, cmap="jet",
                                      plot_cbar=True, cbar_style="clickpoints",
                                      boundary_resolution=3, cbar_title_pad=1, outer_cb_color="grey",
                                      outer_cb_style="-",
                                      **kwargs):
    '''
    plotting the line stresses (total transmitted force of cell boundaries), colored by their absolute values
    as continous lines.
    :param shape:
    :param edge_lines:
    :param lines_interpol:
    :param min_v:
    :param max_v:
    :param mask_boundaries:
    :param plot_t_vecs:
    :param plot_n_arrows:
    :param figsize:
    :param scale_ratio:
    :param arrow_filter:
    :param cbar_str:
    :param vmin:  overwrites max_v and min_v if provided
    :param vmax:  overwrites max_v and min_v if provided
    :return:
    '''

    if not isinstance(plot_values[0], (list)):
        plot_values = [plot_values]

    min_v = np.min([pv[3] for pv in plot_values])  # minimum over all objects
    max_v = np.max([pv[4] for pv in plot_values])  # maximum over all objects
    shape = plot_values[0][0]  # image shape, should be the same for all objects
    print("plotting cell border stresses")
    min_v = vmin if isinstance(vmin, (float, int)) else min_v
    max_v = vmax if isinstance(vmax, (float, int)) else max_v
    mask_boundaries = np.zeros(shape) if not isinstance(mask_boundaries, np.ndarray) else mask_boundaries

    fig = plt.figure(figsize=figsize)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    background_color = matplotlib.cm.get_cmap(cmap)(0) if background_color == "cmap_0" else background_color
    fig.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_axis_off()
    fig.add_axes(ax)

    for shape, edge_lines, lines_interpol, *rest in plot_values:
        all_t_vecs = np.vstack([subdict["t_vecs"] for subdict in lines_interpol.values()])
        if plot_t_vecs:
            scale = scale_for_quiver(all_t_vecs[:, 0], all_t_vecs[:, 1], dims=mask_boundaries.shape,
                                     scale_ratio=scale_ratio, return_scale=True)
        for line_id, interp in tqdm(lines_interpol.items(), total=len(lines_interpol.values())):
            p_new = interp["points_new"]
            x_new = p_new[:, 0]
            y_new = p_new[:, 1]
            t_norm = interp["t_norm"]
            t_vecs = interp["t_vecs"]
            n_vecs = interp["n_vecs"]
            # plotting line segments

            c = matplotlib.cm.get_cmap(cmap)(
                (t_norm - min_v) / (max_v - min_v))  # normalization and creating a color range
            ## see how well that works
            if line_id in edge_lines:  # plot lines at the edge
                plt.plot(x_new, y_new, outer_cb_style, color=outer_cb_color, linewidth=linewidth)
            else:
                for i in range(0, len(x_new) - boundary_resolution, boundary_resolution):
                    plt.plot([x_new[i], x_new[i + boundary_resolution]], [y_new[i], y_new[i + boundary_resolution]],
                             color=c[i], linewidth=linewidth)

            # plotting stress vectors
            if plot_t_vecs:
                t_vecs_scale = t_vecs * scale
                for i, (xn, yn, t) in enumerate(zip(x_new, y_new, t_vecs_scale)):
                    if i % border_arrow_filter == 0:
                        plt.arrow(xn, yn, t[0], t[1], head_width=0.5)
            # plotting normal vectors
            if plot_n_arrows:
                for i in range(len(x_new) - 1):
                    if i % border_arrow_filter == 0:
                        plt.arrow(x_new[i], y_new[i], n_vecs[i][0], n_vecs[i][1], head_width=0.5)

    plt.gca().invert_yaxis()  # to get the typicall imshow orientation
    plt.xlim(0, shape[1])
    plt.ylim(shape[0], 0)
    # background_color=matplotlib.cm.get_cmap(cmap)(0) if background_color=="cmap_0" else background_color
    # ax.set_facecolor(background_color)
    if plot_cbar:
        add_colorbar(min_v, max_v, cmap, ax=ax, cbar_style=cbar_style, cbar_width=cbar_width, cbar_height=cbar_height,
                     cbar_borderpad=cbar_borderpad, v=cbar_tick_label_size, cbar_str=cbar_str,
                     cbar_axes_fraction=cbar_axes_fraction, cbar_title_pad=cbar_title_pad)
    return fig, ax


def add_colorbar(vmin, vmax, cmap="rainbow", ax=None, cbar_style="not-clickpoints", cbar_width="2%",
                 cbar_height="50%", cbar_borderpad=0.1, cbar_tick_label_size=15, cbar_str="",
                 cbar_axes_fraction=0.2, shrink=0.8, aspect=20, cbar_title_pad=1, **kwargs):
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=matplotlib.cm.get_cmap(cmap), norm=norm)
    sm.set_array([])  # bug fix for lower matplotlib version
    if cbar_style == "clickpoints":  # colorbar inside of the plot
        cbaxes = inset_axes(ax, width=cbar_width, height=cbar_height, loc=5, borderpad=cbar_borderpad * 30)
        cb0 = plt.colorbar(sm, cax=cbaxes)
        with suppress(TypeError, AttributeError):
            cbaxes.set_title(cbar_str, color="white", pad=cbar_title_pad)
        cbaxes.tick_params(colors="white", labelsize=cbar_tick_label_size)
    else:  # colorbar outide of the plot
        cb0 = plt.colorbar(sm, aspect=aspect, shrink=shrink, fraction=cbar_axes_fraction,
                           pad=cbar_borderpad)  # just exploiting the axis generation by a plt.colorbar
        cb0.outline.set_visible(False)
        cb0.ax.tick_params(labelsize=cbar_tick_label_size)
        with suppress(TypeError, AttributeError):
            cb0.ax.set_title(cbar_str, color="black", pad=cbar_title_pad)
    return cb0


def check_order(mask, coords):
    plt.figure()
    plt.imshow(mask, origin="lower")
    plt.plot(coords[:, 1], coords[:, 0], "o")
    for i, (x1, y1) in enumerate(zip(coords[:, 0], coords[:, 1])):
        plt.text(y1, x1, str(i))


def check_normal_vectors(mask, coords, n):
    plt.figure()
    plt.imshow(mask, origin="lower", cmap="viridis")
    plt.plot(coords[:, 1], coords[:, 0], "o")
    for i, (x1, y1) in enumerate(zip(coords[:, 1], coords[:, 0])):
        plt.text(x1, y1, str(i))
    for vec, p in zip(n, coords):
        plt.arrow(p[1], p[0], vec[0], vec[1], head_width=0.2)


def check_normal_vectors_graph(mask, n, points):
    '''
    n is dictionary point_id:normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    plt.figure()
    plt.imshow(mask, origin="lower", cmap="viridis")

    for p, vec in n.items():
        plt.arrow(points[p][1], points[p][0], vec[0], vec[1], head_width=0.2)


def check_normal_vectors_array(mask, n_array, origin="lower"):
    '''
    n is a 3d array.axis 0 and 1 represetn xy coordinates and axis 2 is the normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    plt.figure()
    plt.imshow(mask, origin=origin, cmap="viridis")
    coords = np.where(~np.logical_and(n_array[:, :, 0] == 0, n_array[:, :, 1] == 0))

    for px, py in zip(coords[0], coords[1]):
        plt.arrow(py, px, n_array[px, py, 0], n_array[px, py, 1], head_width=0.2)


def plot_stress_vectors(mask, arrows_array, origin="lower"):
    '''
    n is a 3d array.axis 0 and 1 represetn xy coordinates and axis 2 is the normal vector
    :param mask:
    :param n:
    :param points:
    :return:
    '''
    mask = mask.astype("bool")
    plt.figure()
    abs_im = np.linalg.norm(arrows_array, axis=2)
    abs_im[~mask] = np.nan
    im = plt.imshow(abs_im, origin=origin, cmap="viridis")
    plt.colorbar(im)

    coords = np.where(~np.logical_and(arrows_array[:, :, 0] == 0, arrows_array[:, :, 1] == 0))
    for px, py in zip(coords[0], coords[1]):
        plt.arrow(py, px, arrows_array[px, py, 0], arrows_array[px, py, 1], head_width=0.2)


def plot_nodes(nodes):
    plt.figure(figsize=(7, 7))
    plt.xlim(0, np.sqrt(len(nodes)))
    plt.ylim(0, np.sqrt(len(nodes)))
    for n in nodes:
        # plt.scatter(n[1],n[2],s=100,facecolors="none",edgecolors="black")
        plt.text(n[1], n[2], str(int(n[0])), bbox={"boxstyle": "circle", "edgecolor": "black", "facecolor": "none"},
                 color="black")


def plot_grid(nodes, elements, inverted_axis=False, symbol_size=4, arrows=False, image=0):  # a connecting lines

    plt.figure(figsize=(7, 7))

    if inverted_axis:
        nodes[:, [1, 2]] = nodes[:, [2, 1]]
        # plt.xlim(np.min(nodes[:, 1]), np.max(nodes[:, 1]))
        # plt.ylim(np.max(nodes[:, 1]), np.min(nodes[:, 2]))

    plt.xlim(np.min(nodes[:, 1]) - 10, np.max(nodes[:, 1]) + 10)
    plt.ylim(np.min(nodes[:, 1]) - 10, np.max(nodes[:, 2]) + 10)
    if inverted_axis:
        plt.gca().invert_yaxis()

    for n in nodes:
        # plt.scatter(n[1],n[2],s=100,facecolors="none",edgecolors="black")
        plt.text(n[1], n[2], str(int(n[0])), fontsize=symbol_size, alpha=1, ha="center", va="center",
                 bbox={"boxstyle": "circle", "edgecolor": "black", "facecolor": "white", "alpha": 0.7},
                 color="black")
    coms_x = (nodes[elements[:, 3], 1] + nodes[elements[:, 4], 1] + nodes[elements[:, 5], 1] + nodes[
        elements[:, 6], 1]) / 4
    coms_y = (nodes[elements[:, 3], 2] + nodes[elements[:, 4], 2] + nodes[elements[:, 5], 2] + nodes[
        elements[:, 6], 2]) / 4
    for i, com in enumerate(zip(coms_x, coms_y)):
        plt.text(com[0], com[1], str(i), fontsize=symbol_size, ha="center", va="center",
                 bbox={"boxstyle": "circle", "edgecolor": "blue", "facecolor": "white", "alpha": 0.7},
                 color="blue")
    if arrows:
        cols = ["black", "red", "green", "blue"]
        for e in elements:
            if e[0] % 10 == 0:
                col = cols[e[0] % 4]
                plt.arrow(nodes[e[3]][1], nodes[e[3]][2], nodes[e[4]][1] - nodes[e[3]][1],
                          nodes[e[4]][2] - nodes[e[3]][2], color=col, head_width=0.2)
                plt.arrow(nodes[e[4]][1], nodes[e[4]][2], nodes[e[5]][1] - nodes[e[4]][1],
                          nodes[e[5]][2] - nodes[e[4]][2], color=col, head_width=0.2)
                plt.arrow(nodes[e[5]][1], nodes[e[5]][2], nodes[e[6]][1] - nodes[e[5]][1],
                          nodes[e[6]][2] - nodes[e[5]][2], color=col, head_width=0.2)
                plt.arrow(nodes[e[6]][1], nodes[e[6]][2], nodes[e[3]][1] - nodes[e[6]][1],
                          nodes[e[3]][2] - nodes[e[6]][2], color=col, head_width=0.2)

    else:
        for e in elements:
            plt.plot([nodes[e[3]][1], nodes[e[4]][1]], [nodes[e[3]][2], nodes[e[4]][2]], color="black")
            plt.plot([nodes[e[4]][1], nodes[e[5]][1]], [nodes[e[4]][2], nodes[e[5]][2]], color="black")
            plt.plot([nodes[e[5]][1], nodes[e[6]][1]], [nodes[e[5]][2], nodes[e[6]][2]], color="black")
            plt.plot([nodes[e[6]][1], nodes[e[3]][1]], [nodes[e[6]][2], nodes[e[3]][2]], color="black")


def show_quiver(fx, fy, filter=[0, 1], scale_ratio=0.2, headwidth=None, headlength=None, headaxislength=None,
                width=None, cmap="rainbow",
                figsize=None, cbar_str="", ax=None, fig=None
                , vmin=None, vmax=None, cbar_axes_fraction=0.2, cbar_tick_label_size=15
                , cbar_width="2%", cbar_height="50%", cbar_borderpad=0.1,
                cbar_style="not-clickpoints", plot_style="not-clickpoints", cbar_title_pad=1, plot_cbar=True, alpha=1,
                ax_origin="upper", filter_method="regular", filter_radius=5, **kwargs):
    # list of all necessary quiver parameters
    quiver_parameters = {"headwidth": headwidth, "headlength": headlength, "headaxislength": headaxislength,
                         "width": width, "scale_units": "xy", "angles": "xy", "scale": None}
    quiver_parameters = {key: value for key, value in quiver_parameters.items() if not value is None}

    fx = fx.astype("float64")
    fy = fy.astype("float64")
    dims = fx.shape  # needed for scaling
    if not isinstance(ax, matplotlib.axes.Axes):
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
    map_values = np.sqrt(fx ** 2 + fy ** 2)
    vmin, vmax = set_vmin_vmax(map_values, vmin, vmax)
    im = plt.imshow(map_values, cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha, origin=ax_origin)  # imshowing
    if plot_style == "clickpoints":
        ax.set_position([0, 0, 1, 1])
    ax.set_axis_off()
    # plotting arrows
    # filtering every n-th value and every value smaller then x
    fx, fy, xs, ys = filter_values(fx, fy, abs_filter=filter[0], f_dist=filter[1],filter_method=filter_method, radius=filter_radius)
    if scale_ratio:  # optional custom scaling with the image axis lenght
        fx, fy = scale_for_quiver(fx, fy, dims=dims, scale_ratio=scale_ratio)
        quiver_parameters["scale"] = 1  # disabeling the auto scaling behavior of quiver
    plt.quiver(xs, ys, fx, fy, **quiver_parameters)  # plotting the arrows
    if plot_cbar:
        add_colorbar(vmin, vmax, cmap, ax=ax, cbar_style=cbar_style, cbar_width=cbar_width, cbar_height=cbar_height,
                     cbar_borderpad=cbar_borderpad, v=cbar_tick_label_size, cbar_str=cbar_str,
                     cbar_axes_fraction=cbar_axes_fraction, cbar_title_pad=cbar_title_pad)
    return fig, ax


def show_edgeline(mask, ax, color="#696969", alpha=0.5, n=6, plot_inner_line=False):
    '''
    colors a region close to the edge of mask.
    :param values: boolean mask
    :param ax: matplotlib axis object
    :param color: color of the edge coloring
    :param alpha: imshow alpha
    :param n: pixels distance from the ask edge
    :return:
    '''
    colony_area = mask != 0
    edge_area = np.logical_and(colony_area, ~binary_erosion(colony_area, iterations=n))
    edge_show = np.zeros(edge_area.shape)
    edge_show.fill(np.nan)
    edge_show[edge_area] = 1

    cmap_custom1 = matplotlib.colors.ListedColormap([color for i in range(3)])
    bounds = [0, 1, 2]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom1.N)
    ax.imshow(edge_show, cmap=cmap_custom1, norm=norm, alpha=alpha)

    if plot_inner_line:
        edge_line_inner = np.logical_and(~binary_erosion(colony_area, iterations=n),
                                         binary_erosion(colony_area, iterations=n - 1))
        edge_inner_show = np.zeros(edge_line_inner.shape)
        edge_inner_show.fill(np.nan)
        edge_inner_show[edge_line_inner] = 1

        cmap_custom2 = matplotlib.colors.ListedColormap(["black" for i in range(3)])
        bounds = [0, 1, 2]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap_custom2.N)
        ax.imshow(edge_inner_show, cmap=cmap_custom2, norm=norm)


def set_vmin_vmax(x, vmin, vmax):
    if not isinstance(vmin, (float, int)):
        vmin = np.nanmin(x)
    if not isinstance(vmax, (float, int)):
        vmax = np.nanmax(x)
    if isinstance(vmax, (float, int)) and not isinstance(vmin, (float, int)):
        vmin = vmax - 1 if vmin > vmax else None
    return vmin, vmax


def set_background(ax, fig, shape, background_color, cmap=None, cbar_style=""):
    ## setting background for the whoel plot or for the ax area
    background_color = matplotlib.cm.get_cmap(cmap)(0) if background_color == "cmap_0" else background_color
    if cbar_style == "clickpoints":
        fig.set_facecolor(background_color)  # setting background for the whole image, otherwise there are small white s
        # stripes at the edges
        ax.set_facecolor(background_color)
    else:
        custom_cmap = matplotlib.colors.ListedColormap([background_color, "white"])
        ax.imshow(np.zeros(shape), cmap=custom_cmap, vmin=0, vmax=1)


def show_map_clickpoints(values, figsize=(6.4, 4.8), cbar_str="", ax=None
                         , cmap="rainbow", vmin=None, vmax=None, background_color="cmap_0", cbar_width="2%",
                         cbar_height="50%",
                         cbar_borderpad=0.15, cbar_tick_label_size=15, cbar_axes_fraction=0.2
                         , cbar_style="clickpoints", plot_style="not-clickpoints", cbar_title_pad=1, plot_cbar=True,
                         show_mask=None, **kwargs):
    values = values.astype("float64")
    dims = values.shape  # save dims for use in scaling, otherwise porblems, because filtering will return flatten array
    values_show = np.zeros(dims)
    values_show.fill(np.nan)
    show_mask = show_mask.astype(bool) if isinstance(show_mask, np.ndarray) else values != 0
    values_show[show_mask] = values[show_mask]
    fig = plt.figure(figsize=figsize)
    if not isinstance(ax, matplotlib.axes.Axes):
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
    if plot_style == "clickpoints":
        ax.set_position([0, 0, 1, 1])

    ax.set_axis_off()
    set_background(ax, fig, dims, background_color, cmap=cmap, cbar_style=cbar_style)

    # other wise the program behaves unexpectedly if the automatically set vmin is bigger the manually set vmax
    vmin, vmax = set_vmin_vmax(values, vmin, vmax)
    im = ax.imshow(values_show, cmap=cmap, vmin=vmin, vmax=vmax)
    if plot_cbar:
        add_colorbar(vmin, vmax, cmap, ax=ax, cbar_style=cbar_style, cbar_width=cbar_width, cbar_height=cbar_height,
                     cbar_borderpad=cbar_borderpad, v=cbar_tick_label_size, cbar_str=cbar_str,
                     cbar_axes_fraction=cbar_axes_fraction, cbar_title_pad=cbar_title_pad)

    return fig, ax


def plot_map(ar1, cbar_str="", origin="upper", title="", mask=0, v_range=0, mask_overlay=0):
    '''
    imshow of an array with a color bar
    :return:
    '''

    ar = copy.deepcopy(ar1)
    if isinstance(mask, np.ndarray):
        mask = mask.astype(bool)
        ar[~mask] = np.nan

    fig, ax1 = plt.subplots(1, 1)
    axs = [ax1]
    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1] += 0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    if isinstance(v_range, tuple):
        min_v = v_range[0]
        max_v = v_range[1]
    else:
        min_v = np.nanmin(ar)
        max_v = np.nanmax(ar)

    axs[0].imshow(ar, origin=origin, cmap="jet", vmin=min_v, vmax=max_v)
    axs[0].set_title(title)

    if isinstance(mask_overlay, np.ndarray):  # showing overlay mask
        mask_overlay = mask_overlay.astype("float")  # cnat set to np.nan otherwise
        mask_overlay[mask_overlay == 0] = np.nan
        axs[0].imshow(mask_overlay, alpha=0.75)
    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    # norm=make_discrete_colorbar()
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str,
                                           orientation='vertical')
    format = ScalarFormatter()
    format.set_powerlimits((0, 0))  # setting scientific notation for color bar ## disable at some point
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str, format=format,
                                           orientation='vertical')

    return fig


def plot_fields(nodes, fields=[], titles=[], cbar_str=[], grid_lines=False, grid_lines_intervall=10, dims=None,
                origin="lower", mask=0, mask_overlay=0):
    if not dims:
        l = np.sqrt(len(nodes)).astype(int)
        fs = [np.zeros((l, l)) for i in
              range(len(fields))]  # note []*x produces x copies of the object with same reference
    else:
        fs = [np.zeros(dims) for i in
              range(len(fields))]  # note []*x produces x copies of the object with same reference

    for i in range(len(fs)):
        fs[i][nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = fields[i]  # inverted because of imshow

        # print(fs[i][nodes[:, 2].astype(int), nodes[:, 1].astype(int)])

    if isinstance(mask, np.ndarray):  # filling non displayed values with nan
        for i, f in enumerate(fs):
            fs[i][~mask] = np.nan

    if isinstance(mask_overlay, np.ndarray):  # preparing overlay mask
        mask_overlay = mask_overlay.astype("float")  # cant set to np.nan otherwise
        mask_overlay[mask_overlay == 0] = np.nan

    # setting plot subplots
    if len(fs) == 1:
        fig, ax1 = plt.subplots(1, 1)
        axs = [ax1]

    if len(fs) == 2:
        fig, axs = plt.subplots(1, 2)

    if len(fs) == 3:
        fig, axs = plt.subplots(1, 3)

    if len(fs) == 4:
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()

    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1] += 0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    plt.grid(grid_lines)

    min_v = np.nanmin(np.nanmin(np.dstack(fs), axis=2))
    max_v = np.nanmax(np.nanmax(np.dstack(fs), axis=2))

    # imshowing every array
    for a, f, t in zip(axs, fs, titles):
        a.imshow(f, origin=origin, cmap="jet", vmin=min_v, vmax=max_v)
        a.set_title(t)
        if grid_lines:
            show_grid(a)
            hide_ticks(a, interval=grid_lines_intervall)
        if isinstance(mask_overlay, np.ndarray):
            a.imshow(mask_overlay, alpha=0.75)
    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    # norm=make_discrete_colorbar()
    cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label=cbar_str,
                                           )
    return fig


def plot_all_sigmas(sigma_max, sigma_min, tau_max, sigma_avg, nodes):
    s_max = np.zeros((int(np.sqrt(len(nodes))), int(np.sqrt(len(nodes)))))
    s_min = copy.deepcopy((s_max))
    t_max = copy.deepcopy((s_max))
    sig_avg = copy.deepcopy((s_max))

    s_max[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_max  # should be the correct assignement
    s_min[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_min
    t_max[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = tau_max
    sig_avg[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = sigma_avg

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    plt.subplots_adjust(right=0.8)
    ax5 = fig.add_axes([0, 0, 0, 0])
    n_b = list(ax4.get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    ax5 = fig.add_axes(n_b)

    min_v = np.min(np.min(np.dstack([s_max, s_min, t_max, sig_avg]), axis=2))
    max_v = np.max(np.max(np.dstack([s_max, s_min, t_max, sig_avg]), axis=2))

    ax1.imshow(s_max, origin='lower', cmap="jet", vmin=min_v, vmax=max_v)  ## without any interpolation and such
    ax2.imshow(s_min, origin='lower', cmap="jet", vmin=min_v, vmax=max_v)
    ax3.imshow(t_max, origin='lower', cmap="jet", vmin=min_v, vmax=max_v)  #
    ax4.imshow(sig_avg, origin='lower', cmap="jet", vmin=min_v, vmax=max_v)
    ax1.set_title("sigma_max")
    ax2.set_title("sigma_min")
    ax3.set_title("tau_max")
    ax4.set_title("sigma_avg")
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    cb1 = matplotlib.colorbar.ColorbarBase(ax5, cmap=matplotlib.cm.get_cmap("jet"),
                                           norm=norm, label="stress",
                                           orientation='vertical')


def construct_elements(sqr1, nodes, elements):
    l = np.zeros(len(elements))
    for i, p in enumerate(zip(sqr1[0], sqr1[1])):
        l[i] = np.where((nodes[:, 1] == p[0]) * (nodes[:, 2] == p[1]))[0]  # finding corresponding point
    return l


def vizualize_forces_on_areas(tx_resize, ty_resize, areas, mask):
    fig = plt.figure()
    mask_show = np.zeros(mask.shape) + np.nan
    plt.imshow(np.sqrt((tx_resize / 1000) ** 2 + (ty_resize / 1000) ** 2), cmap="rainbow")
    cbar = plt.colorbar()
    cbar.set_label("traktion forces in kPa")
    pixx = np.arange(mask.shape[0])
    pixy = np.arange(mask.shape[1])
    xv, yv = np.meshgrid(pixy, pixx)
    select_x = ((xv - 1) % 50) == 0
    select_y = ((yv - 1) % 50) == 0
    y = select_x[0, :].sum()
    x = select_y[:, 0].sum()  ## this is correct
    select = select_x * select_y
    tx_show = tx_resize[select]
    ty_show = ty_resize[select]
    x1 = np.where(select)[1].reshape((x, y))
    y1 = np.where(select)[0].reshape((x, y))
    ## maek beter xs...
    scale_ratio = 0.2
    scale = scale_ratio * np.max(np.shape(tx_resize)) / np.max(
        np.sqrt((tx_show / 1000) ** 2 + (ty_show / 1000) ** 2))  # automatic sacleing in dependace of the image size #
    plt.quiver(x1, y1, (tx_show / 1000) * scale, (ty_show / 1000) * scale, scale=1, scale_units='xy', angles="xy",
               width=0.002)
    for i, (key, value) in enumerate(areas.items()):
        mask_show[value[0]] = i
    plt.imshow(mask_show, alpha=0.8, cmap="magma")
    scale = 0.5 * 10 ** -4
    for key, value in areas.items():
        plt.arrow(value[4][0], value[4][1], value[5][0] * scale, value[5][1] * scale, head_width=20, color="red")


def find_areas(start_line, lines, i, com_all, invert_direction=False):
    circ_line = []
    line_ids = []
    id = i
    line = np.array(start_line)

    # fixing direction:
    v1 = line[1] - line[0]
    v2 = com_all - line[1]  # vector from line end to center of mass
    cross = (np.cross(v2, v1) > 0) * 2 - 1  # gets angle direction
    angle = np.arccos(np.dot(v2, v1) / (np.linalg.norm(v2) * np.linalg.norm(v1))) * cross
    direction_factor = (angle < 0) * 2 - 1  # minus one if positive angle towards center, else negative value
    if invert_direction: direction_factor *= -1  # used
    # plt.figure()
    # plt.imshow(mask)
    # plt.arrow(line[0][0], line[0][1], v1[0], v1[1], head_width=20)
    # plt.arrow(line[1][0], line[1][1],  v2[0],  v2[1], head_width=20)

    check = False

    while not check:
        circ_line.append(line)
        line_ids.append(id)

        # logical and operation to find where both coordinates of a point are close
        child_lines1 = np.where(
            np.isclose(line[1], lines)[:, :, 0] * np.isclose(line[1], lines)[:, :, 1])  # lines at the end
        child_lines2 = np.unique(child_lines1[0])
        child_lines3 = lines[child_lines2, :, :]
        # reoreintating child lines to get uniform direction
        child_lines = np.array(
            [[l[int(~np.isclose(l[0], line[1]).all())], l[int(np.isclose(l[0], line[1]).all())]] for l in child_lines3])
        # finding left most line
        child_vec = child_lines[:, 1, :] - child_lines[:, 0, :]
        line_vec = line[1] - line[0]
        # filtering
        # angles = np.arcsin(np.abs(child_vec[:,0] * line_vec[1] - child_vec[:,1] * line_vec[0]))
        cross = (np.cross(child_vec, line_vec) > 0) * 2 - 1  # gets angle direction
        angles = np.arccos(
            np.dot(child_vec, line_vec) / (np.linalg.norm(child_vec, axis=1) * np.linalg.norm(line_vec))) * cross
        #####
        angles[np.isclose(-angles, np.pi)] = np.nan
        angles[np.isclose(angles, np.pi)] = np.nan
        print(angles)
        id = child_lines2[np.nanargmin(angles * direction_factor)]
        line = child_lines[np.nanargmin(angles * direction_factor)]  ## maybe exclude zero??
        ##reoreintate line if necessary, new point goes first
        check = np.array([np.array([np.isclose(l[0], line[0]), np.isclose(l[1], line[1])]).all() for l in
                          circ_line]).any()  # chekcing if finisched by comparing sum of points ## doesnt work............
        # plt.arrow(line[0][0], line[0][1], line[1][0] - line[0][0], line[1][1] - line[0][1], head_width=20)
    return circ_line, line_ids
    # plt.text(line[0][0], line[0][1],str(np.round(angles[np.nanargmin(angles)],4)))


def plot_arrows(nodes, x, y, cbar_str=[], scale_ratio=0.05, grid_lines=False, dims=None, origin="lower", title="",
                mask=0, filter=0, overlay_mask=0):
    '''

    :param nodes:
    :param x:
    :param y:
    :param cbar_str:
    :param scale_ratio:
    :param grid_lines:
    :param dims:
    :param origin:
    :param title:
    :param mask:
    :param filter: list of pararmeters to filter dispplayed values:[minimal abs value, spacing between arrows]
    :return:
    '''

    if not dims:
        l = np.sqrt(len(nodes)).astype(int)
        fi_x = np.zeros((l, l))
        fi_y = np.zeros((l, l))
        dims = (l, l)

    else:
        fi_x = np.zeros(dims)  # note []*x produces x copies of the object with same reference
        fi_y = np.zeros(dims)

    fi_x[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = x  # inverted because of imshow
    fi_y[nodes[:, 2].astype(int), nodes[:, 1].astype(int)] = y
    if isinstance(mask, np.ndarray):
        fi_x[~mask] = np.nan
        fi_y[~mask] = np.nan

    if isinstance(filter, list):
        fi_x_f, fi_y_f, xs, ys = filter_values(fi_x, fi_y, abs_filter=filter[0], f_dist=filter[1])
        fi_x_f, fi_y_f = scale_for_quiver(fi_x_f, fi_y_f, fi_y_f.shape, scale_ratio=scale_ratio)
        min_v = np.nanmin(np.sqrt(fi_x ** 2 + fi_y ** 2))
        max_v = np.nanmax(np.sqrt(fi_x ** 2 + fi_y ** 2))
    else:
        fi_x, fi_y = scale_for_quiver(fi_x, fi_y, fi_x.shape, scale_ratio=scale_ratio)
        min_v = np.nanmin(np.sqrt(fi_x ** 2 + fi_y ** 2))
        max_v = np.nanmax(np.sqrt(fi_x ** 2 + fi_y ** 2))

    fig, ax1 = plt.subplots(1, 1)
    axs = [ax1]
    plt.subplots_adjust(right=0.8)
    ax_cbar = fig.add_axes([0, 0, 0, 0])
    n_b = list(axs[-1].get_position().bounds)
    n_b[0] += n_b[2] + 0.02
    n_b[2] = 0.03
    n_b[1] += 0.2
    n_b[3] -= 0.4

    ax_cbar = fig.add_axes(n_b)
    plt.grid(grid_lines)

    axs[0].imshow(np.sqrt(fi_x ** 2 + fi_y ** 2), origin=origin, cmap="jet")
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)

    if isinstance(filter, list):
        axs[0].quiver(xs, ys, fi_x_f, fi_y_f, scale_units="xy", scale=1, angles="xy")
    else:
        axs[0].quiver(fi_x, fi_y, scale_units="xy", scale=1, angles="xy")
    axs[0].set_title(title)

    # plotting additional mask if desired
    if isinstance(overlay_mask, np.ndarray):
        ov_mask_show = np.zeros(overlay_mask.shape) + np.nan
        ov_mask_show[overlay_mask] = 1
        axs[0].imshow(ov_mask_show, alpha=0.5)

    # plotting one color bar
    norm = matplotlib.colors.Normalize(vmin=min_v, vmax=max_v)
    # norm=make_discrete_colorbar()
    if len(np.unique(np.sqrt(fi_x ** 2 + fi_y ** 2))) > 1:  # no color bar if only one absolute value
        cb1 = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=matplotlib.cm.get_cmap("jet"),
                                               norm=norm, label=cbar_str,
                                               orientation='vertical')
    else:
        ax_cbar.remove()

    return fig



def find_maxima(ar1,ar2,radius=5,shape="circle"):
    # generating circle

    ys,xs = np.indices((radius*2 + 1,radius*2+1))
    xs = (xs - radius).astype(float)
    ys = (ys - radius).astype(float)
    if shape=="circle":
        out = np.sqrt(xs ** 2 + ys ** 2) <= radius
        xs[~out] = np.nan
        ys[~out] = np.nan
    abs = np.sqrt(ar1**2+ar2**2)
    lmax = np.unravel_index(np.nanargmax(abs),shape=abs.shape)
    maxis  = [lmax]
    while True:
        x_exclude = (lmax[1] + xs).flatten()
        y_exclude = (lmax[0] + ys).flatten()
        outside_image = (x_exclude>=abs.shape[1]) | (x_exclude<0) |  (y_exclude>=abs.shape[0]) | (y_exclude<0) | (np.isnan(x_exclude)) | (np.isnan(y_exclude))
        x_exclude = x_exclude[~outside_image]
        y_exclude = y_exclude[~outside_image]
        abs[y_exclude.astype(int),x_exclude.astype(int)] = np.nan
        try:
            lmax = np.unravel_index(np.nanargmax(abs), shape=abs.shape)
        except ValueError:
            break
        maxis.append(lmax)

    maxis_y = [i[0] for i in maxis]
    maxis_x = [i[1] for i in maxis]
    return maxis_y, maxis_x




def filter_values(ar1, ar2, abs_filter=0, f_dist=3, filter_method="regular", radius=5):
    '''
    function to filter out values from an array for better display
    :param ar1:
    :param ar2:
    :param ar:
    :param f_dist: distance betweeen filtered values
    :return:
    '''


    if filter_method == "regular":
        pixx = np.arange(np.shape(ar1)[0])
        pixy = np.arange(np.shape(ar1)[1])
        xv, yv = np.meshgrid(pixy, pixx)

        def_abs = np.sqrt((ar1 ** 2 + ar2 ** 2))
        select_x = ((xv - 1) % f_dist) == 0
        select_y = ((yv - 1) % f_dist) == 0
        select_size = def_abs > abs_filter
        select = select_x * select_y * select_size
        s1 = ar1[select]
        s2 = ar2[select]
        x_ind = xv[select]
        y_ind = yv[select]
    if filter_method == "local_maxima":
        y_ind,x_ind = find_maxima(ar1, ar2, radius=radius,shape="circle")
        s1 = ar1[y_ind, x_ind]
        s2 = ar2[y_ind, x_ind]
    if filter_method == "local_maxima_square":
        y_ind,x_ind = find_maxima(ar1, ar2, radius=radius,shape="square")
        s1 = ar1[y_ind, x_ind]
        s2 = ar2[y_ind, x_ind]
    return s1, s2, x_ind, y_ind


def check_closet_neigbours(points1, points2, assign, mask1=None, mask2=None):
    plt.figure()
    if isinstance(mask1, np.ndarray) and isinstance(mask2, np.ndarray):
        plt.imshow(mask1.astype(int) + mask2.astype(int))

    for p1_ind, p1 in enumerate(points1):
        plt.text(p1[1], p1[0], str(p1_ind))
    for p2_ind, p1_ind in enumerate(assign):
        plt.text(points2[p2_ind][1], points2[p2_ind][0], str(p1_ind))


def get_line(start, end):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end
    start end as tupels of points
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points


def scale_for_quiver(ar1, ar2, dims, scale_ratio=0.2, return_scale=False):
    scale = scale_ratio * np.max(dims) / np.nanmax(np.sqrt((ar1) ** 2 + (ar2) ** 2))
    if return_scale:
        return scale
    return ar1 * scale, ar2 * scale
