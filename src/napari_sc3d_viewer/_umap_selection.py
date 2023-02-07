"""
This file is subject to the terms and conditions defined in
file 'LICENCE', which is part of this source code package.
Author: Leo Guignard (leo.guignard...@AT@...univ-amu.fr)
"""

from qtpy.QtWidgets import QWidget, QVBoxLayout
from ._utils import error_points_selection, safe_toarray
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
)
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path as PathMPL
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# This class is directly copied from there:
# https://matplotlib.org/stable/gallery/widgets/lasso_selector_demo_sgskip.html
# I got the info from there:
# https://github.com/BiAPoL/napari-clusters-plotter
# <3 <3 <3
class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError("Collection must have a facecolor")
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect, button=1)
        self.ind = []

    def onselect(self, verts):
        path = PathMPL(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.ind_mask = path.contains_points(self.xys)
        self.canvas.draw_idle()
        self.selected_coordinates = self.xys[self.ind].data

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


class UmapSelection:
    """
    Display a pre-computed umap and allow for selection within the umap that
    will be reflected on the parent napari viewer
    """

    def get_stats(self, indices):
        """
        Given a stat function (standard deviation for example) and cell indices,
        compute this for each gene expression over the set of selected cells
        """
        # If only variable are asked, do not take the raw data
        if self.variable_genes.value:
            data = self.embryo.anndata
            gene_stat = self.stat_func(data[indices].X, axis=0)
        else:
            data = self.embryo.anndata.raw
            gene_stat = self.stat_func(safe_toarray(data[indices].X), axis=0)

        top_genes = np.argsort(gene_stat)[-3:]
        top_gene_names = data.var_names[top_genes]
        return top_genes, top_gene_names

    def show_cells(self, event):
        """
        Handle the showing of cells according to a lasso selection
        """
        points = self.points
        # when an event is captured
        if event:
            points.shown = np.zeros_like(points.shown, dtype=bool)
            # If no cells where within the lasso selection, reset the view
            if sum(len(s.ind) for s in self.selectors) == 0:
                points.shown = points.features["current_view"]
                [pt.set_alpha(1) for pt in self.pts]
                self.ax_G.set_title(f"Gene: {self.gene.value}")
                for i, (ax_hist, (gene_hist, maxi)) in enumerate(
                    self.ax_hists
                ):
                    vals = self.embryo.anndata.raw[:, gene_hist]
                    ax_hist.clear()
                    vals = safe_toarray(vals.X)[
                        points.features["current_view"], 0
                    ]
                    hist = ax_hist.hist(vals, bins=50)
                    ax_hist.set_yticks([])
                    if self.tissues.value:
                        ax_hist.set_yticks([])
                        ax_hist.set_xlabel("Gene expression")
                        if i == 0:
                            ax_hist.set_ylabel("#cells")
                    else:
                        ax_hist.set_ylabel("#cells")
                        if i == 2:
                            ax_hist.set_xlabel("Gene expression")
                    ax_hist.set_title(f"{gene_hist} distribution")
                    ax_hist.set_xlim(0, max(0.01, maxi))
            # If a set of cells where within the lasso selection, display them
            else:
                # If both gene and tissues are displayed, there are 2 selectors
                for s in self.selectors:
                    if len(s.ind) != 0:
                        # Get, count and show the selected cells
                        indices = self.corres_to_mask[self.sorted_vals][s.ind]
                        nb_init = np.sum(points.features["current_view"])
                        points.shown[indices] = points.features[
                            "current_view"
                        ][indices]
                        nb = np.sum(points.shown[indices])
                        self.ax_G.set_title(
                            f"Gene: {self.gene.value} ({nb} cells "
                            f"({100*nb/nb_init:.1f}% of the initial))"
                        )
                        (top_genes, top_gene_names) = self.get_stats(indices)
                        # Display the distribution of the highest `stat` distribution
                        # for the selected cells
                        for j, (ax_hist, _) in enumerate(self.ax_hists):
                            gene_name = top_gene_names[j]
                            gene_id = top_genes[j]
                            vals = self.embryo.anndata.raw[indices, gene_name]
                            ax_hist.clear()
                            hist = ax_hist.hist(
                                safe_toarray(vals.X)[:, 0], bins=50
                            )
                            ax_hist.set_yticks([])
                            if self.tissues.value:
                                ax_hist.set_xlabel("Gene expression")
                                if j == 0:
                                    ax_hist.set_ylabel("#cells")
                            else:
                                ax_hist.set_ylabel("#cells")
                                if j == 2:
                                    ax_hist.set_xlabel("Gene expression")
                            ax_hist.set_title(f"{gene_name} distribution")
                            ax_hist.set_xlim(
                                0, max(0.01, self.maximums[gene_id])
                            )
                        alpha = np.zeros_like(self.sorted_vals) + 0.1
                        alpha[s.ind] = 1
                        [pt.set_alpha(alpha) for pt in self.pts]
                        s.ind = []
        self.points.refresh()

    def build_figure(self):
        """
        Build the canvas and the figures that will held the umaps and histograms
        """
        mask = self.mask
        # The main figure
        static_canvas = FigureCanvas()
        fig = static_canvas.figure
        static_canvas.toolbar = NavigationToolbar(
            static_canvas, static_canvas.parent()
        )
        fig.set_figwidth(10)
        fig.set_figheight(8)
        self.ax_hists = []

        # Set the grid and axes according to whether the tissues
        # have to be shown or not
        if self.tissues.value:
            gs0 = GridSpec(3, 1, figure=fig)
            gs00 = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[:2, 0])
            self.ax_G = fig.add_subplot(gs00[0, 0])
            ax_T = fig.add_subplot(gs00[0, 1])
            gs01 = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[2, 0])
            for i in range(3):
                ax = fig.add_subplot(gs01[0, i])
                ax.set_xlabel("Gene expression")
                if i == 0:
                    ax.set_ylabel("#cells")
                self.ax_hists.append([ax, self.gene.value])
        else:
            gs = GridSpec(3, 3, figure=fig)
            self.ax_G = fig.add_subplot(gs[:, :2])
            for i in range(3):
                ax = fig.add_subplot(gs[i, -1])
                if i == 2:
                    ax.set_xlabel("Gene expression")
                ax.set_ylabel("#cells")
                self.ax_hists.append([ax, self.gene.value])

        # Scatter of the umap for gene expression
        val_g = safe_toarray(self.embryo.anndata.raw[:, self.gene.value].X)[
            mask, 0
        ]
        self.maximums = safe_toarray(self.embryo.anndata.raw.X)[mask].max(
            axis=0
        )
        colors = val_g
        pos_cells = self.embryo.anndata.obsm[self.embryo.umap_id][mask, :2]
        self.sorted_vals = np.argsort(val_g)
        min_c, max_c = colors.min(), colors.max()
        colors = (colors - min_c) / (max_c - min_c)
        pts_G = self.ax_G.scatter(
            *pos_cells[self.sorted_vals, :].T,
            marker=".",
            c=colors[self.sorted_vals],
        )
        self.pts = [pts_G]
        (top_genes, top_gene_names) = self.get_stats(
            self.points.features["current_view"]
        )
        # Histograms of the most variable genes
        for j, (ax_hist, _) in enumerate(self.ax_hists):
            gene_hist = top_gene_names[j]
            gene_id = top_genes[j]
            vals = self.embryo.anndata.raw[:, gene_hist]
            ax_hist.set_yticks([])
            ax_hist.hist(safe_toarray(vals.X)[mask, 0], bins=50)
            self.ax_hists[j][1] = (gene_hist, self.maximums[gene_id])
            ax_hist.set_title(f"{gene_hist} distribution")
            ax_hist.set_xlim(0, max(0.001, self.maximums[gene_id]))
        pts_G.set_edgecolor("none")
        self.ax_G.set_xticks([])
        self.ax_G.set_yticks([])
        self.ax_G.set_xlabel("umap 1")
        self.ax_G.set_ylabel("umap 2")
        self.ax_G.set_title(f"Gene: {self.gene.value}")
        self.ax_G.set_aspect("equal")

        # Create and add the lasso selector
        self.selectors = []
        self.selectors.append(SelectFromCollection(self.ax_G, pts_G))

        # Scatter plot of the umap according to tissue type if asked for
        if self.tissues.value:
            colors_T = [
                self.color_map_tissues.get(self.embryo.tissue[c], [0, 0, 0])
                for c in self.corres_to_mask
            ]
            colors_T = np.array(colors_T)
            pts_T = ax_T.scatter(
                *pos_cells[self.sorted_vals, :].T,
                marker=".",
                color=colors_T[self.sorted_vals],
            )
            self.pts.append(pts_T)
            pts_T.set_edgecolor("none")
            ax_T.set_xticks([])
            ax_T.set_yticks([])
            ax_T.set_ylabel("umap 2")
            ax_T.set_title(f"Tissues")
            tissues_found = set(
                [self.embryo.tissue[c] for c in self.corres_to_mask]
            )
            for t in tissues_found:
                ax_T.plot(
                    [],
                    [],
                    "o",
                    color=self.color_map_tissues.get(t, [0, 0, 0]),
                    label=self.embryo.corres_tissue.get(t, f"{t}"),
                )
            self.selectors.append(SelectFromCollection(ax_T, pts_T))
            ax_T.set_aspect("equal")
            ax_T.legend(fontsize="xx-small", frameon=False, shadow=False)

        fig.canvas.mpl_connect("button_release_event", self.show_cells)
        # That tight_layout does not work, I am not sure why ...
        fig.tight_layout()

        return static_canvas

    def run(self):
        """
        Build the umap figure and add it to the correct tab in the viewer
        """
        # Retrieve the parameters
        if self.stats.value == "Standard Deviation":
            self.stat_func = np.std
        elif self.stats == "Mean":
            self.stat_func = np.mean
        elif self.stats == "Median":
            self.stat_func = np.median
        else:
            self.stat_func = np.max

        # Make sure that points and parameters are actually correct
        if (
            self.points is None
            or self.points.as_layer_data_tuple()[-1] != "points"
        ):
            error_points_selection()
            return

        if not hasattr(self.points.features, "current_view"):
            self.points.features["current_view"] = self.points.shown.copy()

        if (
            self.points is None
            or not self.gene.value in self.embryo.anndata.raw.var_names
        ):
            return f"'{self.gene.value}' not found"

        if 0 < len(self.points.selected_data):
            mask = np.zeros_like(self.points.shown)
            mask[list(self.points.selected_data)] = True
            self.points.shown[~mask] = False
            self.points.refresh()
        else:
            mask = self.points.shown
        self.mask = mask
        self.corres_to_mask = np.where(mask)[0]

        # Create the figure and the widget containing it
        static_canvas = self.build_figure()
        fig_can = self.viewer.window.add_dock_widget(
            static_canvas, name="umap"
        )
        V_box = QWidget()
        V_box.setLayout(QVBoxLayout())
        V_box.layout().addWidget(fig_can)
        V_box.layout().addWidget(static_canvas.toolbar)
        self.tab2.removeTab(self.tab2.nb_tabs + 1)
        self.tab2.addTab(V_box, "umap graph")

    def __init__(
        self,
        viewer,
        embryo,
        gene,
        tissues,
        stats,
        variable_genes,
        color_map_tissues,
        tab2,
    ):
        """
        Creation of the umap widget

        Args:
            viewer (napari.Viewer): napari viewer containing the `sc3D` points
            embryo (sc3D.Embryo): embryo to display
            gene (str): gene to display on the umap
            tissues (bool): whether or not to display the umap with the tissues
            stats (str ['Standard Deviation' | 'Mean' | 'Median']): Stat to compute
            variable_genes (bool): whether to display only variable genes or not
            color_map_tissues (dict): dictionnary mapping tissues id to their colors,
                necessary if `tissues` is `True`
            tab2 (qtpy.QtWidgets.QTabWidget): tab where to insert the produced
                scatter plot. Need `tab2.nb_tabs` parameter
        """
        super().__init__()
        self.viewer = viewer
        self.points = self.viewer.layers.selection.active
        self.embryo = embryo
        self.stats = stats
        self.tissues = tissues
        self.gene = gene
        self.variable_genes = variable_genes
        self.color_map_tissues = color_map_tissues
        self.tab2 = tab2
        mpl.rcParams["font.size"] = 6
