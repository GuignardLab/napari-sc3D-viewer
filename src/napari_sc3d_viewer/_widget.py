"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/plugins/guides.html?#widgets

Replace code below according to your needs.
"""
import json
from sc3D import Embryo
from qtpy.QtWidgets import (QWidget,
                            QVBoxLayout,
                            QHBoxLayout,
                            QTabWidget,
                            QFileDialog,
                            QMessageBox,
                            QPushButton)
from magicgui import magicgui
from magicgui.widgets import FileEdit, LineEdit
from pathlib import Path
from napari import Viewer
from napari.utils.colormaps import ALL_COLORMAPS
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import cm, colors
from collections.abc import Iterable
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import LassoSelector, TextBox
from matplotlib.path import Path as PathMPL
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

def error_points_selection():
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText('Point cloud selection error')
    msg.setInformativeText(('Please select an adequate point cloud\n'
                            'You can select point clouds on the left hand side of the viewer'))
    msg.setWindowTitle('Point cloud selection error')
    msg.exec_()

def safe_toarray(tab):
    if not isinstance(tab, np.ndarray):
        return tab.toarray()
    else:
        return tab

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
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect, button=1)
        self.ind = []

    def onselect(self, verts):
        path = PathMPL(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.ind_mask = path.contains_points(self.xys)
        # self.fc[:, -1] = self.alpha_other
        # self.fc[self.ind, -1] = 1
        # self.collection.set_alpha(self.fc[:, -1])
        self.canvas.draw_idle()
        self.selected_coordinates = self.xys[self.ind].data

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

def display_embryo(viewer, embryo):
    color_map_tissues = {
        5: [0.7411764705882353, 0.803921568627451, 1.0],
        6: [0.19607843137254902, 0.35294117647058826, 0.6078431372549019],
        7: [0.996078431372549, 0.6862745098039216, 0.08627450980392157],
        9: [0.7686274509803922, 0.27058823529411763, 0.10980392156862745],
        10: [0.10980392156862745, 1.0, 0.807843137254902],
        12: [0.7529411764705882, 0.4588235294117647, 0.6509803921568628],
        13: [0.9647058823529412, 0.13333333333333333, 0.1803921568627451],
        14: [0.7411764705882353, 0.43529411764705883, 0.6705882352941176],
        15: [0.9686274509803922, 0.8823529411764706, 0.6274509803921569],
        16: [1.0, 0.9803921568627451, 0.9803921568627451],
        18: [0.47058823529411764, 0.16470588235294117, 0.7137254901960784],
        20: [0.5019607843137255, 0.5019607843137255, 0.5019607843137255],
        21: [0.9803921568627451, 0.0, 0.5294117647058824],
        22: [0.5098039215686274, 0.1803921568627451, 0.10980392156862745],
        23: [0.5215686274509804, 0.4, 0.050980392156862744],
        24: [0.803921568627451, 0.1607843137254902, 0.5647058823529412],
        27: [0.6588235294117647, 0.6588235294117647, 0.6588235294117647],
        29: [0.0, 0.0, 0.5450980392156862],
        30: [0.5450980392156862, 0.2784313725490196, 0.36470588235294116],
        31: [1.0, 0.7568627450980392, 0.1450980392156863],
        32: [0.8705882352941177, 0.6274509803921569, 0.9921568627450981],
        33: [0.19607843137254902, 0.5137254901960784, 0.996078431372549],
        34: [0.9725490196078431, 0.6313725490196078, 0.6235294117647059],
        35: [0.7098039215686275, 0.9372549019607843, 0.7098039215686275],
        36: [0.1803921568627451, 0.8509803921568627, 1.0],
        39: [0.10980392156862745, 0.5137254901960784, 0.33725490196078434],
        40: [1.0, 0.6470588235294118, 0.30980392156862746],
        41: [0.8470588235294118, 0.7490196078431373, 0.8470588235294118]
    }
    tissues_to_plot = [18, 21, 30, 31, 34]
    tissues_to_plot = [t for t in tissues_to_plot if t in embryo.all_tissues]
    if len(tissues_to_plot)<1:
        tissues_to_plot = list(embryo.all_tissues)
    cells = sorted(embryo.all_cells)
    positions = [embryo.pos_3D[c] for c in cells]
    shown = [embryo.tissue[c] in tissues_to_plot for c in cells]
    if not any(shown):
        shown = [True]*len(cells)
    properties = {'cells': cells}

    properties['gene'] = [0 for _ in cells]
    if 0<len(embryo.all_tissues.difference(color_map_tissues)):
        nb_tissues = len(embryo.all_tissues)
        color_map_tissues = {v:cm.tab20(i/nb_tissues)
                                for i, v in enumerate(embryo.all_tissues)}

    colors_rgb = [color_map_tissues.get(embryo.tissue[c], [0, 0, 0]) for c in cells]

    viewer.dims.ndisplay=3
    points = viewer.add_points(positions, face_color=colors_rgb,
                               properties=properties,
                               metadata={'gene': None, '2genes': None}, shown=shown)

    all_tissues = [embryo.corres_tissue.get(t, f'{t}')
                   for t in embryo.all_tissues]
    @magicgui(call_button='Select tissues',
              tissues={'widget_type': 'Select',
                         'choices': all_tissues,
                         'value': [embryo.corres_tissue.get(t, f'{t}')
                                   for t in tissues_to_plot]})
    def select_tissues(viewer: Viewer, tissues):
        tissue_to_num = {v:k for k, v in embryo.corres_tissue.items()}
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        tissues_to_plot = []
        for t in tissues:
            if t in tissue_to_num:
                tissues_to_plot.append(tissue_to_num[t])
            else:
                tissues_to_plot.append(int(t))
        shown = [embryo.tissue[c] in tissues_to_plot for c in points.properties['cells']]
        points.shown = shown
        points.features['current_view'] = shown
        if points.metadata['gene'] is None and points.metadata['2genes'] is None:
            show_tissues(viewer)
        elif points.metadata['2genes'] is None:
            show_gene(viewer, points.metadata['gene'])
        else:
            show_two_genes(viewer, *points.metadata['2genes'])

    @magicgui(auto_call=True,
              cmap={'label': 'colormap',
                    'choices': ALL_COLORMAPS.keys()})
    def apply_cmap(viewer: Viewer, cmap: str):
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        if len(points.properties) == 0:
            return
        if points.face_color_mode.lower() != 'colormap':
            points.face_color = 'gene'
            points.face_color_mode = 'colormap'
        points.face_colormap = cmap
        points.refresh()

    @magicgui(call_button='Show gene/metric',
              metric={'widget_type': 'ComboBox',
                      'label': 'What do you want to display:',
                      'choices': (['Gene']+
                                  [c for c in list(embryo.anndata.obs.columns)
                                   if embryo.anndata.obs[c].dtype in [float, int]]),
                      'value': 'Gene'},
              gene={'label': 'Which gene (if necessary)', 'value': 'T'},
              result_widget=True)
    def show_gene(viewer: Viewer, metric: str, gene: str):
        points = viewer.layers.selection.active
        cell_list = list(embryo.all_cells)
        is_metric = metric in embryo.anndata.obs.columns
        if is_metric:
            gene = metric
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return f'Wrong point selection'
        if (not gene in embryo.anndata.obs.columns and
            not gene in embryo.anndata.raw.var_names):
            return f"'{gene}' not found"
        if gene != points.metadata['gene']:
            if 'current_view' in points.features:
                mask = points.features['current_view']
            else:
                mask = points.shown
            if is_metric:
                colors = embryo.anndata.obs[metric].to_numpy()
                try:
                    mask &= ~np.isnan(colors)
                except Exception as e:
                    print(colors.dtype)
                    return('Failed')
                points.shown = mask
            else:
                colors = safe_toarray(embryo.anndata.raw[:, gene].X)[:, 0]
            min_c, max_c = colors[mask].min(), colors[mask].max()
            colors = (colors-min_c)/(max_c-min_c)
            colors[~mask] = 0
            points.features['gene'] = colors
            points.metadata['gene_min_max'] = min_c, max_c
            points.metadata['gene'] = gene
        points.metadata['2genes'] = None
        points.edge_color = 'black'
        points.face_color = 'gene'
        points.face_color_mode = 'colormap'
        points.face_contrast_limits = (0, 1)
        points.refresh()
        return f"{points.metadata['gene']} displayed"

    @magicgui(call_button='Show gene',
              gene1={'label': 'First gene (main)', 'value': 'T'},
              gene2={'label': 'Second gene', 'value': 'Sox2'},
              low_th={'label': 'Low threshold',
                      'value': 2, 'min':0, 'max':100},
              high_th={'label': 'High threshold',
                       'value': 98, 'min':0, 'max':100},
              main_bi_color={'label': 'Main color', 'widget_type': 'ComboBox',
                             'choices': ['Red', 'Green', 'Blue'],
                             'value': 'Red'},
              result_widget=True)

    def show_two_genes(viewer: Viewer, gene1: str, gene2: str,
                       low_th: float, high_th: float,
                       main_bi_color: str):
        points = viewer.layers.selection.active
        cell_list = list(embryo.all_cells)
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return 'A point layout should be selected'
        if not gene1 in embryo.anndata.raw.var_names:
            return f"'{gene1}' not found"
        if not gene2 in embryo.anndata.raw.var_names:
            return f"'{gene2}' not found"
        if (not points.metadata['2genes'] or
            (gene1, gene2, low_th,
             high_th, main_bi_color) != points.metadata['2genes']):
            if 'current_view' in points.features:
                mask = points.features['current_view']
            else:
                mask = points.shown
            colors1 = safe_toarray(embryo.anndata.raw[:, gene1].X)[:, 0]
            colors2 = safe_toarray(embryo.anndata.raw[:, gene2].X)[:, 0]
            C = np.array([colors1, colors2])
            min_g1 = np.percentile(C[0][mask], low_th)
            min_g2 = np.percentile(C[1][mask], low_th)
            max_g1 = np.percentile(C[0][mask], high_th)
            max_g2 = np.percentile(C[1][mask], high_th)
            norm = lambda C: (C-[[min_g1], [min_g2]]) / [[max_g1-min_g1], [max_g2-min_g2]]
            V = norm(C)
            V[V<0] = 0
            V[1<V] = 1
            final_C = np.zeros((len(colors1), 3))
            on_channel = (np.array(['Red', 'Green', 'Blue'])!=main_bi_color).astype(int)
            final_C[:,0] = V[on_channel[0]]
            final_C[:,1] = V[on_channel[1]]
            final_C[:,2] = V[on_channel[2]]
            points.face_color = final_C
            points.face_color_mode = 'direct'
            points.features['2genes'] = colors1
            points.metadata['2genes'] = (gene1, gene2, low_th, high_th, main_bi_color)
            points.metadata['2genes_params'] = (max_g1, max_g2, norm, on_channel)
            points.metadata['gene'] = None
            points.edge_color = 'black'
        return ', '.join(points.metadata['2genes'][:2])

    @magicgui(call_button='Color according to cell types')
    def show_tissues(viewer: Viewer):
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        if (points.metadata['gene'] is not None or
            points.metadata['2genes'] is not None):
            points.face_color = [color_map_tissues[embryo.tissue[c]]
                                 for c in points.properties['cells']]
            points.face_color_mode = 'direct'
            points.metadata['gene'] = None
            points.metadata['2genes'] = None
        points.refresh()

    @magicgui(call_button='Display color legend')
    def disp_legend(viewer: Viewer):
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        with plt.style.context('dark_background'):
            fig, ax = plt.subplots()
            if (points.metadata['gene'] is None and
                points.metadata['2genes'] is None):
                tissues = set([embryo.tissue[c] for c in
                               points.properties['cells'][points.shown]])
                for t in tissues:
                    ax.plot([], 'o', c=color_map_tissues[t], label=embryo.corres_tissue.get(t, f'{t}'))
                ax.legend()
                ax.set_axis_off()
            elif points.metadata['2genes'] is None:
                if points.face_contrast_limits is None:
                    m, M = 0, 1
                else:
                    m, M = points.face_contrast_limits
                if points.face_colormap.name in plt.colormaps():
                    plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(m, M),
                                                   cmap=points.face_colormap.name))
                    min_, max_ = points.metadata['gene_min_max']
                    min_ = (max_-min_)*m+min_
                    max_ = (max_-min_)*M+min_
                    plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(min_, max_),
                                                   cmap=points.face_colormap.name))
                else:
                    plt.text(0, 0, ( 'Could not find the colormap '
                                    f'`{points.face_colormap.name}` '
                                     'to plot the legend'))
                ax.set_axis_off()
            else:
                scale_square = np.zeros((256, 256, 3))
                max_g1, max_g2, norm, on_channel = points.metadata['2genes_params']
                V1 = np.linspace(0, max_g1, 256)
                V2 = np.linspace(0, max_g2, 256)
                VS = np.array([V1, V2])
                VS = norm(VS)
                VS[VS<0] = 0
                VS[1<VS] = 1
                scale_square[...,np.where(1-on_channel)[0][0]] = VS[0]
                for axes in np.where(on_channel)[0]:
                    scale_square[...,axes] = VS[1].reshape(-1, 1)
                ax.imshow(scale_square.swapaxes(1, 0), origin='lower')
                recap_g1 = lambda x: x*255/max_g1
                recap_g2 = lambda x: x*255/max_g2
                vals_g1 = np.arange(np.floor(max_g1)+1, dtype=int)
                vals_g2 = np.arange(np.floor(max_g2)+1, dtype=int)
                ax.set_xticks(recap_g1(vals_g1))
                ax.set_yticks(recap_g2(vals_g2))
                ax.set_xticklabels(vals_g1)
                ax.set_yticklabels(vals_g2)
                ax.set_xlabel(points.metadata['2genes'][1])
                ax.set_ylabel(points.metadata['2genes'][0])
            fig.tight_layout()
            plt.show()

    @magicgui(call_button='Adjust contrast',
              min={'widget_type': 'FloatSlider', 'max': 1, 'min': 0, 'label': ''},
              max={'widget_type': 'FloatSlider', 'max': 1, 'min': 0, 'label': ''})
    def adj_int(viewer: Viewer, min: float=0, max: float=1):
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        if points.face_color_mode.upper() != 'COLORMAP':
            return
        if max < min:
            max, min = min, max
        points.face_contrast_limits = (min, max)
        points.refresh()

    @magicgui(call_button='Threshold cells',
              min={"widget_type": "FloatSlider", 'max': 1, 'min': 0, 'label': ''},
              max={"widget_type": "FloatSlider", 'max': 1, 'min': 0, 'label': ''},
              result_widget=True)
    def threshold(viewer: Viewer, min: float=0, max: float=1):
        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return
        if not hasattr(points.features, 'current_view'):
            points.features['current_view'] = points.shown.copy()
        if max < min:
            max, min = min, max
        mask = (points.features['current_view'] &
                 (min<=points.features['gene'])&(points.features['gene']<=max))
        points.shown = (mask)
        points.refresh()
        nb_selected = np.sum(mask)
        overall = np.sum(points.features['current_view'])
        return (f'{nb_selected} cells '
                f'({100*nb_selected/overall:.1f}% of the initial)')

    def get_stats(stat_func, indices, variable_genes):
        if variable_genes:
            data = embryo.anndata
            gene_stat = stat_func(data[indices].X, axis=0)
        else:
            data = embryo.anndata.raw
            gene_stat = stat_func(safe_toarray(data[indices].X), axis=0)

        top_genes = np.argsort(gene_stat)[-3:]
        top_gene_names = data.var_names[top_genes]
        return top_genes, top_gene_names

    @magicgui(call_button='Show gene on Umap',
              gene={'label': 'Choose gene', 'value': 'T'},
              tissues={'label': 'Display tissues umap', 'value': False},
              variable_genes={'label': 'Take only variable genes', 'value': True},
              stats={'widget_type': 'RadioButtons',
                     'label': 'Stat for choosing distributions',
                     'choices': ['Standard Deviation', 'Mean', 'Median'],
                     'value': 'Standard Deviation'},
              result_widget=True)
    def umap_gene(viewer: Viewer, gene: str, variable_genes: bool,
                  tissues: bool, stats: str):
        if not hasattr(embryo.anndata, 'obsm') or not embryo.umap_id in embryo.anndata.obsm_keys():
            return 'No umap coordinates found in the dataset'
        mpl.rcParams['font.size'] = 6
        def show_cells(event):
            if event:
                points.shown = np.zeros_like(points.shown, dtype=bool)
                if sum(len(s.ind) for s in selectors)==0:
                    points.shown = points.features['current_view']
                    [pt.set_alpha(1) for pt in pts]
                    ax_G.set_title(f'Gene: {gene}')
                    for i, (ax_hist, (gene_hist, maxi)) in enumerate(ax_hists):
                        vals = embryo.anndata.raw[:, gene_hist]
                        ax_hist.clear()
                        vals = safe_toarray(vals.X)[points.features['current_view'], 0]
                        hist = ax_hist.hist(vals, bins=50)
                        if i==2:
                            ax_hist.set_xlabel('Gene expression')
                        ax_hist.set_ylabel('#cells')
                        ax_hist.set_title(f'{gene_hist} distribution')
                        ax_hist.set_xlim(0, maxi)
                else:
                    for s in selectors:
                        if len(s.ind)!=0:
                            indices = corres_to_mask[sorted_vals][s.ind]
                            nb_init = np.sum(points.features['current_view'])
                            points.shown[indices] = points.features['current_view'][indices]
                            nb = np.sum(points.shown[indices])
                            ax_G.set_title(f'Gene: {gene} ({nb} cells '
                                           f'({100*nb/nb_init:.1f}% of the initial))')
                            (top_genes,
                             top_gene_names) = get_stats(stat_func, indices, variable_genes)
                            for j, (ax_hist, _) in enumerate(ax_hists):
                                gene_name = top_gene_names[j]
                                gene_id = top_genes[j]
                                vals = embryo.anndata.raw[indices, gene_name]
                                ax_hist.clear()
                                hist = ax_hist.hist(safe_toarray(vals.X)[:, 0], bins=50)
                                if tissues:
                                    ax_hist.set_yticks([])
                                    ax_hist.set_xlabel('Gene expression')
                                    if j==0:
                                        ax_hist.set_ylabel('#cells')
                                if not tissues:
                                    ax_hist.set_ylabel('#cells')
                                    if j==2:
                                        ax_hist.set_xlabel('Gene expression')
                                ax_hist.set_title(f'{gene_name} distribution')
                                ax_hist.set_xlim(0, maximums[gene_id])
                            alpha = np.zeros_like(sorted_vals)+.1
                            alpha[s.ind] = 1
                            [pt.set_alpha(alpha) for pt in pts]
                            s.ind = []
                points.refresh()
        if stats == 'Standard Deviation':
            stat_func = np.std
        elif stats == 'Mean':
            stat_func = np.mean
        elif stats == 'Median':
            stat_func = np.median
        else:
            stat_func = np.max

        points = viewer.layers.selection.active
        if points is None or points.as_layer_data_tuple()[-1]!='points':
            error_points_selection()
            return 'Please select a point layer'
        if not hasattr(points.features, 'current_view'):
            points.features['current_view'] = points.shown.copy()
        cell_list = list(embryo.all_cells)
        if points is None or not gene in embryo.anndata.raw.var_names:
            return f"'{gene}' not found"
        if 0<len(points.selected_data):
            mask = np.zeros_like(points.shown)
            mask[list(points.selected_data)] = True
            points.shown[~mask] = False
            points.refresh()
        else:
            mask = points.shown

        static_canvas = FigureCanvas()
        fig = static_canvas.figure
        static_canvas.toolbar = NavigationToolbar(static_canvas,
                                                  static_canvas.parent())
        fig.set_figwidth(10)
        fig.set_figheight(8)
        ax_hists = []
        if tissues:
            gs0 = GridSpec(3, 1, figure=fig)
            gs00 = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[:2, 0])
            ax_G = fig.add_subplot(gs00[0, 0])
            ax_T = fig.add_subplot(gs00[0, 1])
            gs01 = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[2, 0])
            for i in range(3):
                ax = fig.add_subplot(gs01[0, i])
                ax.set_xlabel('Gene expression')
                if i==0:
                    ax.set_ylabel('#cells')
                ax_hists.append([ax, gene])
        else:
            gs = GridSpec(3, 3, figure=fig)
            ax_G = fig.add_subplot(gs[:, :2])
            for i in range(3):
                ax = fig.add_subplot(gs[i, -1])
                if i==2:
                    ax.set_xlabel('Gene expression')
                ax.set_ylabel('#cells')
                ax_hists.append([ax, gene])

        corres_to_mask = np.where(mask)[0]
        val_g = safe_toarray(embryo.anndata.raw[:, gene].X)[mask, 0]
        maximums = safe_toarray(embryo.anndata.raw.X)[mask].max(axis=0)
        colors = val_g
        pos_cells = embryo.anndata.obsm[embryo.umap_id][mask, :2]
        sorted_vals = np.argsort(val_g)
        min_c, max_c = colors.min(), colors.max()
        colors = (colors-min_c)/(max_c-min_c)
        pts_G = ax_G.scatter(*pos_cells[sorted_vals, :].T, marker='.',
                              c=colors[sorted_vals])
        pts = [pts_G]
        (top_genes,
         top_gene_names) = get_stats(stat_func, points.features['current_view'], variable_genes)
        for j, (ax_hist, _) in enumerate(ax_hists):
            gene_hist = top_gene_names[j]
            gene_id = top_genes[j]
            vals = embryo.anndata.raw[:, gene_hist]
            ax_hist.hist(safe_toarray(vals.X)[mask, 0], bins=50)
            ax_hists[j][1] = (gene_hist, maximums[gene_id])
            ax_hist.set_title(f'{gene_hist} distribution')
            ax_hist.set_xlim(0, maximums[gene_id])
        pts_G.set_edgecolor('none')
        ax_G.set_xticks([])
        ax_G.set_yticks([])
        ax_G.set_xlabel('umap 1')
        ax_G.set_ylabel('umap 2')
        ax_G.set_title(f'Gene: {gene}')
        ax_G.set_aspect('equal')
        fig.tight_layout()
        selectors = []
        selectors.append(SelectFromCollection(ax_G, pts_G))
        if tissues:
            colors_T = [color_map_tissues.get(embryo.tissue[c], [0, 0, 0])
                        for c in corres_to_mask]
            colors_T = np.array(colors_T)
            pts_T = ax_T.scatter(*pos_cells[sorted_vals, :].T, marker='.',
                                 color=colors_T[sorted_vals])
            pts.append(pts_T)
            pts_T.set_edgecolor('none')
            ax_T.set_xticks([])
            ax_T.set_yticks([])
            ax_T.set_ylabel('umap 2')
            ax_T.set_title(f'Tissues')
            tissues_found = set([embryo.tissue[c] for c in corres_to_mask])
            for t in tissues_found:
                ax_T.plot([], [], 'o',
                          color=color_map_tissues.get(t, [0, 0, 0]),
                          label=embryo.corres_tissue.get(t, f'{t}'))
            selectors.append(SelectFromCollection(ax_T, pts_T))
            ax_T.set_aspect('equal')
            ax_T.legend(fontsize='xx-small', frameon=False, shadow=False)

        @magicgui(call_button='Save')
        def saving_fig(viewer: Viewer):
            file_path = QFileDialog.getSaveFileName()
            fig.savefig(file_path[0])
        @magicgui(call_button='Tight Layout')
        def tight_layout(viewer: Viewer):
            fig.tight_layout()
        fig.canvas.mpl_connect("button_release_event", show_cells)
        fig_can = viewer.window.add_dock_widget(static_canvas, name='umap')
        V_box = QWidget()
        V_box.setLayout(QVBoxLayout())
        V_box.layout().addWidget(fig_can)
        V_box.layout().addWidget(static_canvas.toolbar)
        viewer.window.add_dock_widget(V_box)

    @magicgui(call_button='Compute and show surface',
              tissue={'widget_type': 'ComboBox',
                    'label': 'Choose tissue',
                    'choices': all_tissues,
                    'value': all_tissues[0]},
              result_widget=True)
    def show_surf(viewer: Viewer, tissue: str):
        try:
            from pyvista import PolyData
        except Exception as e:
            raise(('pyvista should be install to run that command\n'
                   'Try pip install pyvista to install it'))
        curr_layer = viewer.layers.selection.active
        tissue_to_num = {v:k for k, v in embryo.corres_tissue.items()}
        t_id = tissue_to_num[tissue]
        points = [embryo.pos_3D[c] for c in embryo.cells_from_tissue[t_id]]
        pd = PolyData(points)
        mesh = pd.delaunay_3d().extract_surface()
        face_list = list(mesh.faces.copy())
        face_sizes = {}
        faces = []
        while 0<len(face_list):
            nb_P = face_list.pop(0)
            if not nb_P in face_sizes:
                face_sizes[nb_P] = 0
            face_sizes[nb_P] += 1
            curr_face = []
            for _ in range(nb_P):
                curr_face.append(face_list.pop(0))
            faces.append(curr_face)
        faces = np.array(faces)
        viewer.add_surface((mesh.points, faces),
                           colormap=(color_map_tissues.get(t_id, [0, 0, 0]),),
                           name=tissue, opacity=.6)
        viewer.layers.selection.select_only(curr_layer)

    sel_t = viewer.window.add_dock_widget(select_tissues, name='Tissue selection')
    legend = viewer.window.add_dock_widget(disp_legend, name='Legend')
    show_t = viewer.window.add_dock_widget(show_tissues, name='Tissue colormap')
    g1_cmp = viewer.window.add_dock_widget(show_gene, name='Gene colormap')
    g2_cmp = viewer.window.add_dock_widget(show_two_genes, name='Two genes colormap')
    umap = viewer.window.add_dock_widget(umap_gene, name='umap selection')
    g_th = viewer.window.add_dock_widget(threshold, name='Gene threshold')
    contrast = viewer.window.add_dock_widget(adj_int, name='Contrast')
    cmap = viewer.window.add_dock_widget(apply_cmap, name='Colormap')
    surf = viewer.window.add_dock_widget(show_surf, name='Surface')
    
    tab1 = QTabWidget()
    tab1.addTab(sel_t, sel_t.name)
    tab1.addTab(legend, legend.name)
    tab1.addTab(show_t, show_t.name)
    tab1.addTab(surf, surf.name)

    tab2 = QTabWidget()
    w_box = QWidget()
    w_box.setLayout(QVBoxLayout())
    w_box.layout().addWidget(g1_cmp)
    tab3 = QTabWidget()
    tab3.addTab(g_th, g_th.name)
    tab3.addTab(contrast, contrast.name)
    tab3.addTab(cmap, cmap.name)
    w_box.layout().addWidget(tab3)
    tab2.addTab(w_box, 'Single Value')
    tab2.addTab(g2_cmp, g2_cmp.name)
    tab2.addTab(umap, 'umap')
    
    viewer.window.add_dock_widget(tab1, name='Tissue visualization')
    viewer.window.add_dock_widget(tab2, name='Metric visualization')

    return viewer

class Startsc3D(QWidget):
    def _on_click(self, event):
        data_path = self.h5ad_file.line_edit.value
        tissue_names = self.json_file.line_edit.value
        tissue_names = Path(tissue_names)
        if tissue_names.suffix == '.json':
            with open(tissue_names) as f:
                corres_tissues = json.load(f)
                corres_tissues = {k if isinstance(k, int) else eval(k): v
                                    for k, v in corres_tissues.items()}
        else:
            corres_tissues = {}
        self.embryo = Embryo(data_path, store_anndata=True, corres_tissue=corres_tissues,
                             tissue_id = self.tissue_id.value,
                             pos_reg_id = self.pos_reg_id.value,
                             gene_name_id = self.gene_name_id.value,
                             umap_id = self.umap_id.value)
        self.viewer.window.remove_dock_widget('all')
        display_embryo(self.viewer, self.embryo)
        return

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        self.tissue_id = LineEdit(label='Column name for Tissue id', value='predicted.id')
        self.pos_reg_id = LineEdit(label='Column name for 3D position', value='X_spatial_registered')
        self.gene_name_id = LineEdit(label='Column name for gene names', value='feature_name')
        self.umap_id = LineEdit(label='Column name for umap coordinates', value='X_umap')

        self.h5ad_file = FileEdit(label='h5ad file', value=Path('.').absolute(), filter='*.h5ad')
        self.json_file = FileEdit(label='Tissue names', value=Path('.').absolute(), filter='*.json')
        load_atlas = QPushButton('Load Atlas')
        self.setLayout(QVBoxLayout())
        load = QWidget()
        load.setLayout(QVBoxLayout())
        load.layout().addWidget(self.h5ad_file.native)
        load.layout().addWidget(self.json_file.native)

        params = QWidget()
        params.setLayout(QVBoxLayout())
        params.layout().addWidget(self.tissue_id.native)
        params.layout().addWidget(self.pos_reg_id.native)
        params.layout().addWidget(self.gene_name_id.native)
        tab = QTabWidget()
        tab.addTab(load, 'Loading files')
        tab.addTab(params, 'Parameters')
        self.layout().addWidget(tab)
        self.layout().addWidget(load_atlas)
        load_atlas.clicked.connect(self._on_click)


        # btn = QPushButton("Click me!")
        # btn.clicked.connect(self._on_click)

        # self.setLayout(QHBoxLayout())
        # self.layout().addWidget(self.get_parameters.native)
        # self.layout().addWidget(self.load_file.native)
        # self.viewer.window.add_dock_widget(self)

        
        # self.loading_embryo()

        # tab = QTabWidget()
        # tab.addTab(self.load_file.native, 'Test1')
        # tab.addTab(self.get_parameters.native, 'Test2')
        # self.viewer.window.add_dock_widget(tab, name='Atlas loading')