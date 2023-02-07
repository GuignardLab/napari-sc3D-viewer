"""
This file is subject to the terms and conditions defined in
file 'LICENCE', which is part of this source code package.
Author: Leo Guignard (leo.guignard...@AT@...univ-amu.fr)
"""
import json
from sc3D import Embryo
from qtpy.QtWidgets import QWidget, QVBoxLayout, QTabWidget, QPushButton
from magicgui import widgets
from ._display_embryo import DisplayEmbryo
from ._utils import error_json_format
from pathlib import Path


class RegisterSc3D(QWidget):
    """
    Build the initial widget to load the spatial single cell data
    """

    def _load_data(self, weights, tissue_ignore):
        data_path = Path(self.h5ad_file.line_edit.value)
        tissue_names = Path(self.json_file.line_edit.value)
        if not data_path.exists() or not data_path.suffix in [
            ".h5ad",
            ".h5",
            ".csv",
        ]:
            self.out_read.value = "file not found:\n{}".format(data_path)
            return
        if tissue_names.suffix == ".json" and tissue_names.exists():
            with open(tissue_names, "r") as f:
                corres_tissues = json.load(f)
                try:
                    corres_tissues = {
                        k if isinstance(k, int) else eval(k): v
                        for k, v in corres_tissues.items()
                    }
                except Exception as e:
                    error_json_format(show=self.show)
        else:
            corres_tissues = {}
        sample_list = self.sample_list_value.value
        if 0 < len(sample_list):
            delim = None
            if "," in sample_list:
                delim = ","
            elif ";" in sample_list:
                delim = ";"
            sample_list = [s.strip() for s in sample_list.split(delim)]
        else:
            sample_list = None

        # Loading embryo as a sc3D object
        if weights:
            tissue_weight = self.tissue_weight.value
            try:
                tissue_weight = eval(tissue_weight)
            except Exception as e:
                print("Could not read the tissue weights")
                tissue_weight = {}
            if not isinstance(tissue_weight, dict):
                print("Could not read the tissue weights")
        else:
            tissue_weight = {}
        if tissue_ignore:
            delim = None
            tissues_to_ignore = self.tissues_to_ignore.value
            if "," in tissues_to_ignore:
                delim = ","
            elif ";" in tissues_to_ignore:
                delim = ";"
            tissues_to_ignore = [
                eval(t.strip()) for t in tissues_to_ignore.split(delim)
            ]
        else:
            tissues_to_ignore = []
        self.embryo = Embryo(
            data_path,
            store_anndata=True,
            corres_tissue=corres_tissues,
            tissue_id=self.tissue_id.value,
            gene_name_id=self.gene_name_id.value,
            umap_id=self.umap_id.value,
            sample_list=sample_list,
            pos_id=self.pos_id.value,
            tissues_to_ignore=tissues_to_ignore,
            tissue_weight=tissue_weight,
        )

    def _on_click_sc3D(self):
        self._load_data(weights=True, tissue_ignore=True)
        try:
            self.embryo.registration_3d(th_d=self.th_d.value)
        except:
            self.out_read.value = "The registration failed :/"
        # Clearing the viewer and running the viewer plugin
        self.viewer.window.remove_dock_widget("all")
        return DisplayEmbryo(self.viewer, self.embryo, show=self.show)

    def _on_click_PASTE(self):
        self._load_data(weights=False, tissue_ignore=True)
        try:
            self.embryo.registration_3d(
                th_d=self.th_d.value,
                method="paste",
                min_counts_cells=self.min_counts_cells.value,
                min_counts_genes=self.min_counts_genes.value,
                alpha=self.alpha.value,
                pre_registration=self.pre_reg.value,
            )
        except:
            self.out_read.value = "The registration failed :/"
        # Clearing the viewer and running the viewer plugin
        self.viewer.window.remove_dock_widget("all")
        return DisplayEmbryo(self.viewer, self.embryo, show=self.show)

    def __init__(self, napari_viewer, *, show=False):
        """
        Build the containers for the loading widget

        Args:
            napari_viewer (napari.Viewer): the parent napari viewer
            show (bool): a parameter for testing the plugin since some function
                wait for user input
        """
        super().__init__()
        self.viewer = napari_viewer
        self.show = show

        # File path widget
        h5ad_label = widgets.Label(value="h5ad file")
        self.h5ad_file = widgets.FileEdit(
            value=Path(".").absolute(), filter="*.h5*"
        )
        h5ad = widgets.Container(
            widgets=[h5ad_label, self.h5ad_file], labels=False
        )
        json_label = widgets.Label(value="Tissue names")
        self.json_file = widgets.FileEdit(
            value=Path(".").absolute(), filter="*.json"
        )
        self.out_read = widgets.Label(value="")
        json = widgets.Container(
            widgets=[json_label, self.json_file], labels=False
        )
        sample_list_label = widgets.Label(
            value="List of samples (if multiple h5 files)"
        )
        self.sample_list_value = widgets.LineEdit(value="")
        sample_list = widgets.Container(
            widgets=[sample_list_label, self.sample_list_value], labels=False
        )
        load = widgets.Container(
            widgets=[h5ad, json, sample_list, self.out_read], labels=False
        )
        load.native.layout().addStretch(1)

        # Registration widget
        tissues_to_ignore_label = widgets.Label(
            value="Tissues to ignore\nfor registration"
        )
        self.tissues_to_ignore = widgets.LineEdit(
            value=""
        )  # ('13, 15, 16, 22, 27,'
        #      '29, 32, 36, 40, 41'))
        tissues_to_ignore = widgets.Container(
            widgets=[tissues_to_ignore_label, self.tissues_to_ignore],
            labels=False,
        )
        nb_CS_begin_ignore_label = widgets.Label(
            value="Number of slices to ignore\nfrom the start"
        )
        self.nb_CS_begin_ignore = widgets.SpinBox(value=0, min=0, max=50)
        nb_CS_begin_ignore = widgets.Container(
            widgets=[nb_CS_begin_ignore_label, self.nb_CS_begin_ignore],
            labels=False,
        )
        nb_CS_end_ignore_label = widgets.Label(
            value="Number of slices to ignore\nfrom the end"
        )
        self.nb_CS_end_ignore = widgets.SpinBox(value=0, min=0, max=50)
        nb_CS_end_ignore = widgets.Container(
            widgets=[nb_CS_end_ignore_label, self.nb_CS_end_ignore],
            labels=False,
        )
        xy_resolution_label = widgets.Label(value="XY Resolution in µm")
        self.xy_resolution = widgets.FloatSpinBox(value=0.6)
        xy_resolution = widgets.Container(
            widgets=[xy_resolution_label, self.xy_resolution], labels=False
        )

        global_params_reg = widgets.Container(
            widgets=[
                tissues_to_ignore,
                nb_CS_begin_ignore,
                nb_CS_end_ignore,
                xy_resolution,
                self.out_read,
            ],
            labels=False,
        )
        global_params_reg.native.layout().addStretch(1)

        # sc3D parameters
        tissue_weight_label = widgets.Label(
            value="Weights for tissues\nfor the registration"
        )
        self.tissue_weight = widgets.LineEdit(value="{21:2000, 18:2000}")
        tissue_weight = widgets.Container(
            widgets=[tissue_weight_label, self.tissue_weight], labels=False
        )
        th_d_label = widgets.Label(value="Max distance between paired beads")
        self.th_d = widgets.FloatSpinBox(value=150)
        th_d = widgets.Container(widgets=[th_d_label, self.th_d], labels=False)
        register_sc3D = QPushButton("Register with sc3D")
        register_sc3D.native = register_sc3D
        register_sc3D.name = "register_sc3D"
        sc3D_tab = widgets.Container(
            widgets=[tissue_weight, th_d, register_sc3D], labels=False
        )
        sc3D_tab.native.layout().addStretch(1)
        register_sc3D.clicked.connect(self._on_click_sc3D)

        # PASTE parameters
        min_counts_genes_label = widgets.Label(
            value="Minimum count for\ngene filtering"
        )
        self.min_counts_genes = widgets.SpinBox(value=15)
        min_counts_genes = widgets.Container(
            widgets=[min_counts_genes_label, self.min_counts_genes],
            labels=False,
        )

        min_counts_cells_label = widgets.Label(
            value="Minimum count for\ncell filtering"
        )
        self.min_counts_cells = widgets.SpinBox(value=100)
        min_counts_cells = widgets.Container(
            widgets=[min_counts_cells_label, self.min_counts_cells],
            labels=False,
        )

        work_with_raw_label = widgets.Label(value="Use raw data")
        self.work_with_raw = widgets.CheckBox(value=True)
        work_with_raw = widgets.Container(
            widgets=[work_with_raw_label, self.work_with_raw],
            labels=False,
            layout="horizontal",
        )

        alpha_label = widgets.Label(value="Alpha")
        self.alpha = widgets.FloatSpinBox(value=0.1, min=0, max=1)
        alpha = widgets.Container(
            widgets=[alpha_label, self.alpha], labels=False
        )

        pre_reg_label = widgets.Label(
            value="Heuristics based pre-registration"
        )
        self.pre_reg = widgets.CheckBox(value=True)
        pre_reg = widgets.Container(
            widgets=[pre_reg_label, self.pre_reg],
            labels=False,
            layout="horizontal",
        )

        register_paste = QPushButton("Register with PASTE")
        register_paste.native = register_paste
        register_paste.name = "register_paste"
        paste_tab = widgets.Container(
            widgets=[
                min_counts_genes,
                min_counts_cells,
                work_with_raw,
                alpha,
                pre_reg,
                register_paste,
                self.out_read,
            ],
            labels=False,
        )
        paste_tab.native.layout().addStretch(1)
        register_paste.clicked.connect(self._on_click_PASTE)

        params_reg = QTabWidget()
        params_reg.addTab(global_params_reg.native, "Global reg. params.")
        params_reg.addTab(sc3D_tab.native, "sc3D")
        params_reg.addTab(paste_tab.native, "PASTE")

        # Parameters information widget
        tissue_id_label = widgets.Label(value="Column name for Tissue id")
        self.tissue_id = widgets.LineEdit(value="predicted.id")
        tissue_id = widgets.Container(
            widgets=[tissue_id_label, self.tissue_id], labels=False
        )
        pos_id_label = widgets.Label(value="Column name for 2D position")
        self.pos_id = widgets.LineEdit(value="X_spatial")
        pos_id = widgets.Container(
            widgets=[pos_id_label, self.pos_id], labels=False
        )
        gene_name_id_label = widgets.Label(value="Column name for gene names")
        self.gene_name_id = widgets.LineEdit(value="feature_name")
        gene_name_id = widgets.Container(
            widgets=[gene_name_id_label, self.gene_name_id], labels=False
        )
        umap_id_label = widgets.Label(value="Column name for umap coordinates")
        self.umap_id = widgets.LineEdit(value="X_umap")
        umap_id = widgets.Container(
            widgets=[umap_id_label, self.umap_id], labels=False
        )
        params = widgets.Container(
            widgets=[tissue_id, pos_id, gene_name_id, umap_id], labels=False
        )
        params.native.layout().addStretch(1)

        tab_reg = QTabWidget()
        tab_reg.addTab(load.native, "Paths")
        tab_reg.addTab(params.native, "Column labels")
        tab_reg.addTab(params_reg, "Reg. params.")

        # Slight improvement of the layout
        layout = QVBoxLayout()
        layout.addStretch(1)
        self.setLayout(layout)
        self.layout().addWidget(tab_reg)
        tab_reg.adjustSize()
