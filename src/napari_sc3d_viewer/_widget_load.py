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


class LoadAtlas(QWidget):
    """
    Build the initial widget to load the spatial single cell data
    """

    def _load_data(self):
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

        sample_list = self.sample_list.value

        if 0 < len(sample_list):
            delim = None
            if "," in sample_list:
                delim = ","
            elif ";" in sample_list:
                delim = ";"
            sample_list = [s.strip() for s in sample_list.split(delim)]
        else:
            sample_list = None

        self.embryo = Embryo(
            data_path,
            store_anndata=True,
            corres_tissue=corres_tissues,
            tissue_id=self.tissue_id.value,
            pos_reg_id=self.pos_reg_id.value,
            gene_name_id=self.gene_name_id.value,
            umap_id=self.umap_id.value,
            sample_list=sample_list,
        )

    def _on_click_load(self):
        """
        Function to load and call the browsing widget
        """
        # When clicking on the loading button

        # Clearing the viewer and running the viewer plugin
        self._load_data()
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

        # Parameters information widget
        tissue_id_label = widgets.Label(value="Column name for Tissue id")
        self.tissue_id = widgets.LineEdit(value="predicted.id")
        tissue_id = widgets.Container(
            widgets=[tissue_id_label, self.tissue_id], labels=False
        )
        pos_reg_id_label = widgets.Label(value="Column name for 3D position")
        self.pos_reg_id = widgets.LineEdit(value="X_spatial_registered")
        pos_reg_id = widgets.Container(
            widgets=[pos_reg_id_label, self.pos_reg_id], labels=False
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
            widgets=[tissue_id, pos_reg_id, gene_name_id, umap_id],
            labels=False,
        )
        params.native.layout().addStretch(1)

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
        self.sample_list = widgets.LineEdit(value="")
        sample_list = widgets.Container(
            widgets=[sample_list_label, self.sample_list], labels=False
        )
        load = widgets.Container(
            widgets=[h5ad, json, sample_list, self.out_read], labels=False
        )
        load.native.layout().addStretch(1)

        # Tabifying
        tab_atlas = QTabWidget()
        tab_atlas.addTab(load.native, "File paths")
        tab_atlas.addTab(params.native, "Parameters")
        tab_atlas.native = tab_atlas

        # Loading button
        load_atlas = QPushButton("Load Atlas")
        load_atlas.native = load_atlas
        atlas_load_W = widgets.Container(
            widgets=[tab_atlas, load_atlas], labels=False
        )

        # Slight improvement of the layout
        layout = QVBoxLayout()
        layout.addStretch(1)
        self.setLayout(layout)
        self.layout().addWidget(atlas_load_W.native)
        tab_atlas.adjustSize()
        load_atlas.clicked.connect(self._on_click_load)
