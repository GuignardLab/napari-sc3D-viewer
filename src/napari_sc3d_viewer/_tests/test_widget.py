from napari_sc3d_viewer import LoadAtlas, RegisterSc3D
import numpy as np

try:
    import pyvista

    pyvista_install = True
except:
    pyvista_install = False

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_example_q_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()

    # create our widget, passing in the viewer
    my_widget = LoadAtlas(viewer, show=False)

    # call our widget method
    my_widget.h5ad_file.value = "test_data/data_test.h5ad"
    my_widget.json_file.value = "test_data/wrong_format.json"
    my_widget.sample_list.value = ""
    displayed_embryo = my_widget._on_click_load()
    my_widget.json_file.value = "test_data/corresptissues.json"
    displayed_embryo = my_widget._on_click_load()
    p = displayed_embryo.viewer.layers.selection.active

    for points in [None, p]:
        displayed_embryo.viewer.layers.selection.active = points

        displayed_embryo.select_tissues_choices.value = [
            "Pharyngeal arch",
            "Somites",
            "Heart",
            "Presomitic mesoderm (PSM)",
        ]
        displayed_embryo.select_tissues()
        if points:
            assert (
                np.unique(points.face_color[points.shown], axis=0).shape[0]
                == 4
            )

        displayed_embryo.cmap.value = "magma"
        displayed_embryo.apply_cmap()

        displayed_embryo.disp_legend()

        displayed_embryo.metric.value = "Gene"
        displayed_embryo.gene.value = "T"
        displayed_embryo.show_gene()
        if points:
            assert (
                np.unique(points.face_color[points.shown], axis=0).shape[0]
                == 43
            )

        displayed_embryo.metric.value = "Gene"
        displayed_embryo.gene.value = "b"
        displayed_embryo.show_gene()
        if points:
            assert displayed_embryo.gene_output.value == "Gene 'b' not found"

        displayed_embryo.show_tissues()
        if points:
            assert (
                np.unique(points.face_color[points.shown], axis=0).shape[0]
                == 4
            )

        displayed_embryo.metric.value = "xcoord"
        displayed_embryo.gene.value = ""
        displayed_embryo.show_gene()
        if points:
            assert (
                np.unique(points.face_color[points.shown], axis=0).shape[0]
                == 103
            )

        displayed_embryo.threshold_low.value = 0.1
        displayed_embryo.threshold_high.value = 0.9
        displayed_embryo.threshold()
        if points:
            assert points.shown.sum() == 91

        displayed_embryo.threshold_low.value = 0.0
        displayed_embryo.threshold_high.value = 1.0
        displayed_embryo.threshold()
        if points:
            assert points.shown.sum() == 103

        displayed_embryo.adj_int_low.value = 0.2
        displayed_embryo.adj_int_high.value = 0.8
        displayed_embryo.adj_int()

        displayed_embryo.cmap.value = "turbo"
        displayed_embryo.apply_cmap()

        displayed_embryo.disp_legend()

        displayed_embryo.gene1.value = "T"
        displayed_embryo.gene2.value = "Sox2"
        displayed_embryo.threhold_low_2g.value = 2
        displayed_embryo.threhold_high_2g.value = 80
        displayed_embryo.main_bi_color.value = "Red"
        displayed_embryo.show_two_genes()
        if points:
            assert (
                np.unique(points.face_color[points.shown], axis=0).shape[0]
                == 62
            )

        displayed_embryo.disp_legend()

        displayed_embryo.umap_selec.gene.value = "T"
        displayed_embryo.umap_selec.tissues.value = False
        displayed_embryo.umap_selec.tissues.stats = "Mean"
        displayed_embryo.umap_selec.run()

        displayed_embryo.umap_selec.gene.value = "Sox2"
        displayed_embryo.umap_selec.tissues.value = True
        displayed_embryo.umap_selec.tissues.stats = "Mean"
        displayed_embryo.umap_selec.run()

        displayed_embryo.umap_selec.show_cells(False)
        displayed_embryo.umap_selec.show_cells(True)

        if pyvista_install:
            displayed_embryo.surf_threshold.value = 0
            displayed_embryo.show_surf()
            displayed_embryo.surf_threshold.value = 5
            displayed_embryo.show_surf()

        my_widget = LoadAtlas(viewer, show=False)

    # create our widget, passing in the viewer
    my_widget = RegisterSc3D(viewer, show=False)

    # call our widget method
    my_widget.h5ad_file.value = "test_data/data_test.h5ad"
    my_widget.json_file.value = "test_data/corresptissues.json"
    displayed_embryo = my_widget._on_click_sc3D()

    my_widget = RegisterSc3D(viewer, show=False)

    # call our widget method
    my_widget.h5ad_file.value = "test_data/data_test.h5ad"
    my_widget.json_file.value = "test_data/corresptissues.json"
    displayed_embryo = my_widget._on_click_PASTE()
