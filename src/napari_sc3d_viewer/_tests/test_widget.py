from napari_sc3d_viewer import Startsc3D
import numpy as np
import inspect

try:
    import pyvista
    pyvista_loaded = True
except:
    pyvista_loaded = False

def run_all_functions(displayed_embryo):
    for name, function in inspect.getmembers(displayed_embryo):
        if callable(function) and not name.startswith('__'):
            if 'show' in inspect.getargspec(function).args:
                function(show=False)
            else:
                function()

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_example_q_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()

    # create our widget, passing in the viewer
    my_widget = Startsc3D(viewer, show=False)

    # call our widget method
    my_widget.h5ad_file.value = 'test_data/data_test.h5ad'
    my_widget.json_file.value = 'test_data/corresptissues.json'
    displayed_embryo = my_widget._on_click()
    # run_all_functions(displayed_embryo)
    points = displayed_embryo.viewer.layers.selection.active

    displayed_embryo.select_tissues_choices.value=['Pharyngeal arch',
                                                   'Somites', 'Heart',
                                                   'Presomitic mesoderm (PSM)']
    displayed_embryo.select_tissues()
    assert np.unique(points.face_color[points.shown], axis=0).shape[0] == 4

    displayed_embryo.disp_legend()

    displayed_embryo.metric.value = 'Gene'
    displayed_embryo.gene.value = 'T'
    displayed_embryo.show_gene()
    assert np.unique(points.face_color[points.shown], axis=0).shape[0] == 43

    displayed_embryo.metric.value = 'Gene'
    displayed_embryo.gene.value = 'b'
    displayed_embryo.show_gene()
    assert displayed_embryo.gene_output.value == "Gene 'b' not found"

    displayed_embryo.show_tissues()
    assert np.unique(points.face_color[points.shown], axis=0).shape[0] == 4

    displayed_embryo.metric.value = 'xcoord'
    displayed_embryo.gene.value = ''
    displayed_embryo.show_gene()
    assert np.unique(points.face_color[points.shown], axis=0).shape[0] == 103

    displayed_embryo.threshold_low.value = .1
    displayed_embryo.threshold_high.value = .9
    displayed_embryo.threshold()
    assert points.shown.sum() == 91

    displayed_embryo.threshold_low.value = .0
    displayed_embryo.threshold_high.value = 1.
    displayed_embryo.threshold()
    assert points.shown.sum() == 103

    displayed_embryo.gene1.value = 'T'
    displayed_embryo.gene2.value = 'Sox2'
    displayed_embryo.threhold_low_2g.value = 2
    displayed_embryo.threhold_high_2g.value = 80
    displayed_embryo.main_bi_color.value = 'Red'
    displayed_embryo.show_two_genes()
    assert np.unique(points.face_color[points.shown], axis=0).shape[0] == 62

    displayed_embryo.umap_selec.gene.value = 'T'
    displayed_embryo.umap_selec.tissues.value = False
    displayed_embryo.umap_selec.tissues.stats = 'Mean'
    displayed_embryo.umap_selec.run()

    displayed_embryo.umap_selec.gene.value = 'Sox2'
    displayed_embryo.umap_selec.tissues.value = True
    displayed_embryo.umap_selec.tissues.stats = 'Mean'
    displayed_embryo.umap_selec.run()

    if pyvista_loaded:
        displayed_embryo.surf_threshold.value = 0
        displayed_embryo.show_surf()
        displayed_embryo.surf_threshold.value = 5
        displayed_embryo.show_surf()