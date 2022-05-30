from napari_sc3d_viewer import Startsc3D

to_test = ['Select tissues',
           'Show gene/metric',
           'Show gene',
           'Color according to cell types',
           'Adjust contrast',
           'Threshold cells',
           'Show gene on Umap',
           'Save',
           'Tight Layout',
           'Compute and show surface']

# make_napari_viewer is a pytest fixture that returns a napari viewer object
# capsys is a pytest fixture that captures stdout and stderr output streams
def test_example_q_widget(make_napari_viewer, capsys):
    # make viewer and add an image layer using our fixture
    viewer = make_napari_viewer()

    # create our widget, passing in the viewer
    my_widget = Startsc3D(viewer)

    # call our widget method
    my_widget.h5ad_file.value = 'test_data/data_test.h5ad'
    my_widget.json_file.value = 'test_data/corresptissues.json'
    out_widget = my_widget._on_click()
    t1 = out_widget.tab1
    t2 = out_widget.tab2
    
    to_treat = t1.children() + t2.children()
    # Basically clicking everywhere that is a button ...
    while 0<len(to_treat):
        current = to_treat.pop()
        if hasattr(current, 'click'):
            if current.text() in to_test:
                current.click()
            # current.click()
        to_treat += current.children()