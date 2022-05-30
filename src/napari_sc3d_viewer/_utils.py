from qtpy.QtWidgets import QMessageBox
import numpy as np

def safe_toarray(tab):
    if not isinstance(tab, np.ndarray):
        return tab.toarray()
    else:
        return tab

def error_points_selection():
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText('Point cloud selection error')
    msg.setInformativeText(('Please select an adequate point cloud\n'
                            'You can select point clouds on the left hand side of the viewer'))
    msg.setWindowTitle('Point cloud selection error')
    msg.exec_()