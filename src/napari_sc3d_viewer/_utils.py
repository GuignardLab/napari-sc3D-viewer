from qtpy.QtWidgets import QMessageBox
import numpy as np


def safe_toarray(tab):
    """
    Returns an ndarray from a ndarray or a sparse matrix
    """
    if not isinstance(tab, np.ndarray):
        return tab.toarray()
    else:
        return tab


def error_json_format(show=True):
    """
    Print a message error in a box
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("json file format error")
    msg.setInformativeText(
        (
            "The provided json file did not have the expected format\n"
            "Please refer to the documentation"
        )
    )
    msg.setWindowTitle("json file format error")
    if show:
        msg.exec_()


def error_points_selection(show=True):
    """
    Print a message error in a box
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Point cloud selection error")
    msg.setInformativeText(
        (
            "Please select an adequate point cloud\n"
            "You can select point clouds on the left hand side of the viewer"
        )
    )
    msg.setWindowTitle("Point cloud selection error")
    if show:
        msg.exec_()
