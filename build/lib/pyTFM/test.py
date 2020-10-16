from PyQt5.QtWidgets import QWidget,QApplication
from PyQt5.QtCore import pyqtSignal,QThread
from PyQt5 import QtWidgets
import sys

class Main(QWidget):

    def __init__(self):
        super().__init__()
        self.button_start = QtWidgets.QPushButton("start")
        self.button_start.clicked.connect(self.StartButtonEvent)
        self.layout = QtWidgets.QGridLayout(self)
        self.layout.addWidget(self.button_start, 0, 0)
        #self.setStyle(QtWidgets.QStyleFactory.create("Windows"))
    def StartButtonEvent(self):
        self.test = ExecuteThread()
        self.test.start()
        self.test.finished.connect(self.thread_finished)
        self.test.my_signal.connect(self.my_event)

    def thread_finished(self):
        # gets executed if thread finished
        pass

    def my_event(self):
        # gets executed on my_signal
        pass


class ExecuteThread(QThread):
    my_signal = pyqtSignal()

    def run(self):
        # do something here
        self.my_signal.emit()
        pass
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle('fusion')
    mainWin = Main()
    mainWin.show()
    sys.exit( app.exec_() )