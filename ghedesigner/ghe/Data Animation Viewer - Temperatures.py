# standard PyQt5 imports
import sys

from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtGui import QCursor
from PyQt5.QtCore import Qt

from OpenGL_2D_class import gl2D

# the ui created by Designer and pyuic
from DataAnimation_ui import Ui_Dialog

# import the Problem Specific class
# !!!!!!!!! this is the first of three custom lines in this file!!!!!!
from District_system_class import GHEHPSystem

from ghedesigner.ghe.run_simulation import read_data_from_json_file
import time

class main_window(QDialog):
    def __init__(self):
        super(main_window, self).__init__()
        self.ui = Ui_Dialog()
        # setup the GUI
        self.ui.setupUi(self)

        #   !!!!!!!!! this is the second of three custom lines in this file!!!!!!
        # Connect to your custom Animation-ready Class
        self.myAnimatorClass = GHEHPSystem # No parentheses here, this is not an Instance of the class

        #   !!!!!!!!! this is the third of three custom lines in this file!!!!!!
        # Allow a file to be opened and displayed on program startup
        self.defaultFilename = "1-pipe_3ghe-6hp_system_w_pumping_station_input.txt"  # Could be None


        self.myAnimator = None  # a new Animator instance will be created each time a file is read
        # The Animator class must have these three methods:
            # self.myAnimator.DrawPicture()
            # self.myAnimator.PrepareNextAnimationFrameData(current frame,number of frames)
            # self.myAnimator.ProcessFileData(data string) # interprets the data string read from the file
        # After  ProcessFileData() is called, the self.Animator class must have meaningful values in
            # the following drawing size class attributes (data items):
                # self.xmin, self.xmax, self.ymin,self. ymax    - Used to set the window working space
                # self.allowDistortion  - Will circles display as round or elliptical?
            # And the following animation control class attributes:
                # self.numberOfAnimationFrames   - total number of animation frames
                # self.AnimDelayTime  - delay time between frames
                # self.AnimReverse, self.AnimRepeat, self.AnimReset


        # create and setup the GL window object
        self.glwindow1 = None

        self.setupGLWindows()

        # and define any Widget callbacks (buttons, etc) or other necessary setup
        self.assign_widgets()

        # show the GUI
        self.show()

    def DrawingCallback(self):
        # this is what actually draws the picture
        if self.myAnimator is None: return
        self.myAnimator.DrawPicture()  # drawing is done by the DroneCatcher object


    def AnimationCallback(self, frame, nframes):
        # calculations handled by DroneCapture class
        self.myAnimator.PrepareNextAnimationFrameData(frame, nframes)
        self.ui.horizontalSlider_frame.setValue(frame)
        self.ui.Frame_Number.setText(str(frame))
        # the next line is absolutely required for pause, resume, stop, etc !!!
        app.processEvents()
        pass


    def assign_widgets(self):  # callbacks for Widgets on your GUI
        self.ui.pushButton_Exit.clicked.connect(self.ExitApp)
        self.ui.pushButton_Animate.clicked.connect(self.StartAnimation)
        self.ui.pushButton_StopAnimation.clicked.connect(self.StopAnimation)
        self.ui.pushButton_PauseResumeAnimation.clicked.connect(self.PauseResumeAnimation)
        self.ui.horizontalSlider_zoom.valueChanged.connect(self.glZoomSlider)
        self.ui.horizontalSlider_frame.valueChanged.connect(self.glFrameSlider)
        self.ui.pushButton_GetFile.clicked.connect(self.ReadFile)
        self.ui.checkBox_Repeat.stateChanged.connect(self.CheckBoxRepeat)
        self.ui.checkBox_Reverse.stateChanged.connect(self.CheckBoxReverse)


    def ReadFile(self, filename = False):

        if filename is False:
        # get the filename using the OPEN dialog
            filename = QFileDialog.getOpenFileName()[0]
            if len(filename) == 0:
                no_file()
                return
        self.ui.textEdit_filename.setText(filename)
        app.processEvents()
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))

        # Read the file
        f1 = open(filename, 'r')  # open the file for reading
        data = f1.readlines()  # read the entire file as a list of strings
        f1.close()  # close the file  ... very important

        #try:
        self.myAnimator = self.myAnimatorClass()
        anim = self.myAnimator
        anim.read_GHEHPSystem_data(data)
        self.solve()

        self.glwindow1.setViewSize(anim.xmin,anim.xmax,anim.ymin,anim.ymax, anim.allowDistortion)

        self.ui.Drawing_bounds.setText("X: " + str(anim.xmin)+ ", "+ str(anim.xmax)+ ",      Y: " + str(anim.ymin)+ ", "+str(anim.ymax))

        QApplication.restoreOverrideCursor()
        self.glwindow1.glUpdate()
        self.ui.horizontalSlider_frame.setValue(0)
        self.ui.horizontalSlider_frame.setMaximum(self.myAnimator.numberOfAnimationFrames)
        self.ui.Frame_Number.setText(str(0))
        self.ui.horizontalSlider_zoom.setValue(100)
        self.ui.checkBox_Repeat.setChecked(anim.AnimRepeat)
        self.ui.checkBox_Reverse.setChecked(anim.AnimReverse)


# Widget callbacks start here

    def glZoomSlider(self):  # I used a slider to control GL zooming
        zoomval = float((self.ui.horizontalSlider_zoom.value()) / 200 + .5)
        self.glwindow1.glZoom(zoomval)  # set the zoom value
        self.glwindow1.glUpdate()  # update the GL image


    def glFrameSlider(self):  # I used a slider to control manual animation
        frameval = int(self.ui.horizontalSlider_frame.value())
        self.myAnimator.PrepareNextAnimationFrameData(frameval, 120)
        self.ui.Frame_Number.setText(str(frameval))
        self.glwindow1.glUpdate()  # update the GL image
        #self.setAngleSliderAndText()


    def StartAnimation(self):  # a button to start GL Animation
        anim = self.myAnimator
        self.glwindow1.glStartAnimation(self.AnimationCallback, anim.numberOfAnimationFrames,
                                        delaytime= anim.AnimDelayTime, reverseDelayTime = 0.5,
                                        reverse=anim.AnimReverse, repeat=anim.AnimRepeat, reset = anim.AnimReset)


    def StopAnimation(self):  # a button to Stop GL Animati0n
        self.glwindow1.glStopAnimation()

    def PauseResumeAnimation(self):  # a button to Resume GL Animation
        self.glwindow1.glPauseResumeAnimation()

    def CheckBoxRepeat(self):  # used a checkbox to Enable drawing of construction lines
        if self.ui.checkBox_Repeat.isChecked():  # it is on
            self.myAnimator.AnimRepeat = True
        else:  # stop
            self.myAnimator.AnimRepeat = False

        self.glwindow1.glUpdate()

    def CheckBoxReverse(self):  # used a checkbox to Enable drawing of construction lines
        if self.ui.checkBox_Reverse.isChecked():  # it is on
            self.myAnimator.AnimReverse = True
        else:  # stop
            self.myAnimator.AnimReverse = False
        self.glwindow1.glUpdate()

    def ExitApp(self):
        app.exit()

    # Essential if using mouse information with GL
    def eventFilter(self, source, event):  # allow GL to handle Mouse Events
        self.glwindow1.glHandleMouseEvents(event)  # let GL handle the event
        return super(QDialog, self).eventFilter(source, event)

    # Setup OpenGL Drawing and Viewing
    def setupGLWindows(self):  # setup all GL windows
        # send it the   GL Widget     and the drawing Callback function
        self.glwindow1 = gl2D(self.ui.openGLWidget, self.DrawingCallback)

        # set the drawing space:    xmin  xmax  ymin   ymax
        self.glwindow1.setViewSize(0,1,0,1, allowDistortion=False)

        # Optional: Setup GL Mouse Functionality
        self.ui.openGLWidget.installEventFilter(self)  # to read mouse events
        self.ui.openGLWidget.setMouseTracking(True)  # to enable mouse events

        # OPTIONAL: to display the mouse location  - the name of the TextBox
        self.glwindow1.glMouseDisplayTextBox(self.ui.MouseLocation)

    def solve(self):

        fluid, pipe, grout, soil, borehole, sim_params = read_data_from_json_file()
        self.myAnimator.solveSystem(fluid, pipe, grout, soil, borehole, sim_params)
        self.myAnimator.createOutput()
        self.myAnimator.current_frame = 1



def no_file():
    msg = QMessageBox()
    msg.setText('There was no file selected')
    msg.setWindowTitle("No File")
    retval = msg.exec_()
    return None

def bad_file():
    msg = QMessageBox()
    msg.setText('Unable to process the selected file')
    msg.setWindowTitle("Bad File")
    retval = msg.exec_()
    return None


if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    if main_win.defaultFilename is  not  None:  #read the default file
        main_win.ReadFile(main_win.defaultFilename)
    sys.exit(app.exec_())


# def AnimationCallback(frame, nframes):
#     # calculations needed to configure the picture
#     # these could be done here or by calling a class method
#     myAnimator.current_frame = (frame + 1) * 145


