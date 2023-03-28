import signal

from PyQt5.QtGui import QIcon

from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QFileDialog, QStyleFactory

# ------------
from PyQt5.QtCore import QTimer, QUrl
from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
# from PyQt5.QtWidgets import QApplication, QMainWindow, QStyleFactory, QVBoxLayout, QFileDialog


import time
import matplotlib

from filtering import ProcessFunction


matplotlib.use("Qt5Agg")  # 声明使用QT5
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import os




# 继承父类QMainWindows和Ui_MainWindow
class MySignal(QMainWindow, signal.Ui_MainWindow):
    # 重写父类的init函数
    def __init__(self):
        # 方法中调用第一个父类的init函数
        super(MySignal, self).__init__()

        self.setupUi(self)

        # 设置窗口的主标题
        self.setWindowTitle("Signal Filter")

        # 设置窗口的默认尺寸
        self.resize(1200, 800)

        # # 设置窗口状态栏
        # self.status = self.statusBar()
        # # 显示状态信息，数字为时间
        # self.status.showMessage('hello', 5000)

        # 设置窗口图标
        self.setWindowIcon(QIcon('logo.ico'))

        # 浏览文件，获取音频文件名，及路径
        # self.loadAudioPath()
        # self.pushButton1.clicked.connect(self.loadAudio)

    # def openFile(self):
            # layout = QVBoxLayout()
            # self.pushButton1.clicked.connect(self.loadAudio)

            # layout.addWidget(self.pushButton1)
            # self.imageLabel = QLabel()
            # layout.addWidget(self.imageLabel)
            # self.setLayout(layout)
        # 默认打开第一个页面
        self.stackedWidget.setCurrentIndex(0)
        # pushButton绑定堆栈窗口
        # 通过 connect 建立信号/槽连接，点击按钮事件发射 triggered 信号，执行相应的子程序 click_pushButton
        self.pushButton3.clicked.connect(self.click_pushButton_3)  # 点击 pushButton1 触发
        self.pushButton4.clicked.connect(self.click_pushButton_4)  # 点击 pushButton2 触发
        self.pushButton5.clicked.connect(self.click_pushButton_5)  # 点击 pushButton3 触发

        # --------------
        self.Process = ProcessFunction()#process对象包含了所有的信号处理函数及其画图
        self.saveDatepath_IIR = os.getcwd().replace('\\','/')+"/ProcessedSignal/sweep.wav"##处理过后的数据保存位置
        self.saveDatepath_FIR = os.getcwd().replace('\\','/')+"/ProcessedSignal/sweepfir.wav"  ##处理过后的数据保存位置

        print(self.saveDatepath_IIR)
        #self.saveDatepath_IIR = os.getcwd()+"/ProcessedSignal/sweep.wav"##处理过后的数据保存位置
        #self.saveDatepath_FIR = os.getcwd()+"/sweepfir.wav"  ##处理过后的数据保存位置
        #print(self.saveDatepath_IIR)

        self.progressBar.setValue(0)#进度条初始化为0
        # #***************标志位的初始化*******************
        self.process_flag=0#处理完毕标志位
        self.isPlay = 0#播放器播放标志位
        self.isPlay_IIR = 0  # 播放器播放标志位
        self.isPlay_FIR = 0

        # **************播放器的设定**********************

        self.player = QMediaPlayer(self)  # 这个播放器是播放原声的
        self.player_IIR = QMediaPlayer(self)  # 定义两个对象出来，这个负责播放处理过后的
        self.player_FIR = QMediaPlayer(self)

        # 美化播放进度条
        self.horizontalSlider.sliderMoved[int].connect(
            lambda: self.player.setPosition(self.horizontalSlider.value()))
        self.horizontalSlider.setStyle(QStyleFactory.create('Fusion'))

        #analyse
        self.horizontalSlider_2.sliderMoved[int].connect(
            lambda: self.player.setPosition(self.horizontalSlider_2.value()))
        self.horizontalSlider_2.setStyle(QStyleFactory.create('Fusion'))
        ##IIR
        self.horizontalSlider_4.sliderMoved[int].connect(
            lambda: self.player_IIR.setPosition(self.horizontalSlider_4.value()))
        self.horizontalSlider_4.setStyle(QStyleFactory.create('Fusion'))
        ##FIR
        self.horizontalSlider_5.sliderMoved[int].connect(
            lambda: self.player_FIR.setPosition(self.horizontalSlider_3.value()))
        self.horizontalSlider_5.setStyle(QStyleFactory.create('Fusion'))

        self.timer = QTimer(self)
        self.timer.start(1000)  ##定时器设定为1s，超时过后链接到playRefresh刷新页面
        self.timer.timeout.connect(self.playRefresh)  ##

        # **************菜单栏的事件绑定*******************
        self.pushButton1.clicked.connect(self.onFileOpen)  ##菜单栏的打开文件
        self.pushButton2.clicked.connect(self.close)  # 菜单栏的退出action
        self.Timelayout_()  ##时间域的四个图窗布局
        self.Iirlayout_()  ##IIR设计界面的四个图窗布局
        self.Firlayout_()  ##FIR设计界面的四个图窗布局


        #**************第一个界面的事件绑定配置*************
        self.dial.setValue(20)#默认音量大小为20
        self.dial_2.setValue(20)
        self.dial.valueChanged.connect(self.changeVoice0)##音量圆盘控制事件绑定,如果值被改变就调起事件
        self.dial_2.valueChanged.connect(self.changeVoice01)
        self.pushButton_3.clicked.connect(self.Analyse_btn_start)  # 给pushButton_3添加一个点击事件，分析音频
        self.pushButton.clicked.connect(self.palyMusic)  # 播放原音频
        self.pushButton_6.clicked.connect(self.palyMusic)  # 播放分析音频

        #**************第二个界面的事件绑定配置*************
        # self.dial_4.setValue(20)  # 默认音量大小为20
        # self.dial_4.valueChanged.connect(self.changeVoice1)  ##音量圆盘控制事件绑定,如果值被改变就调起事件
        self.dial_4.setValue(20)  # 默认音量大小为20
        self.dial_4.valueChanged.connect(self.changeVoice)
        self.pushButton_7.clicked.connect(self.desigenIIR)#点击开始设计IIR滤波器按钮之后，调用函数
        self.pushButton_8.clicked.connect(self.applyIIR)#点击应用滤波器
        # self.pushButton_4.clicked.connect(self.palyMusic)
        self.pushButton_4.clicked.connect(self.playIIRaudio)

        #**************第三个界面的事件绑定配置*************
        # self.dial_10.setValue(20)  # 默认音量大小为20
        # self.dial_10.valueChanged.connect(self.changeVoice2)  ##音量圆盘控制事件绑定,如果值被改变就调起事件
        self.dial_5.setValue(20)  # 默认音量大小为20
        self.dial_5.valueChanged.connect(self.changeVoice)
        self.pushButton_17.clicked.connect(self.designFIR)#点击开始设计IIR滤波器按钮之后，调用函数
        self.pushButton_16.clicked.connect(self.applyFIR)#点击应用滤波器
        # self.pushButton_16.clicked.connect(self.palyMusic)
        self.pushButton_5.clicked.connect(self.playFIRaudio)



    def click_pushButton_3(self):  # 点击 pushButton_1 触发
        # self.textEdit.append("当前动作：click_pushButton_1")
        # self.textEdit.append("选择堆叠布局页面：page_0")
        self.stackedWidget.setCurrentIndex(1)  # 打开 stackedWidget > page_0
        # self.label_1.setPixmap(QtGui.QPixmap("../image/fractal01.png"))
        return

    def click_pushButton_4(self):  # 点击 pushButton_1 触发
        # self.textEdit.append("当前动作：click_pushButton_1")
        # self.textEdit.append("选择堆叠布局页面：page_0")
        self.stackedWidget.setCurrentIndex(2)  # 打开 stackedWidget > page_0
        # self.label_1.setPixmap(QtGui.QPixmap("../image/fractal01.png"))
        return

    def click_pushButton_5(self):  # 点击 pushButton_1 触发
        # self.textEdit.append("当前动作：click_pushButton_1")
        # self.textEdit.append("选择堆叠布局页面：page_0")
        self.stackedWidget.setCurrentIndex(3)  # 打开 stackedWidget > page_0
        # self.label_1.setPixmap(QtGui.QPixmap("../image/fractal01.png"))
        return


    #
    # def loadAudio(self):
    #     fname, _ = QFileDialog.getOpenFileName(self, '打开文件', '.', "音频文件 (*.wav *.mp4)")
    #     self.play(self, fname)
    #
    # def play(self, path):
    #     self.sound = QSound(path, self)
    #     self.play_btn = QPushButton('Play Sound', self)
    #     self.play_btn.clicked.connect(self.sound.play)

    def Timelayout_(self):
        self.fig5 = plt.figure()
        self.canvas5 = FigureCanvas(self.fig5)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas5)
        self.graphicsView.setLayout(layout)  # 设置好布局之后调用函数

        self.fig6 = plt.figure()
        self.canvas6 = FigureCanvas(self.fig6)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas6)
        self.graphicsView_2.setLayout(layout)  # 设置好布局之后调用函数

        self.fig7 = plt.Figure()
        self.canvas7 = FigureCanvas(self.fig7)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas7)
        self.graphicsView_3.setLayout(layout)  # 设置好布局之后调用函数

        self.fig8 = plt.Figure()
        self.canvas8 = FigureCanvas(self.fig8)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas8)
        self.graphicsView_4.setLayout(layout)  # 设置好布局之后调用函数

    def Iirlayout_(self):
        self.fig1 = plt.figure()
        self.canvas1 = FigureCanvas(self.fig1)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas1)
        self.graphicsView_5.setLayout(layout)  # 设置好布局之后调用函数

        self.fig2 = plt.figure()
        self.canvas2 = FigureCanvas(self.fig2)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas2)
        self.graphicsView_7.setLayout(layout)  # 设置好布局之后调用函数

        self.fig3 = plt.Figure()
        self.canvas3 = FigureCanvas(self.fig3)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas3)
        self.graphicsView_6.setLayout(layout)  # 设置好布局之后调用函数

        self.fig4 = plt.Figure()
        self.canvas4 = FigureCanvas(self.fig4)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas4)
        self.graphicsView_8.setLayout(layout)  # 设置好布局之后调用函数

    def Firlayout_(self):
        self.fig25 = plt.figure()
        self.canvas25 = FigureCanvas(self.fig25)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas25)
        self.graphicsView_9.setLayout(layout)  # 设置好布局之后调用函数

        self.fig26 = plt.figure()
        self.canvas26 = FigureCanvas(self.fig26)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas26)
        self.graphicsView_10.setLayout(layout)  # 设置好布局之后调用函数

        self.fig27 = plt.Figure()
        self.canvas27 = FigureCanvas(self.fig27)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas27)
        self.graphicsView_11.setLayout(layout)  # 设置好布局之后调用函数

        self.fig28 = plt.Figure()
        self.canvas28 = FigureCanvas(self.fig28)
        layout = QVBoxLayout()  # 垂直布局
        layout.addWidget(self.canvas28)
        self.graphicsView_12.setLayout(layout)  # 设置好布局之后调用函数



    ###############对应的一些触发方法######################
    def onFileOpen(self): ##打开文件
        self.path, _ = QFileDialog.getOpenFileName(self, '打开文件', '', '音乐文件 (*.wav)')

        if self.path:##选中文件之后就选中了需要播放的音乐，并同时显示出来
            self.isPlay=0#每次打开文件的时候就需要暂停播放，无论是否在播放与否
            self.isPlay_IIR=0

            self.player.pause()
            self.player_IIR.pause()

            self.player.setMedia(QMediaContent(QUrl(self.path)))  ##选中需要播放的音乐
            # 绑定到菜单栏的进度条
            self.horizontalSlider.setMinimum(0)
            self.horizontalSlider.setMaximum(self.player.duration())
            self.horizontalSlider.setValue(self.horizontalSlider.value() + 1000)
            self.horizontalSlider.setSliderPosition(0)
            # 绑定到analyse的进度条
            self.horizontalSlider_2.setMinimum(0)
            self.horizontalSlider_2.setMaximum(self.player.duration())
            self.horizontalSlider_2.setValue(self.horizontalSlider_2.value() + 1000)
            self.horizontalSlider_2.setSliderPosition(0)

            self.label_2.setText(os.path.basename(self.path))
            self.label_3.setText("File:"+os.path.basename(self.path))

    def Analyse_btn_start(self):##这里对应的是打开文件，并点击按钮
        try:
            if self.path:##要必须在打开文件之后才允许进行处理
                self.textBrowser.append("*********This file :"+str(os.path.basename(self.path))+"*********")
                self.progressBar.setValue(0)  ##每次允许处理时进度条归0
                self.Process.Audio_TimeDomain(self)  ##把实例传入进去
                self.Process.Audio_FrequencyDomain(self)
                self.Process.Audio_YuPuDomain(self)
                self.progressBar.setValue(100)
                self.textBrowser.append("Analyse Succeed!")
                self.textBrowser.append("---------  "+str(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))+"  ---------")

        except Exception as e:
            print(e)
            self.textBrowser.setText("There are some errors occuring when programme trying to open file")

    def palyMusic(self):
        try:
            if self.path:#这个path是当前的路径，如果path变了，那么就意味着更换了文件
                if not self.isPlay:##如果isPaly=0，那就说明播放器并没有打开，且此时按下了播放按钮，就开始播放
                    self.player.play()
                    self.isPlay=1##播放之后同时置为1，代表播放器目前正在播放
                else:
                    self.player.pause()
                    self.isPlay = 0  ##暂停之后同时置为0，代表播放器目前没有播放
        except Exception as e:
            print(e)
            self.textBrowser_2.setText("There are some errors occuring when playing audio")

    def playRefresh(self):
        if self.isPlay:
            #print(self.player.duration())
            self.horizontalSlider.setMinimum(0)
            self.horizontalSlider.setMaximum(self.player.duration())
            self.horizontalSlider.setValue(self.horizontalSlider.value() + 1000)

            self.horizontalSlider_2.setMinimum(0)
            self.horizontalSlider_2.setMaximum(self.player.duration())
            self.horizontalSlider_2.setValue(self.horizontalSlider_2.value() + 1000)
        elif self.isPlay_IIR:
            self.horizontalSlider_4.setMinimum(0)
            self.horizontalSlider_4.setMaximum(self.player_IIR.duration())
            self.horizontalSlider_4.setValue(self.horizontalSlider_4.value() + 1000)
        elif self.isPlay_FIR:
            self.horizontalSlider_5.setMinimum(0)
            self.horizontalSlider_5.setMaximum(self.player_FIR.duration())
            self.horizontalSlider_5.setValue(self.horizontalSlider_5.value() + 1000)
        #ORIGINAL AUDIO
        self.label_66.setText(time.strftime('%M:%S', time.localtime(self.player.position() / 1000)))
        self.label_67.setText(time.strftime('%M:%S', time.localtime(self.player.duration() / 1000)))

        self.label_68.setText(time.strftime('%M:%S', time.localtime(self.player.position() / 1000)))
        self.label_69.setText(time.strftime('%M:%S', time.localtime(self.player.duration() / 1000)))
        #
        # self.label_68.setText(time.strftime('%M:%S', time.localtime(self.player.position() / 1000)))
        # self.label_69.setText(time.strftime('%M:%S', time.localtime(self.player.duration() / 1000)))

        #IIR
        self.label_71.setText(time.strftime('%M:%S', time.localtime(self.player_IIR.position() / 1000)))
        self.label_70.setText(time.strftime('%M:%S', time.localtime(self.player_IIR.duration() / 1000)))
        #FIR
        self.label_73.setText(time.strftime('%M:%S', time.localtime(self.player_FIR.position() / 1000)))
        self.label_72.setText(time.strftime('%M:%S', time.localtime(self.player_FIR.duration() / 1000)))

    def changeVoice(self):
        #print(self.dial.value())
        self.player_IIR.setVolume(self.dial_4.value())
        self.player_FIR.setVolume(self.dial_5.value())

    def changeVoice0(self):
        self.player.setVolume(self.dial.value())
        self.dial.setValue(self.dial.value())
        # self.dial_2.setValue(self.dial.value())
        # self.dial_10.setValue(self.dial.value())
    def changeVoice01(self):
        self.player.setVolume(self.dial_2.value())
        self.dial_2.setValue(self.dial_2.value())

    def changeVoice1(self):
        self.player.setVolume(self.dial_4.value())
        self.dial.setValue(self.dial_4.value())
        # self.dial_10.setValue(self.dial_2.value())
    def changeVoice2(self):
        self.player.setVolume(self.dial_10.value())
        self.dial.setValue(self.dial_10.value())
        self.dial_2.setValue(self.dial_10.value())


    def desigenIIR(self):
        ###获取到输入参数：滤波器四个指标
        self.progressBar.setValue(0)
        try:
            self.An_wp=self.lineEdit.text()
            self.An_wst=self.lineEdit_2.text()
            self.Rp=self.lineEdit_3.text()
            self.As=self.lineEdit_4.text()

            self.fs=self.lineEdit_5.text()

            self.filterType=self.comboBox.currentText()
            self.iirType=self.comboBox_2.currentText()
            self.progressBar.setValue(10)
            self.Process.IIR_Designer(self)
        except Exception as e:
            print(e)
        self.progressBar.setValue(100)
    def applyIIR(self):
        self.progressBar.setValue(0)
        try:
            self.process_flag=0
            self.player_IIR.pause()
            self.player_IIR.setMedia(QMediaContent(QUrl(self.path)))  # 先把绑定改过去，不然文件占用
            self.progressBar.setValue(20)
            self.Process.apply_IIR(self)
            if self.process_flag:  # 如果处理好了
                try:
                    self.isPlay = 0
                    self.isPlay_IIR = 0
                    self.isPlay_FIR = 0
                    self.player.pause()  # 暂停另外的播放器
                    self.player_FIR.pause()
                    self.horizontalSlider_4.setMinimum(0)
                    self.horizontalSlider_4.setMaximum(self.player_IIR.duration())
                    self.horizontalSlider_4.setValue(self.horizontalSlider_4.value() + 1000)
                    self.horizontalSlider_4.setSliderPosition(0)
                    self.label_74.setText("Processed Audio: " + os.path.basename(self.path))

                    self.player_IIR.setMedia(QMediaContent(QUrl(self.saveDatepath_IIR)))  ##选中需要播放的音乐
                    #self.player_IIR.setMedia(QMediaContent(), buf)  # 从缓存里面读出来的

                #self.player_IIR.setMedia(QMediaContent(QUrl(self.path)))
                except Exception as e:
                    print(e)
            else:#没有处理好，也就是没有进行滤波操作
                self.textBrowser_2.setText("please choose a audio to design filter and apply before previewing")
        except Exception as e:
            print(e)
        self.progressBar.setValue(100)
    def playIIRaudio(self):
        try:
            self.isPlay=0#点按任意一个播放器的播放暂停按钮都会停止
            self.player.pause()
            if self.process_flag:  # 如果处理好了
                if not self.isPlay_IIR:
                    self.horizontalSlider_4.setValue(self.player_IIR.position())
                    self.player_IIR.play()
                    self.isPlay_IIR=1
                    self.textBrowser_2.append("play")
                else:##如果发现播放器正在播放
                    self.player_IIR.pause()
                    self.isPlay_IIR=0
                    self.textBrowser_2.append("pause")
            else:#没有处理好，也就是没有进行滤波操作
                self.textBrowser_2.setText("please choose a audio to design filter and apply before previewing")
        except Exception as e:
            print(e)

    def designFIR(self):
        try:
            self.progressBar.setValue(0)
            self.f1=float(self.lineEdit_10.text())

            self.f2=float(self.lineEdit_11.text())

            self.filter_length=float(self.lineEdit_13.text())

            # self.As= float(self.lineEdit_12.text())

            self.fs_FIR = float(self.lineEdit_23.text())

            self.filterType_FIR=self.comboBox_7.currentText()

            self.firType=self.comboBox_5.currentText()

            self.progressBar.setValue(20)
            self.Process.FIR_Designer(self)
            self.progressBar.setValue(100)
        except Exception as e:
            print(e)

    def applyFIR(self):
        self.progressBar.setValue(0)
        try:
            #处理之前全部关掉播放器对音乐的链接，不然会导致文件写不进去
            self.isPlay = 0
            self.isPlay_IIR = 0
            self.isPlay_FIR = 0
            self.player.pause()  # 暂停另外的播放器
            self.player_IIR.pause()
            self.player_FIR.pause()
            self.player_FIR.setMedia(QMediaContent(QUrl(self.path)))  # 先绑定到其他地方去，不然文件占用会导致写不进去文件
            self.process_flag=0

            self.progressBar.setValue(10)

            self.Process.apply_FIR(self)

            print(self.process_flag)
            if self.process_flag:  # 如果处理好了
                self.horizontalSlider_5.setMinimum(0)
                self.horizontalSlider_5.setMaximum(self.player_FIR.duration())
                self.horizontalSlider_5.setValue(self.horizontalSlider_5.value() + 1000)
                self.horizontalSlider_5.setSliderPosition(0)
                self.label_75.setText("Processed Audio: " + os.path.basename(self.path))

                self.player_FIR.setMedia(QMediaContent(QUrl(self.saveDatepath_FIR)))  ##选中需要播放的音乐
                #self.player_IIR.setMedia(QMediaContent(QUrl(self.path)))

            else:#没有处理好，也就是没有进行滤波操作
                self.textBrowser_3.setText("please choose a audio to design filter and apply before previewing")
        except Exception as e:
            print(e)
        self.progressBar.setValue(100)
    def playFIRaudio(self):
        try:
            self.isPlay=0#点按任意一个播放器的播放暂停按钮都会停止
            self.isPlay_IIR = 0
            self.player.stop()
            self.player_IIR.stop()
            print(self.process_flag)
            if self.process_flag:  # 如果处理好了
                if not self.isPlay_FIR:
                    self.horizontalSlider_5.setValue(self.player_FIR.position())#从停止位置继续播放
                    self.player_FIR.play()
                    self.isPlay_FIR=1
                    self.textBrowser_3.append("play")
                else:##如果发现播放器正在播放
                    self.player_FIR.pause()
                    self.isPlay_FIR=0
                    self.textBrowser_3.append("pause")
            else:#没有处理好，也就是没有进行滤波操作
                self.textBrowser_3.setText("please choose a audio to design filter and apply before previewing")
        except Exception as e:
            print(e)