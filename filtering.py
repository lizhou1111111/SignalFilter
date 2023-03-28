import wave

import scipy.signal as signal
import numpy as np
import pylab as plt


##信号处理函数类，里面的是一些信号处理的函数
from scipy.io import wavfile


class ProcessFunction(object):  ##这里负责写一些数字信号处理的方法
    def Audio_TimeDomain(self,feature):  ##时域
        f = wave.open(feature.path,"rb")
        params = f.getparams()
        nchannels, sampwidth, framerate, nframes = params[:4]
        # nchannels通道数
        # sampwidth量化位数
        # framerate采样频率
        # nframes采样点数
        str_data = f.readframes(nframes)
        f.close()
        # 将字符串转换为数组，得到一维的short类型的数组
        wave_data = np.fromstring(str_data, dtype=np.short)
        # 赋值的归一化
        wave_data = wave_data * 1.0 / (max(abs(wave_data)))
        # 整合左声道和右声道的数据
        wave_data = np.reshape(wave_data, [nframes, nchannels])
        # 最后通过采样点数和取样频率计算出每个取样的时间
        time = np.arange(0, nframes) * (1.0 / framerate)

        feature.textBrowser_2.append("AUDIO INFO:   Number of channel: " + str(nchannels))
        feature.textBrowser_2.append("AUDIO INFO:   Sampling Frequency: " + str(framerate)+" Hz")
        feature.textBrowser_2.append("AUDIO INFO:   Sampling number: " + str(nframes))
        feature.textBrowser_2.append("AUDIO INFO:   Sampling duration: " + str(nframes/framerate)+" seconds")

        ax = feature.fig5.add_subplot(111)
        ###进度条显示******
        feature.progressBar.setValue(10)
        #***************
        #调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        #ax.clear()
        ax.plot(time, wave_data[:, 0])
        ax.set_title('Normalized Magnitude')
        ax.set_xlabel('Time [sec]')

        feature.fig5.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas5.draw()  # TODO:这里开始绘制

        feature.progressBar.setValue(20)

    def Audio_FrequencyDomain(self,feature):
        #*********************STFT图像绘制*****************************
        sampling_freq, audio = wavfile.read(feature.path)
        T = 20  # 短时傅里叶变换的时长 单位 ms
        fs = sampling_freq  # 采样频率是这么多，那么做出来的频谱宽度就是这么多
        N = len(audio)  # 采样点的个数
        audio = audio * 1.0 / (max(abs(audio)))

        # 计算并绘制STFT的大小
        f, t, Zxx = signal.stft(audio, fs, nperseg=T * fs / 1000)

        feature.progressBar.setValue(30)

        ax = feature.fig7.add_subplot(111)
        feature.fig7.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        # ax=plt.figure()
        ax.pcolormesh(t, f, np.abs(Zxx), vmin=0, vmax=0.1)
        ax.set_title('STFT Magnitude')
        #feature.fig6.colorbar(ax=ax)
        #feature.fig6.colorbar(feature.fig6)
        ####还存在的问题是colorbar显示不了的问题####
        ax.set_xlabel('Time [sec]')
        ax.set_ylabel('Frequency [Hz]')
        feature.canvas7.draw()  # TODO:这里开始绘制
        feature.progressBar.setValue(40)
        #**************************************************************

        # *******************FFT图像绘制*********************************
        fft_signal = np.fft.fft(audio)
        fft_signal = abs(fft_signal)
        # 建立频率轴

        # 建立频率轴

        fft_signal = np.fft.fftshift(fft_signal)

        fft_signal = fft_signal[int(fft_signal.shape[0] / 2):]

        freqInteral = (sampling_freq / len(fft_signal))  ###频率轴的间隔

        Freq = np.arange(0, sampling_freq / 2, sampling_freq / (2*len(fft_signal)))

        feature.progressBar.setValue(50)
        #
        highFreq=  (np.argmax(fft_signal[int(len(fft_signal) / 2):len(fft_signal)]) )*freqInteral
        feature.textBrowser_2.append("FFT INFO:   Highest frequency: "+str(highFreq))


        ax = feature.fig6.add_subplot(111)
        # 调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        ax.plot(Freq, fft_signal, color='red')
        ax.set_title('FFT Figure')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Am')
        feature.fig6.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas6.draw()  # TODO:这里开始绘制

        feature.progressBar.setValue(60)

    def Audio_YuPuDomain(self, feature):
        #*******************语谱图绘制***********************************
        f = wave.open(feature.path, "rb")
        params = f.getparams()
        nchannels, sampwidth, framerate, nframes = params[:4]
        strData = f.readframes(nframes)  # 读取音频，字符串格式
        waveData = np.fromstring(strData, dtype=np.int16)  # 将字符串转化为int
        waveData = waveData * 1.0 / (max(abs(waveData)))  # wave幅值归一化
        waveData = np.reshape(waveData, [nframes, nchannels]).T
        f.close()
        feature.progressBar.setValue(70)

        ax = feature.fig8.add_subplot(111)
        feature.fig8.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        # ax=plt.figure()
        ax.specgram(waveData[0],Fs = framerate, scale_by_freq = True, sides = 'default')
        ax.set_title('Spectrogram')
        # feature.fig6.colorbar(ax=ax)
        # feature.fig6.colorbar(feature.fig6)
        ####还存在的问题是colorbar显示不了的问题####
        ax.set_xlabel('Time [sec]')
        ax.set_ylabel('Frequency [Hz]')
        feature.canvas8.draw()  # TODO:这里开始绘制
        feature.progressBar.setValue(80)

    def IIR_Designer(self,feature):
        if str(feature.iirType)=='Butterworth':##巴特沃斯 双线性变换法间接设计模拟滤波器
            fs = float(feature.fs)
            feature.textBrowser_3.append(str(feature.filterType))
            feature.textBrowser_3.append(str(feature.iirType))
            if str(feature.filterType)=="Bandpass" or str(feature.filterType)=="bandstop":
                ##如果是带通带阻需要输入四组频率数据
                ##########频率预畸##################
                wp=str(feature.An_wp).split()#切分开之后再转换为array
                wp0=float(wp[0]) * (2 * np.pi / fs)
                wp1=float(wp[1]) * (2 * np.pi / fs)

                #************************************
                wp[0]=(2 * fs) * np.tan(wp0 / 2)#双线性变换
                wp[1]=(2 * fs) * np.tan(wp1 / 2)

                omiga_p=[float(wp[0]),float(wp[1])]
                wst=str(feature.An_wst).split()#切分开之后再转换为array
                wst0=float(wst[0]) * (2 * np.pi / fs)
                wst1=float(wst[1]) * (2 * np.pi / fs)
                #************************************

                wst[0]=(2 * fs) * np.tan(wst0 / 2)#双线性变换
                wst[1]=(2 * fs) * np.tan(wst1 / 2)

                omiga_st = [wst[0], wst[1]]
            else:
                wp = float(feature.An_wp) * (2 * np.pi / fs)
                wst = float(feature.An_wst) * (2 * np.pi / fs)
                ##########频率预畸##################
                omiga_p = (2 * fs) * np.tan(wp / 2)
                omiga_st = (2 * fs) * np.tan(wst / 2)
            feature.Rp=float(feature.Rp)
            feature.As=float(feature.As)
            N, Wn = signal.buttord(omiga_p, omiga_st, feature.Rp, feature.As, True)
            feature.filts = signal.lti(*signal.butter(N, Wn, btype=str(feature.filterType),
                                              analog=True))
            feature.filtz = signal.lti(*signal.bilinear(feature.filts.num, feature.filts.den, fs))

            feature.z, feature.p = signal.bilinear(feature.filts.num, feature.filts.den, fs)

            wz, hz = signal.freqz(feature.filtz.num, feature.filtz.den)

            ax = feature.fig1.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            #ax.clear()
            ax.semilogx(wz * fs / (2 * np.pi), 20 * np.log10(np.abs(hz).clip(1e-15)),
                     label=r'$|H_z(e^{j \omega})|$')
            ax.set_xlabel('Hz')
            ax.set_ylabel('dB')
            ax.set_title('Butterworth')
            feature.fig1.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas1.draw()  # TODO:这里开始绘制

            ###绘制零极点图########
            ax = feature.fig3.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            z1, p1, k1 = signal.tf2zpk(feature.z, feature.p)  # zero, pole and gain
            c = np.vstack((feature.p, feature.z))
            Max = (abs(c)).max()  # find the largest value
            a = feature.p / Max  # normalization
            b = feature.z / Max
            Ra = (a * (2 ** ((N - 1) - 1))).astype(int)  # quantizan and truncate
            Rb = (b * (2 ** ((N - 1) - 1))).astype(int)
            z2, p2, k2 = signal.tf2zpk(Rb, Ra)
            ##参数方程画圆
            theta = np.arange(0, 2 * np.pi, 0.01)
            x = np.cos(theta)
            y = np.sin(theta)
            ax.plot(x,y,color='black')
            for i in p1:
                ax.plot(np.real(i), np.imag(i), 'bx')  # pole before quantization
            for i in z1:
                ax.plot(np.real(i), np.imag(i), 'bo')  # zero before quantization
            for i in p2:
                ax.plot(np.real(i), np.imag(i), 'rx')  # pole after quantization
            for i in z2:
                ax.plot(np.real(i), np.imag(i), 'ro')  # zero after quantization
            ax.set_xlim(-1.8, 1.8)
            ax.set_ylim(-1.2, 1.2)
            ax.grid()
            ax.set_title("%d bit quantization" % N)
            feature.fig3.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas3.draw()  # TODO:这里开始绘制

        if str(feature.iirType) == 'Chebyshev I':  ##切比雪夫一型

            fs = float(feature.fs)
            feature.textBrowser_3.append(str(feature.filterType))
            feature.textBrowser_3.append(str(feature.iirType))
            if str(feature.filterType)=="Bandpass" or str(feature.filterType)=="bandstop":
                ##如果是带通带阻需要输入四组频率数据
                ##########频率预畸##################
                wp=str(feature.An_wp).split()#切分开之后再转换为array
                wp0=float(wp[0]) * (2 * np.pi / fs)
                wp1=float(wp[1]) * (2 * np.pi / fs)

                #************************************
                wp[0]=(2 * fs) * np.tan(wp0 / 2)#双线性变换
                wp[1]=(2 * fs) * np.tan(wp1 / 2)

                omiga_p=[float(wp[0]),float(wp[1])]
                wst=str(feature.An_wst).split()#切分开之后再转换为array
                wst0=float(wst[0]) * (2 * np.pi / fs)
                wst1=float(wst[1]) * (2 * np.pi / fs)
                #************************************

                wst[0]=(2 * fs) * np.tan(wst0 / 2)#双线性变换
                wst[1]=(2 * fs) * np.tan(wst1 / 2)

                omiga_st = [wst[0], wst[1]]
            else:
                wp = float(feature.An_wp) * (2 * np.pi / fs)
                wst = float(feature.An_wst) * (2 * np.pi / fs)
                ##########频率预畸##################
                omiga_p = (2 * fs) * np.tan(wp / 2)
                omiga_st = (2 * fs) * np.tan(wst / 2)

            if len(str(feature.Rp).split())>1: #纹波参数
                Rpinput=str(feature.Rp).split()
                feature.Rp = float(Rpinput[0])
                feature.As = float(feature.As)
                rp_in=float(Rpinput[1])
            else:
                feature.Rp=float(feature.Rp)
                feature.As=float(feature.As)
                rp_in = 0.1*feature.Rp

            #N, Wn = signal.cheb1ord(wp, wst, feature.Rp, feature.As, True)
            N, Wn = signal.cheb1ord(omiga_p, omiga_st, feature.Rp, feature.As, True)
            feature.filts = signal.lti(*signal.cheby1(N, rp_in,Wn, btype=str(feature.filterType),
                                              analog=True))##切比雪夫是还有一个纹波参数
            feature.filtz = signal.lti(*signal.bilinear(feature.filts.num, feature.filts.den, fs))

            feature.z,feature.p=signal.bilinear(feature.filts.num, feature.filts.den, fs)

            wz, hz = signal.freqz(feature.filtz.num, feature.filtz.den)

            ax = feature.fig1.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            #ax.clear()
            ax.semilogx(wz * fs / (2 * np.pi), 20 * np.log10(np.abs(hz).clip(1e-15)),
                     label=r'$|H_z(e^{j \omega})|$')
            ax.set_xlabel('Hz')
            ax.set_ylabel('dB')
            ax.set_title('Chebyshev I')
            feature.fig1.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas1.draw()  # TODO:这里开始绘制

            ###绘制零极点图########
            ax = feature.fig3.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            z1, p1, k1 = signal.tf2zpk(feature.z, feature.p)  # zero, pole and gain
            c = np.vstack((feature.p, feature.z))
            Max = (abs(c)).max()  # find the largest value
            a = feature.p / Max  # normalization
            b = feature.z / Max
            Ra = (a * (2 ** ((N - 1) - 1))).astype(int)  # quantizan and truncate
            Rb = (b * (2 ** ((N - 1) - 1))).astype(int)
            z2, p2, k2 = signal.tf2zpk(Rb, Ra)
            ##参数方程画圆
            theta = np.arange(0, 2 * np.pi, 0.01)
            x = np.cos(theta)
            y = np.sin(theta)
            ax.plot(x,y,color='black')
            for i in p1:
                ax.plot(np.real(i), np.imag(i), 'bx')  # pole before quantization
            for i in z1:
                ax.plot(np.real(i), np.imag(i), 'bo')  # zero before quantization
            for i in p2:
                ax.plot(np.real(i), np.imag(i), 'rx')  # pole after quantization
            for i in z2:
                ax.plot(np.real(i), np.imag(i), 'ro')  # zero after quantization
            ax.set_xlim(-1.8, 1.8)
            ax.set_ylim(-1.2, 1.2)
            ax.grid()
            ax.set_title("%d bit quantization" % N)
            feature.fig3.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas3.draw()  # TODO:这里开始绘制

        if str(feature.iirType) == 'Chebyshev II':  ##切比雪夫二型

            fs = float(feature.fs)
            feature.textBrowser_3.append(str(feature.filterType))
            feature.textBrowser_3.append(str(feature.iirType))
            if str(feature.filterType)=="Bandpass" or str(feature.filterType)=="bandstop":
                ##如果是带通带阻需要输入四组频率数据
                ##########频率预畸##################
                wp=str(feature.An_wp).split()#切分开之后再转换为array
                wp0=float(wp[0]) * (2 * np.pi / fs)
                wp1=float(wp[1]) * (2 * np.pi / fs)

                #************************************
                wp[0]=(2 * fs) * np.tan(wp0 / 2)#双线性变换
                wp[1]=(2 * fs) * np.tan(wp1 / 2)

                omiga_p=[float(wp[0]),float(wp[1])]
                wst=str(feature.An_wst).split()#切分开之后再转换为array
                wst0=float(wst[0]) * (2 * np.pi / fs)
                wst1=float(wst[1]) * (2 * np.pi / fs)
                #************************************

                wst[0]=(2 * fs) * np.tan(wst0 / 2)#双线性变换
                wst[1]=(2 * fs) * np.tan(wst1 / 2)

                omiga_st = [wst[0], wst[1]]
            else:
                wp = float(feature.An_wp) * (2 * np.pi / fs)
                wst = float(feature.An_wst) * (2 * np.pi / fs)
                ##########频率预畸##################
                omiga_p = (2 * fs) * np.tan(wp / 2)
                omiga_st = (2 * fs) * np.tan(wst / 2)
            if len(str(feature.As).split())>1: #纹波参数
                Asinput=str(feature.As).split()
                feature.As = float(Asinput[0])
                feature.Rp = float(feature.Rp)
                rs_in=float(Asinput[1])
            else:
                feature.Rp=float(feature.Rp)
                feature.As=float(feature.As)
                rs_in = 0.1*feature.As
            N, Wn = signal.cheb2ord(omiga_p, omiga_st, feature.Rp, feature.As, True)
            feature.filts = signal.lti(*signal.cheby2(N, rs_in,Wn, btype=str(feature.filterType),
                                              analog=True))##切比雪夫是还有一个纹波参数
            feature.filtz = signal.lti(*signal.bilinear(feature.filts.num, feature.filts.den, fs))

            feature.z,feature.p=signal.bilinear(feature.filts.num, feature.filts.den, fs)

            wz, hz = signal.freqz(feature.filtz.num, feature.filtz.den)

            ax = feature.fig1.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            #ax.clear()
            ax.semilogx(wz * fs / (2 * np.pi), 20 * np.log10(np.abs(hz).clip(1e-15)),
                     label=r'$|H_z(e^{j \omega})|$')
            ax.set_xlabel('Hz')
            ax.set_ylabel('dB')
            ax.set_title('Chebyshev I')
            feature.fig1.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas1.draw()  # TODO:这里开始绘制

            ###绘制零极点图########
            ax = feature.fig3.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            z1, p1, k1 = signal.tf2zpk(feature.z, feature.p)  # zero, pole and gain
            c = np.vstack((feature.p, feature.z))
            Max = (abs(c)).max()  # find the largest value
            a = feature.p / Max  # normalization
            b = feature.z / Max
            Ra = (a * (2 ** ((N - 1) - 1))).astype(int)  # quantizan and truncate
            Rb = (b * (2 ** ((N - 1) - 1))).astype(int)
            z2, p2, k2 = signal.tf2zpk(Rb, Ra)
            ##参数方程画圆
            theta = np.arange(0, 2 * np.pi, 0.01)
            x = np.cos(theta)
            y = np.sin(theta)
            ax.plot(x,y,color='black')
            for i in p1:
                ax.plot(np.real(i), np.imag(i), 'bx')  # pole before quantization
            for i in z1:
                ax.plot(np.real(i), np.imag(i), 'bo')  # zero before quantization
            for i in p2:
                ax.plot(np.real(i), np.imag(i), 'rx')  # pole after quantization
            for i in z2:
                ax.plot(np.real(i), np.imag(i), 'ro')  # zero after quantization
            ax.set_xlim(-1.8, 1.8)
            ax.set_ylim(-1.2, 1.2)
            ax.grid()
            ax.set_title("%d bit quantization" % N)
            feature.fig3.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas3.draw()  # TODO:这里开始绘制

        if str(feature.iirType) == 'Cauer/elliptic':
            fs = float(feature.fs)
            feature.textBrowser_3.append(str(feature.filterType))
            feature.textBrowser_3.append(str(feature.iirType))
            if str(feature.filterType)=="Bandpass" or str(feature.filterType)=="bandstop":
                ##如果是带通带阻需要输入四组频率数据
                ##########频率预畸##################
                wp=str(feature.An_wp).split()#切分开之后再转换为array
                wp0=float(wp[0]) * (2 * np.pi / fs)
                wp1=float(wp[1]) * (2 * np.pi / fs)

                #************************************
                wp[0]=(2 * fs) * np.tan(wp0 / 2)#双线性变换
                wp[1]=(2 * fs) * np.tan(wp1 / 2)

                omiga_p=[float(wp[0]),float(wp[1])]
                wst=str(feature.An_wst).split()#切分开之后再转换为array
                wst0=float(wst[0]) * (2 * np.pi / fs)
                wst1=float(wst[1]) * (2 * np.pi / fs)
                #************************************

                wst[0]=(2 * fs) * np.tan(wst0 / 2)#双线性变换
                wst[1]=(2 * fs) * np.tan(wst1 / 2)

                omiga_st = [wst[0], wst[1]]
            else:
                wp = float(feature.An_wp) * (2 * np.pi / fs)
                wst = float(feature.An_wst) * (2 * np.pi / fs)
                ##########频率预畸##################
                omiga_p = (2 * fs) * np.tan(wp / 2)
                omiga_st = (2 * fs) * np.tan(wst / 2)

            Asinput = str(feature.As).split()
            Rpinput = str(feature.Rp).split()
            feature.As = float(Asinput[0])
            feature.Rp = float(Rpinput[0])
            rs_in = float(Asinput[1])
            rp_in = float(Rpinput[1])

            feature.Rp=float(feature.Rp)
            feature.As=float(feature.As)
            N, Wn = signal.ellipord(omiga_p, omiga_st, feature.Rp, feature.As, True)
            feature.filts = signal.lti(*signal.ellip(N, rp_in,rs_in,Wn, btype=str(feature.filterType),
                                              analog=True))##切比雪夫是还有一个纹波参数
            feature.filtz = signal.lti(*signal.bilinear(feature.filts.num, feature.filts.den, fs))

            feature.z,feature.p=signal.bilinear(feature.filts.num, feature.filts.den, fs)

            wz, hz = signal.freqz(feature.filtz.num, feature.filtz.den)

            ax = feature.fig1.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            #ax.clear()
            ax.semilogx(wz * fs / (2 * np.pi), 20 * np.log10(np.abs(hz).clip(1e-15)),
                     label=r'$|H_z(e^{j \omega})|$')
            ax.set_xlabel('Hz')
            ax.set_ylabel('dB')
            ax.set_title('Cauer/elliptic')
            feature.fig1.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas1.draw()  # TODO:这里开始绘制

            ###绘制零极点图########
            ax = feature.fig3.add_subplot(111)
            ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
            z1, p1, k1 = signal.tf2zpk(feature.z, feature.p)  # zero, pole and gain
            c = np.vstack((feature.p, feature.z))
            Max = (abs(c)).max()  # find the largest value
            a = feature.p / Max  # normalization
            b = feature.z / Max
            Ra = (a * (2 ** ((N - 1) - 1))).astype(int)  # quantizan and truncate
            Rb = (b * (2 ** ((N - 1) - 1))).astype(int)
            z2, p2, k2 = signal.tf2zpk(Rb, Ra)
            ##参数方程画圆
            theta = np.arange(0, 2 * np.pi, 0.01)
            x = np.cos(theta)
            y = np.sin(theta)
            ax.plot(x,y,color='black')
            for i in p1:
                ax.plot(np.real(i), np.imag(i), 'bx')  # pole before quantization
            for i in z1:
                ax.plot(np.real(i), np.imag(i), 'bo')  # zero before quantization
            for i in p2:
                ax.plot(np.real(i), np.imag(i), 'rx')  # pole after quantization
            for i in z2:
                ax.plot(np.real(i), np.imag(i), 'ro')  # zero after quantization
            ax.set_xlim(-1.8, 1.8)
            ax.set_ylim(-1.2, 1.2)
            ax.grid()
            ax.set_title("%d bit quantization" % N)
            feature.fig3.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
            feature.canvas3.draw()  # TODO:这里开始绘制

        feature.textBrowser.setText("PARAMETER OF THIS FILTER")
        feature.textBrowser.append("*********" )
        feature.textBrowser.append("FILTER TPYE=" +str(feature.filterType))
        feature.textBrowser.append("IIR TPYE=" + str(feature.iirType))
        feature.textBrowser.append("ORDER=" + str(N))
        feature.textBrowser.append("b="+str(feature.z))
        feature.textBrowser.append("a="+str(feature.p))
        feature.textBrowser.append()
    def apply_IIR(self,feature):
        f = wave.open(feature.path,"rb")
        params = f.getparams()
        nchannels, sampwidth, framerate, nframes = params[:4]
        # nchannels通道数
        # sampwidth量化位数
        # framerate采样频率
        # nframes采样点数
        str_data = f.readframes(nframes)
        f.close()
        # 将字符串转换为数组，得到一维的short类型的数组
        wave_data = np.fromstring(str_data, dtype=np.short)
        # 赋值的归一化
        maximum=max(abs(wave_data))
        wave_data = wave_data * 1.0 / (maximum)
        # 整合左声道和右声道的数据
        wave_data = np.reshape(wave_data, [nframes, nchannels])

        # 最后通过采样点数和取样频率计算出每个取样的时间
        time = np.arange(0, nframes) * (1.0 / framerate)
        #print(time)
        t = np.linspace(0, nframes/ framerate, nframes, endpoint=False)
        print("p maxium")
        print(feature.p)
        print(feature.z)
        feature.yout=signal.filtfilt(feature.z,feature.p,wave_data[:, 0],method='gust')
        print(max(feature.yout))
        ax = feature.fig2.add_subplot(111)
        #调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        #ax.clear()
        ax.plot(t, feature.yout)
        ax.set_title('Passed Filter')
        ax.set_xlabel('Time [sec]')

        feature.fig2.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas2.draw()  # TODO:这里开始绘制

        ##绘制出时域的图像之后，再到频率分析
        #FFT变换#
        fft_signal = np.fft.fft(feature.yout)
        fft_signal = np.fft.fftshift(abs(fft_signal))
        fft_signal = fft_signal[int(fft_signal.shape[0]/2):]
        # 建立频率轴
        Freq = np.arange(0, framerate / 2, framerate / (2*len(fft_signal)))

        ####绘图######
        ax = feature.fig4.add_subplot(111)
        # 调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        ax.plot(Freq, fft_signal, color='red')
        ax.set_title('FFT Figure')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Am')
        feature.fig4.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas4.draw()  # TODO:这里开始绘制


        #feature.precessed_Audio=feature.filtz.output(wave_data,time,X0=None)#求系统的零状态响应
        # feature.precessed_Audio =feature.precessed_Audio.tostring()
        # feature.process_flag=1#标志位为1，代表处理好了，否则的话就代表没有
        ##写音频文件###

        feature.yout = feature.yout * maximum  # 去归一化
        feature.yout = feature.yout.astype(np.short)
        f = wave.open(feature.saveDatepath_IIR, "wb")  ##
        f.setnchannels(nchannels)
        f.setsampwidth(sampwidth)
        f.setframerate(framerate)
        f.setnframes(nframes)
        f.writeframes(feature.yout.tostring())
        f.close()

        feature.process_flag=1#代表本次处理完毕

    def FIR_Designer(self,feature):
        kaiser_para=0.85
        if feature.filterType_FIR == 'Lowpass':
            numtaps=int(feature.filter_length)
            fcut=feature.f2*2/feature.fs_FIR
            if str(feature.firType)=='kaiser':
                width=kaiser_para
            else:
                width=None
            feature.FIR_b = signal.firwin(numtaps, fcut,width=width, window=str(feature.firType))  #
        if feature.filterType_FIR == 'Highpass':
            numtaps=int(feature.filter_length)
            fcut=feature.f2*2/feature.fs_FIR
            if str(feature.firType)=='kaiser':
                width=kaiser_para
            else:
                width=None
            feature.FIR_b = signal.firwin(numtaps, fcut,width=width,window=str(feature.firType),pass_zero=False)  #
        if feature.filterType_FIR == 'Bandpass':
            numtaps = int(feature.filter_length)
            fcut = [feature.f1*2/feature.fs_FIR,feature.f2*2/feature.fs_FIR]
            if str(feature.firType)=='kaiser':
                width=kaiser_para
            else:
                width=None
            feature.FIR_b = signal.firwin(numtaps, fcut,width=width,window=str(feature.firType), pass_zero=False)  #
        if feature.filterType_FIR == 'Band-stop pass':
            numtaps = int(feature.filter_length)
            fcut = [feature.f1*2/feature.fs_FIR,feature.f2*2/feature.fs_FIR]
            if str(feature.firType)=='kaiser':
                width=kaiser_para
            else:
                width=None
            feature.FIR_b = signal.firwin(numtaps, fcut,width=width,window=str(feature.firType))  #


        feature.textBrowser_3.append("FilterType_FIR:"+feature.filterType_FIR)
        feature.textBrowser_3.append("FirType"+str(feature.firType))
        feature.textBrowser_3.append("INFO:     *****Succeed!*****    ")
        # 绘制频率响应：
        wz, hz = signal.freqz(feature.FIR_b)

        ax = feature.fig25.add_subplot(111)
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        # ax.clear()
        ax.semilogx(wz * feature.fs_FIR / (2 * np.pi), 20 * np.log10(np.abs(hz).clip(1e-15)),
                    label=r'$|H_z(e^{j \omega})|$')
        ax.set_xlabel('Hz')
        ax.set_ylabel('dB')
        ax.set_title(str(feature.firType))
        feature.fig25.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas25.draw()  # TODO:这里开始绘制


        #####绘制零极点图###############
        ##只有零点没有极点###
        ax = feature.fig26.add_subplot(111)
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        fir_a=np.zeros(numtaps)
        fir_a[numtaps-1]=1
        z1, p1, k1 = signal.tf2zpk(feature.FIR_b, fir_a)  # zero, pole and gain
        c = np.vstack((fir_a, feature.FIR_b))
        Max = (abs(c)).max()  # find the largest value
        a = fir_a / Max  # normalization
        b = feature.FIR_b / Max
        Ra = (a * (2 ** ((numtaps - 1) - 1))).astype(int)  # quantizan and truncate
        Rb = (b * (2 ** ((numtaps - 1) - 1))).astype(int)
        z2, p2, k2 = signal.tf2zpk(Rb, Ra)
        ##参数方程画圆
        theta = np.arange(0, 2 * np.pi, 0.01)
        x = np.cos(theta)
        y = np.sin(theta)
        ax.plot(x, y, color='black')
        for i in p1:
            ax.plot(np.real(i), np.imag(i), 'bx')  # pole before quantization
        for i in z1:
            ax.plot(np.real(i), np.imag(i), 'bo')  # zero before quantization
        for i in p2:
            ax.plot(np.real(i), np.imag(i), 'rx')  # pole after quantization
        for i in z2:
            ax.plot(np.real(i), np.imag(i), 'ro')  # zero after quantization
        ax.set_xlim(-1.8, 1.8)
        ax.set_ylim(-1.2, 1.2)
        ax.grid()
        ax.set_title("%d bit quantization" % numtaps)
        feature.fig26.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas26.draw()  # TODO:这里开始绘制
        feature.textBrowser_11.setText("PARAMETER OF THIS FILTER")
        feature.textBrowser_11.append("*********")
        feature.textBrowser_11.append("FILTER TPYE=" + str(feature.filterType_FIR))
        feature.textBrowser_11.append("FIR TPYE=" + str(feature.firType))
        feature.textBrowser_11.append("ORDER(length of filter)=" + str(numtaps))
        feature.textBrowser_11.append("b=" + str(feature.FIR_b))
        feature.textBrowser_11.append("a=" + str(fir_a))
    def apply_FIR(self,feature):
        f = wave.open(feature.path,"rb")
        params = f.getparams()
        nchannels, sampwidth, framerate, nframes = params[:4]
        # nchannels通道数
        # sampwidth量化位数
        # framerate采样频率
        # nframes采样点数
        str_data = f.readframes(nframes)
        f.close()
        # 将字符串转换为数组，得到一维的short类型的数组
        wave_data = np.fromstring(str_data, dtype=np.short)
        # 赋值的归一化
        maximum=max(abs(wave_data))
        wave_data = wave_data * 1.0 / (maximum)
        # 整合左声道和右声道的数据
        wave_data = np.reshape(wave_data, [nframes, nchannels])

        # 最后通过采样点数和取样频率计算出每个取样的时间
        time = np.arange(0, nframes) * (1.0 / framerate)
        #print(time)
        t = np.linspace(0, nframes/ framerate, nframes, endpoint=False)

        feature.yout=signal.filtfilt(feature.FIR_b,1,wave_data[:, 0])

        ax = feature.fig27.add_subplot(111)
        #调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        #ax.clear()
        ax.plot(t, feature.yout)
        ax.set_title('Passed Filter')
        ax.set_xlabel('Time [sec]')

        feature.fig27.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas27.draw()  # TODO:这里开始绘制

        ##绘制出时域的图像之后，再到频率分析
        #FFT变换#
        fft_signal = np.fft.fft(feature.yout)
        fft_signal = np.fft.fftshift(abs(fft_signal))[int(fft_signal.shape[0]/2):]
        # 建立频率轴
        Freq = np.arange(0, framerate / 2, framerate / (2*len(fft_signal)))

        ####绘图######
        ax = feature.fig28.add_subplot(111)
        # 调整图像大小
        ax.cla()  # TODO:删除原图，让画布上只有新的一次的图
        ax.plot(Freq, fft_signal, color='red')
        ax.set_title('FFT Figure')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Am')
        feature.fig28.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
        feature.canvas28.draw()  # TODO:这里开始绘制


        #feature.precessed_Audio=feature.filtz.output(wave_data,time,X0=None)#求系统的零状态响应
        # feature.precessed_Audio =feature.precessed_Audio.tostring()
        # feature.process_flag=1#标志位为1，代表处理好了，否则的话就代表没有
        ##写音频文件###
        feature.yout = feature.yout * maximum  # 去归一化
        feature.yout = feature.yout.astype(np.short)
        f = wave.open(feature.saveDatepath_FIR, "wb")  ##
        f.setnchannels(nchannels)
        f.setsampwidth(sampwidth)
        f.setframerate(framerate)
        f.setnframes(nframes)
        f.writeframes(feature.yout.tostring())
        f.close()
        feature.process_flag=1#代表本次处理完毕


# def filtering():
#     # 思路：滤波器顾名思义是滤除特定频率滤波的工具，把时域信号转化为频域信号，再行滤除。
#
#
#
#
#
#     # # 带通IIR滤波器
#     # # freqz计算的滤波器频谱
#     # b, a = signal.iirdesign([0.2, 0.5], [0.1, 0.6], 2, 40)
#     # w, h = signal.freqz(b, a)
#     # power = 20 * np.log10(np.clip(np.abs(h), 1e-8, 1e100))
#     # plt.plot(w / np.pi * 4000, power)
#     # plt.show()
#     #
#     #
#     # # 频率扫描得到的滤波器频谱
#     # t = np.arange(0, 2, 1 / 8000.0)
#     # sweep = signal.chirp(t, f0=0, t1=2, f1=4000.0)
#     # out = signal.lfilter(b, a, sweep)
#     # out = 20 * np.log10(np.abs(out))
#     # index = np.where(np.logical_and(out[1:-1] > out[:-2], out[1:-1] > out[2:]))[0] + 1
#     # plt.plot(t[index] / 2.0 * 4000, out[index])
#     # plt.show()
#
#     #  巴特沃斯滤波器
#     b, a = signal.butter(4, 100, 'low', analog=True)  # 4阶低通临界频率为100Hz
#     w, h = signal.freqs(b, a)  # h为频率响应,w为频率
#     plt.figure(1)
#     plt.semilogx(w, 20 * np.log10(abs(h)))  # 用于绘制折线图，两个函数的 x 轴、y 轴分别是指数型的，并且转化为分贝
#     plt.title('Butterworth filter frequency response')
#     plt.xlabel('Frequency [radians / second]')
#     plt.ylabel('Amplitude [dB]')
#     plt.margins(0, 0.1)  # 去除画布四周白边
#     plt.grid(which='both',
#              axis='both')  # 生成网格，matplotlin.pyplot.grid(b, which, axis, color, linestyle, linewidth， **kwargs)， which : 取值为'major', 'minor'， 'both'
#     plt.axvline(100, color='green')  # 绘制竖线
#     plt.show()
#     t = np.linspace(0, 1, 1000, False)  # 1秒，1000赫兹刻度
#     sig = np.sin(2 * np.pi * 10 * t) + np.sin(2 * np.pi * 20 * t)  # 合成信号
#     fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)  # 2行1列的图
#     ax1.plot(t, sig)
#     ax1.set_title('10 Hz and 20 Hz sinusoids')
#     ax1.axis([0, 1, -2, 2])  # 坐标范围

