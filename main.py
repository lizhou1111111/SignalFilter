from PyQt5.QtWidgets import QApplication
import sys
from windows import MySignal


def print_hi(name):
    print(f'Hi, {name}')  #


if __name__ == '__main__':

    # filtering.filtering()

    print_hi('PyCharm')
    # 获取运行python文件的时候命令行参数，且以list形式存储参数
    app = QApplication(sys.argv)
    mySignal = MySignal()
    mySignal.show()
    sys.exit(app.exec_())

