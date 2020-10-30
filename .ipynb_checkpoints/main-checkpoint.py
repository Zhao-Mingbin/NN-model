#############################################模块引进#################################################
from nn import Ui_MainWindow
from nn_therm import *
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QFormLayout, QLineEdit,QMessageBox
from PyQt5.QtGui import QPixmap
from scipy.stats import norm
from PyQt5 import QtCore, QtGui, QtWidgets
#############################################功能函数#################################################
def radio_button(btn,btn_line):
    """[用于实现将被选中的流体参数的输入框变为不可输入状态的功能]

    Args:
        btn ([QRadioButton]): [选择参数框中的流体参数名称]
        btn_line ([QLineEdit]): [流体参数编辑框]
    """    
    if btn.isChecked()==True:
        btn_line.setEnabled(False)

    if btn.isChecked() == False:
        btn_line.setEnabled(True)

class UnitConv:
    """[用于实现输入输出结果的单位换算]
    """
    length_output = {}
    pressure_output = {}
    temperature_output = {}
    rho_output = {}
    def length_Unit_conversion(self,input_data,input_line):

        """[用于实现长度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [长度单位选择框]
            input_line ([QLineEdit]): [长度数值编辑框]
        """   

        try:
            if input_data.currentText() == "分米（dm）":
                UnitConv.length_output["分米（dm）"]=float(input_line.text())
                UnitConv.length_output["厘米（cm）"] = float(input_line.text())*10
                UnitConv.length_output["毫米（mm）"] = float(input_line.text())*100
                UnitConv.length_output["米（m）"] = float(input_line.text())/10

            elif input_data.currentText()=="厘米（cm）":
                UnitConv.length_output["分米（dm）"]=float(input_line.text())/10
                UnitConv.length_output["厘米（cm）"] = input_line.text()
                UnitConv.length_output["毫米（mm）"] = float(input_line.text())*10
                UnitConv.length_output["米（m）"] = float(input_line.text())/100

            elif input_data.currentText()=="毫米（mm）":
                UnitConv.length_output["分米（dm）"]=float(input_line.text())/100
                UnitConv.length_output["厘米（cm）"] = float(input_line.text())/10
                UnitConv.length_output["毫米（mm）"] = float(input_line.text())
                UnitConv.length_output["米（m）"] = float(input_line.text())/1000

            elif input_data.currentText()=="米（m）":
                UnitConv.length_output["分米（dm）"]=float(input_line.text())/10
                UnitConv.length_output["厘米（cm）"] =float (input_line.text())/100
                UnitConv.length_output["毫米（mm）"] = float(input_line.text())/1000
                UnitConv.length_output["米（m）"] = float(input_line.text())
        except:
            pass
    def length_refresh(self,input_data,input_line):

        """[用于实现长度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [长度单位选择框]
            input_line ([QLineEdit]): [长度数值编辑框]
        """        

        try:
            if input_data.currentText() == "分米（dm）":
                input_line.setText(str(UnitConv.length_output["分米（dm）"]))

            elif input_data.currentText()=="厘米（cm）":
                input_line.setText(str(UnitConv.length_output["厘米（cm）"]))

            elif input_data.currentText() == "毫米（mm）":
                input_line.setText(str(UnitConv.length_output["毫米（mm）"]))

            elif input_data.currentText() == "米（m）":
                input_line.setText(str(UnitConv.length_output["米（m）"]))

        except:
            pass

    def _pressure_Unit_conversion(self,input_line):

        """[用于实现对计算结果的压力单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """       

        UnitConv.pressure_output["帕斯卡（Pa）"] = float(input_line)
        UnitConv.pressure_output["千帕（kPa）"] = float(input_line) * 1e-3
        UnitConv.pressure_output["兆帕（MPa）"] = float(input_line) * 1e-6
        UnitConv.pressure_output["标准大气压（atm）"] = float(input_line) * 9.8692e-6
        UnitConv.pressure_output["巴（Bar）"] = float(input_line) * 0.00001

    def pressure_Unit_conversion(self,input_data,input_line):

        """[用于实现压力单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [压力单位选择框]
            input_line ([QLineEdit]): [压力数值编辑框]
        """   

        try:
            if input_data.currentText() == "帕斯卡（Pa）":
                UnitConv.pressure_output["帕斯卡（Pa）"]=float(input_line.text())
                UnitConv.pressure_output["千帕（kPa）"] = float(input_line.text())*1e-3
                UnitConv.pressure_output["兆帕（MPa）"] = float(input_line.text())*1e-6
                UnitConv.pressure_output["标准大气压（atm）"] = float(input_line.text())*9.8692e-6
                UnitConv.pressure_output["巴（Bar）"] = float(input_line.text()) *0.00001

            elif input_data.currentText()=="千帕（kPa）":
                UnitConv.pressure_output["帕斯卡（Pa）"]=float(input_line.text())*1000
                UnitConv.pressure_output["千帕（kPa）"] = float(input_line.text())
                UnitConv.pressure_output["兆帕（MPa）"] = float(input_line.text())*1000
                UnitConv.pressure_output["标准大气压（atm）"] = float(input_line.text())*0.0098692
                UnitConv.pressure_output["巴（Bar）"] = float(input_line.text()) *0.01

            elif input_data.currentText()=="兆帕（MPa）":
                UnitConv.pressure_output["帕斯卡（Pa）"]=float(input_line.text())*1e6
                UnitConv.pressure_output["千帕（kPa）"] = float(input_line.text())*1000
                UnitConv.pressure_output["兆帕（MPa）"] = float(input_line.text())
                UnitConv.pressure_output["标准大气压（atm）"] = float(input_line.text())*9.8692327
                UnitConv.pressure_output["巴（Bar）"] = float(input_line.text()) *10

            elif input_data.currentText()=="标准大气压（atm）":
                UnitConv.pressure_output["帕斯卡（Pa）"]=float(input_line.text())*101325
                UnitConv.pressure_output["千帕（kPa）"] = float(input_line.text())*101.325
                UnitConv.pressure_output["兆帕（MPa）"] = float(input_line.text())*0.101325
                UnitConv.pressure_output["标准大气压（atm）"] = float(input_line.text())
                UnitConv.pressure_output["巴（Bar）"] = float(input_line.text()) *1.01325

            elif input_data.currentText()=="巴（Bar）":
                UnitConv.pressure_output["帕斯卡（Pa）"]=float(input_line.text())*100000
                UnitConv.pressure_output["千帕（kPa）"] = float(input_line.text())*100
                UnitConv.pressure_output["兆帕（MPa）"] = float(input_line.text())*0.1
                UnitConv.pressure_output["标准大气压（atm）"] = float(input_line.text())*0.9869233
                UnitConv.pressure_output["巴（Bar）"] = float(input_line.text())

        except:
            pass


    def pressure_refresh (self,input_data,input_line):

        """[用于实现压力单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [压力单位选择框]
            input_line ([QLineEdit]): [压力数值编辑框]
        """     

        try:
            if input_data.currentText() == "巴（Bar）":
                input_line.setText(str(UnitConv.pressure_output["巴（Bar）"]))

            elif input_data.currentText()=="帕斯卡（Pa）":
                input_line.setText(str(UnitConv.pressure_output["帕斯卡（Pa）"]))

            elif input_data.currentText() == "千帕（kPa）":
                input_line.setText(str(UnitConv.pressure_output["千帕（kPa）"]))

            elif input_data.currentText() == "兆帕（MPa）":
                input_line.setText(str(UnitConv.pressure_output["兆帕（MPa）"]))

            elif input_data.currentText() == "标准大气压（atm）":
                input_line.setText(str(UnitConv.pressure_output["标准大气压（atm）"]))

        except:
            pass
    def _temperature_unit_conversion(self, input_line):

        """[用于实现对计算结果的温度单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """    

        UnitConv.temperature_output["摄氏度（℃）"] = float(input_line) - 273.15
        UnitConv.temperature_output["华氏度（F）"] = 32 + 1.8 * (float(input_line) - 273.15)
        UnitConv.temperature_output["开尔文（K）"] = float(input_line)

    def temperature_unit_conversion(self, input_data, input_line):

        """[用于实现温度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """ 

        try:
            if input_data.currentText() == "摄氏度（℃）":
                UnitConv.temperature_output["摄氏度（℃）"] = float(input_line.text())
                UnitConv.temperature_output["华氏度（F）"] = 32+1.8*float(input_line.text())
                UnitConv.temperature_output["开尔文（K）"] = float(input_line.text()) +273.15

            elif input_data.currentText() == "华氏度（F）":
                UnitConv.temperature_output["摄氏度（℃）"] = (float(input_line.text())-32)/1.8
                UnitConv.temperature_output["华氏度（F）"] =  float(input_line.text())
                UnitConv.temperature_output["开尔文（K）"] = (float(input_line.text())-32)/1.8+273.15

            elif input_data.currentText() == "开尔文（K）":
                UnitConv.temperature_output["摄氏度（℃）"] = float(input_line.text())-273.15
                UnitConv.temperature_output["华氏度（F）"] = 32 + 1.8 * (float(input_line.text())-273.15)
                UnitConv.temperature_output["开尔文（K）"] = float(input_line.text())

        except:
            pass

    def temperature_refresh(self, input_data, input_line):

        """[用于实现温度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """   

        try:
            if input_data.currentText() == "摄氏度（℃）":
                input_line.setText(str(UnitConv.temperature_output["摄氏度（℃）"]))

            elif input_data.currentText() == "华氏度（F）":
                input_line.setText(str(UnitConv.temperature_output["华氏度（F）"]))

            elif input_data.currentText() == "开尔文（K）":
                input_line.setText(str(UnitConv.temperature_output["开尔文（K）"]))

        except:
            pass

    def _rho_unit_conversion(self, input_line):

        """[用于实现对计算结果的密度单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """    

        UnitConv.rho_output["千克/立方米（kg/m3）"] = float(input_line)
        UnitConv.rho_output["克/立方厘米（g/cm3）"] = float(input_line) * 1000

    def rho_unit_conversion(self, input_data, input_line):

        """[用于实现密度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [密度单位选择框]
            input_line ([QLineEdit]): [密度数值编辑框]
        """ 

        try:
            if input_data.currentText() == "千克/立方米（kg/m3）":
                UnitConv.rho_output["千克/立方米（kg/m3）"] = float(input_line.text())
                UnitConv.rho_output["克/立方厘米（g/cm3）"] = float(input_line.text())*1000

            elif input_data.currentText() == "克/立方厘米（g/cm3）":
                UnitConv.rho_output["千克/立方米（kg/m3）"] = (float(input_line.text()))/1000
                UnitConv.rho_output["克/立方厘米（g/cm3）"] =  float(input_line.text())
        except:
            pass
    
    def rho_refresh(self, input_data, input_line):

        """[用于实现温度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """   

        try:
            if input_data.currentText() == "千克/立方米（kg/m3）":
                input_line.setText(str(UnitConv.rho_output["千克/立方米（kg/m3）"]))

            elif input_data.currentText() == "克/立方厘米（g/cm3）":
                input_line.setText(str(UnitConv.rho_output["克/立方厘米（g/cm3）"]))

        except:
            pass


def EOS_model():

    """[实现设置好参数后进行计算的功能]

    Raises:
        ValueError: [当三个参数全被激活时则报错，提醒需要选择计算的参数]
    """   

    try:
        if ui.PA_gas_name.currentText() == "氢气":
            sepcies="hydrogen"

        elif ui.PA_gas_name.currentText() == "氦气":
            sepcies="helium"

        EOSS=EOS(sepcies)

        if ui.PA_P_in.isEnabled() == False:
            output=str(EOSS.P(float(ui.PA_T_in.text()), float(ui.PA_rho_in.text())))
            conv._pressure_Unit_conversion(output)
            conv.pressure_refresh(ui.PA_P,ui.PA_P_in)

        if ui.PA_T_in.isEnabled() == False:
            output = str(EOSS.T(float(ui.PA_P_in.text()), float(ui.PA_rho_in.text())))
            conv._temperature_unit_conversion(output)
            conv.temperature_refresh(ui.PA_T, ui.PA_T_in)

        if ui.PA_rho_in.isEnabled() == False:
            output = str(EOSS.rho(float(ui.PA_T_in.text()), float(ui.PA_P_in.text())))
            conv._rho_unit_conversion(output)
            conv.rho_refresh(ui.PA_rho, ui.PA_rho_in)

        if ui.PA_rho_in.isEnabled() and ui.PA_T_in.isEnabled() and ui.PA_P_in.isEnabled() ==True:
            raise ValueError("请选择需要计算的参数")

    except ValueError as Argument:
        QMessageBox.critical(MainWindow,"Error", str(Argument).split(":")[0], QMessageBox.Yes, QMessageBox.Yes)

def NN_model():
    gas_name = ui.NN_gas_name.currentText()
    Dia = ui.NN_ambT_in
    Pres = ui.NN_P_in
    Tamb = ui.NN_ambT_in
    Tgas = ui.NN_gasT_in
    NN_model_name = ui.NN_model


#############################################主循环#################################################
if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
#############################################信号连接#################################################
    """选择参数的按钮，当按钮被激活时对应的输入栏变灰色"""
    ui.parameters_P.toggled.connect(lambda:radio_button(ui.parameters_P,ui.PA_P_in))
    ui.parameters_T.toggled.connect(lambda: radio_button(ui.parameters_T,ui.PA_T_in))
    ui.parameters_rho.toggled.connect(lambda: radio_button(ui.parameters_rho,ui.PA_rho_in))
    
    """实例化单位换算类"""
    conv = UnitConv()

    """直径单位换算及数值更新"""
    ui.NN_dia.currentIndexChanged.connect(lambda: conv.length_refresh(ui.NN_dia,ui.NN_dia_in))
    ui.NN_dia_in.editingFinished.connect(lambda :conv.length_Unit_conversion(ui.NN_dia,ui.NN_dia_in))
    
    """压力单位换算及数值更新"""
    ui.NN_P.currentIndexChanged.connect(lambda :conv.pressure_refresh(ui.NN_P,ui.NN_P_in))
    ui.NN_P_in.editingFinished.connect(lambda :conv.pressure_Unit_conversion(ui.NN_P,ui.NN_P_in))
    ui.PA_P.currentIndexChanged.connect(lambda :conv.pressure_refresh(ui.PA_P,ui.PA_P_in))
    ui.PA_P_in.editingFinished.connect(lambda :conv.pressure_Unit_conversion(ui.PA_P,ui.PA_P_in))
    
    """温度单位换算及数值更新"""
    ui.PA_T.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.PA_T,ui.PA_T_in))
    ui.PA_T_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.PA_T,ui.PA_T_in))
    ui.NN_ambT.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.NN_ambT,ui.NN_ambT_in))
    ui.NN_ambT_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.NN_ambT,ui.NN_ambT_in))
    ui.NN_gasT.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.NN_gasT,ui.NN_gasT_in))
    ui.NN_gasT_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.NN_gasT,ui.NN_gasT_in))
    
    """密度单位换算及数值更新"""
    ui.PA_rho.currentIndexChanged.connect(lambda :conv.rho_refresh(ui.PA_rho,ui.PA_rho_in))
    ui.PA_rho_in.editingFinished.connect(lambda :conv.rho_unit_conversion(ui.PA_rho,ui.PA_rho_in))
    
    """开始计算流体参数"""
    ui.pushButton.clicked.connect(lambda :EOS_model())
    #ui.NN_start_calc.clicked.connect()

    MainWindow.show()
    sys.exit(app.exec_())