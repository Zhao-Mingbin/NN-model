#############################################模块引进#################################################
from nn import Ui_MainWindow
from nn_therm import *
from Other_Nn import *
from PyQt5.QtWidgets import QFileDialog, QFormLayout, QLineEdit,QMessageBox
from PyQt5 import QtCore, QtGui, QtWidgets
import os
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

    global g_length_output
    global g_rho_output
    global PA_pressure_output
    global PA_temperature_output
    global Tamb_temperature_output
    global Tgas_temperature_output
    global NN_pressure_output

    g_length_output = {}

    g_rho_output = {}
    PA_temperature_output = {}
    Tamb_temperature_output = {}
    Tgas_temperature_output = {}

    PA_pressure_output = {}
    NN_pressure_output = {}




    def length_Unit_conversion(self,input_data,input_line,length_output):

        """[用于实现长度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [长度单位选择框]
            input_line ([QLineEdit]): [长度数值编辑框]
        """   

        try:
            if input_data.currentText() == "分米（dm）":
                length_output["分米（dm）"]=float(input_line.text())
                length_output["厘米（cm）"] = float(input_line.text())*10
                length_output["毫米（mm）"] = float(input_line.text())*100
                length_output["米（m）"] = float(input_line.text())/10

            elif input_data.currentText()=="厘米（cm）":
                length_output["分米（dm）"]=float(input_line.text())/10
                length_output["厘米（cm）"] = input_line.text()
                length_output["毫米（mm）"] = float(input_line.text())*10
                length_output["米（m）"] = float(input_line.text())/100

            elif input_data.currentText()=="毫米（mm）":
                length_output["分米（dm）"]=float(input_line.text())/100
                length_output["厘米（cm）"] = float(input_line.text())/10
                length_output["毫米（mm）"] = float(input_line.text())
                length_output["米（m）"] = float(input_line.text())/1000

            elif input_data.currentText()=="米（m）":
                length_output["分米（dm）"]=float(input_line.text())/10
                length_output["厘米（cm）"] =float (input_line.text())/100
                length_output["毫米（mm）"] = float(input_line.text())/1000
                length_output["米（m）"] = float(input_line.text())
        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的长度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global g_length_output
            g_length_output = {}
    def length_refresh(self,input_data,input_line,length_output):

        """[用于实现长度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [长度单位选择框]
            input_line ([QLineEdit]): [长度数值编辑框]
        """        

        try:
            if input_data.currentText() == "分米（dm）":
                input_line.setText(str(length_output["分米（dm）"]))

            elif input_data.currentText()=="厘米（cm）":
                input_line.setText(str(length_output["厘米（cm）"]))

            elif input_data.currentText() == "毫米（mm）":
                input_line.setText(str(length_output["毫米（mm）"]))

            elif input_data.currentText() == "米（m）":
                input_line.setText(str(length_output["米（m）"]))

        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的长度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global g_length_output
            g_length_output = {}
    def _pressure_Unit_conversion(self,input_line,pressure_output):

        """[用于实现对计算结果的压力单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """       
        try:
            pressure_output["帕斯卡（Pa）"] = float(input_line)
            pressure_output["千帕（kPa）"] = float(input_line) * 1e-3
            pressure_output["兆帕（MPa）"] = float(input_line) * 1e-6
            pressure_output["标准大气压（atm）"] = float(input_line) * 9.8692e-6
            pressure_output["巴（Bar）"] = float(input_line) * 0.00001
        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的压力计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global NN_pressure_output, PA_pressure_output
            NN_pressure_output = {}
            PA_pressure_output = {}
    def pressure_Unit_conversion(self,input_data,input_line,pressure_output):

        """[用于实现压力单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [压力单位选择框]
            input_line ([QLineEdit]): [压力数值编辑框]
        """   

        try:
            if input_data.currentText() == "帕斯卡（Pa）":
                pressure_output["帕斯卡（Pa）"]=float(input_line.text())
                pressure_output["千帕（kPa）"] = float(input_line.text())*1e-3
                pressure_output["兆帕（MPa）"] = float(input_line.text())*1e-6
                pressure_output["标准大气压（atm）"] = float(input_line.text())*9.8692e-6
                pressure_output["巴（Bar）"] = float(input_line.text()) *0.00001

            elif input_data.currentText()=="千帕（kPa）":
                pressure_output["帕斯卡（Pa）"]=float(input_line.text())*1000
                pressure_output["千帕（kPa）"] = float(input_line.text())
                pressure_output["兆帕（MPa）"] = float(input_line.text())/1000
                pressure_output["标准大气压（atm）"] = float(input_line.text())*0.0098692
                pressure_output["巴（Bar）"] = float(input_line.text()) *0.01

            elif input_data.currentText()=="兆帕（MPa）":
                pressure_output["帕斯卡（Pa）"]=float(input_line.text())*1e6
                pressure_output["千帕（kPa）"] = float(input_line.text())*1000
                pressure_output["兆帕（MPa）"] = float(input_line.text())
                pressure_output["标准大气压（atm）"] = float(input_line.text())*9.8692327
                pressure_output["巴（Bar）"] = float(input_line.text()) *10

            elif input_data.currentText()=="标准大气压（atm）":
                pressure_output["帕斯卡（Pa）"]=float(input_line.text())*101325
                pressure_output["千帕（kPa）"] = float(input_line.text())*101.325
                pressure_output["兆帕（MPa）"] = float(input_line.text())*0.101325
                pressure_output["标准大气压（atm）"] = float(input_line.text())
                pressure_output["巴（Bar）"] = float(input_line.text()) *1.01325

            elif input_data.currentText()=="巴（Bar）":
                pressure_output["帕斯卡（Pa）"]=float(input_line.text())*100000
                pressure_output["千帕（kPa）"] = float(input_line.text())*100
                pressure_output["兆帕（MPa）"] = float(input_line.text())*0.1
                pressure_output["标准大气压（atm）"] = float(input_line.text())*0.9869233
                pressure_output["巴（Bar）"] = float(input_line.text())


        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的压力计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global NN_pressure_output, PA_pressure_output
            NN_pressure_output = {}
            PA_pressure_output = {}

    def pressure_refresh (self,input_data,input_line, pressure_output):

        """[用于实现压力单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [压力单位选择框]
            input_line ([QLineEdit]): [压力数值编辑框]
        """     

        try:
            if input_data.currentText() == "巴（Bar）":
                input_line.setText(str(pressure_output["巴（Bar）"]))

            elif input_data.currentText()=="帕斯卡（Pa）":
                input_line.setText(str(pressure_output["帕斯卡（Pa）"]))

            elif input_data.currentText() == "千帕（kPa）":
                input_line.setText(str(pressure_output["千帕（kPa）"]))

            elif input_data.currentText() == "兆帕（MPa）":
                input_line.setText(str(pressure_output["兆帕（MPa）"]))

            elif input_data.currentText() == "标准大气压（atm）":
                input_line.setText(str(pressure_output["标准大气压（atm）"]))

        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的压力计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global NN_pressure_output,PA_pressure_output
            NN_pressure_output = {}
            PA_pressure_output = {}


    def _temperature_unit_conversion(self, input_line,temperature_output):

        """[用于实现对计算结果的温度单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """    
        try:
            temperature_output["摄氏度（℃）"] = float(input_line) - 273.15
            temperature_output["华氏度（F）"] = 32 + 1.8 * (float(input_line) - 273.15)
            temperature_output["开尔文（K）"] = float(input_line)
        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的温度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global PA_temperature_output, Tamb_temperature_output, Tgas_temperature_output
            PA_temperature_output = {}
            Tamb_temperature_output = {}
            Tgas_temperature_output = {}
    def temperature_unit_conversion(self, input_data, input_line,temperature_output):

        """[用于实现温度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """ 

        try:
            if input_data.currentText() == "摄氏度（℃）":
                temperature_output["摄氏度（℃）"] = float(input_line.text())
                temperature_output["华氏度（F）"] = 32+1.8*float(input_line.text())
                temperature_output["开尔文（K）"] = float(input_line.text()) +273.15

            elif input_data.currentText() == "华氏度（F）":
                temperature_output["摄氏度（℃）"] = (float(input_line.text())-32)/1.8
                temperature_output["华氏度（F）"] =  float(input_line.text())
                temperature_output["开尔文（K）"] = (float(input_line.text())-32)/1.8+273.15

            elif input_data.currentText() == "开尔文（K）":
                temperature_output["摄氏度（℃）"] = float(input_line.text())-273.15
                temperature_output["华氏度（F）"] = 32 + 1.8 * (float(input_line.text())-273.15)
                temperature_output["开尔文（K）"] = float(input_line.text())

        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的温度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global PA_temperature_output, Tamb_temperature_output, Tgas_temperature_output
            PA_temperature_output = {}
            Tamb_temperature_output = {}
            Tgas_temperature_output = {}

    def temperature_refresh(self, input_data, input_line, temperature_output):

        """[用于实现温度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """   

        try:
            if input_data.currentText() == "摄氏度（℃）":
                input_line.setText(str(temperature_output["摄氏度（℃）"]))

            elif input_data.currentText() == "华氏度（F）":
                input_line.setText(str(temperature_output["华氏度（F）"]))

            elif input_data.currentText() == "开尔文（K）":
                input_line.setText(str(temperature_output["开尔文（K）"]))

        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的温度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global PA_temperature_output, Tamb_temperature_output, Tgas_temperature_output
            PA_temperature_output = {}
            Tamb_temperature_output = {}
            Tgas_temperature_output = {}
    def _rho_unit_conversion(self, input_line,rho_output):

        """[用于实现对计算结果的密度单位转换]

        Args:
            input_line ([str]): [CoolProp的输出结果]
        """    
        try:
            rho_output["千克/立方米（kg/m3）"] = float(input_line)
            rho_output["克/立方厘米（g/cm3）"] = float(input_line) * 1000
        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的密度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global g_rho_output
            g_rho_output = {}
    def rho_unit_conversion(self, input_data, input_line,rho_output):

        """[用于实现密度单位换算，生成对应的{单位：数值}字典，实现单位换算]

        Args:
            input_data ([QComboBox]): [密度单位选择框]
            input_line ([QLineEdit]): [密度数值编辑框]
        """ 

        try:
            if input_data.currentText() == "千克/立方米（kg/m3）":
                rho_output["千克/立方米（kg/m3）"] = float(input_line.text())
                rho_output["克/立方厘米（g/cm3）"] = float(input_line.text())*1000

            elif input_data.currentText() == "克/立方厘米（g/cm3）":
                rho_output["千克/立方米（kg/m3）"] = (float(input_line.text()))/1000
                rho_output["克/立方厘米（g/cm3）"] =  float(input_line.text())
        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的密度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global g_rho_output
            g_rho_output = {}
    def rho_refresh(self, input_data, input_line,rho_output):

        """[用于实现温度单位换算后更新编辑框内的数值]

        Args:
            input_data ([QComboBox]): [温度单位选择框]
            input_line ([QLineEdit]): [温度数值编辑框]
        """   

        try:
            if input_data.currentText() == "千克/立方米（kg/m3）":
                input_line.setText(str(rho_output["千克/立方米（kg/m3）"]))

            elif input_data.currentText() == "克/立方厘米（g/cm3）":
                input_line.setText(str(rho_output["克/立方厘米（g/cm3）"]))

        except:
            QMessageBox.critical(MainWindow, "错误", '请输入合适的密度计算参数', QMessageBox.Yes, QMessageBox.Yes)
            global g_rho_output
            g_rho_output = {}

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

        EOSS = EOS(sepcies)

        if ui.PA_P_in.isEnabled() == False:

            output = str(EOSS.P(float(PA_temperature_output["开尔文（K）"]), float(g_rho_output["千克/立方米（kg/m3）"])))
            conv._pressure_Unit_conversion(output,PA_pressure_output)
            conv.pressure_refresh(ui.PA_P,ui.PA_P_in, PA_pressure_output)

        if ui.PA_T_in.isEnabled() == False:

            output = str(EOSS.T(float(PA_pressure_output["帕斯卡（Pa）"]), float(g_rho_output["千克/立方米（kg/m3）"])))
            conv._temperature_unit_conversion(output,PA_temperature_output)
            conv.temperature_refresh(ui.PA_T, ui.PA_T_in,PA_temperature_output)

        if ui.PA_rho_in.isEnabled() == False:

            output = str(EOSS.rho(float(PA_temperature_output["开尔文（K）"]), float(PA_pressure_output["帕斯卡（Pa）"])))
            conv._rho_unit_conversion(output,g_rho_output)
            conv.rho_refresh(ui.PA_rho, ui.PA_rho_in,g_rho_output)

        if ui.PA_rho_in.isEnabled() and ui.PA_T_in.isEnabled() and ui.PA_P_in.isEnabled() ==True:
            raise ValueError("请选择需要计算的参数")

    except ValueError as Argument:
        QMessageBox.critical(MainWindow,"错误", str(Argument).split(":")[0], QMessageBox.Yes, QMessageBox.Yes)
    except KeyError:
        QMessageBox.critical(MainWindow, "错误", '请输入合适的计算参数', QMessageBox.Yes, QMessageBox.Yes)
def printf(output_message):  # 用于向屏幕输出信息，可以替代print函数
    ui.NN_textbroswer.append(output_message)  # 在指定的区域显示提示信息
    ui.cursor = ui.NN_textbroswer.textCursor()
    ui.NN_textbroswer.moveCursor(ui.cursor.End)  # 光标移到最后，这样就会自动显示出来
    QtWidgets.QApplication.processEvents()  # 一定加上这个功能，不然有卡顿

def NN_model():
    try:
        gas_name = ui.NN_gas_name.currentText()
        Dia = float(g_length_output["米（m）"])
        Pres = float(NN_pressure_output["帕斯卡（Pa）"])
        Tamb = float(Tamb_temperature_output["开尔文（K）"])
        Tgas = float(Tgas_temperature_output["开尔文（K）"])
        NN_model_name = ui.NN_model.currentText()
        if gas_name == '氢气':
            species = 'hydrogen'
        if gas_name == '氦气':
            species = 'helium'
        therm_gas = EOS(species = species)
        therm_air = EOS(species = "air")
        gas_high = Gas(therm_gas, T=Tgas, P=Pres)
        gas_low = Gas(therm_air, T=Tamb, P=1.01325e5)
        orifice = Orifice(d=Dia)
        nn = NotionalNozzle(gas_high, orifice, gas_low)
        # print(type(NN_model_name),NN_model_name)
        # print(nn.calculate())#model=NN_model_name
        class effective:
            def __init__(self, gas_eff=None, orifice_eff=None, Veff=None):
                '''class for effective source'''
                self.gas_eff, self.orifice_eff, self.Veff = gas_eff, orifice_eff, Veff
                self.rho = gas_eff.rho
                self.mdot, self.d = orifice_eff.mdot(self.rho, Veff), orifice_eff.d
                self.P, self.T = gas_eff.P, gas_eff.T

        result = effective(*(nn.calculate(model=NN_model_name)))
        global string
        string = f'''
    ******** { NN_model_name} ********

           ****输入计算参数**** 

    * 环境温度（{ui.NN_ambT.currentText().split('（')[1].split('）')[0]}）= {ui.NN_ambT_in.text()} 

    * 气体温度（{ui.NN_gasT.currentText().split('（')[1].split('）')[0]}）= {ui.NN_gasT_in.text()} 

    * 储存压力（{ui.NN_P.currentText().split('（')[1].split('）')[0]}）= {ui.NN_P_in.text()} 

    * 出口直径（{ui.NN_dia.currentText().split('（')[1].split('）')[0]} = {ui.NN_dia_in.text()} 

        *****虚喷嘴计算参数***** 

    * 出口速度（m/s） = {result.Veff:5.4e}    

    * 质量变化率（kg/s）= {result.mdot:5.4e}  

    * 虚喷嘴半径（m）= {result.d / 2.:5.4e}   

    * 出口温度（K）= {result.T:5.4e}      

    * 出口压力（Pa）= {result.P:5.4e}   

    *******************************
    
        '''
        printf(string)
        QMessageBox.information(MainWindow, "提示", "计算成功！", QMessageBox.Yes, QMessageBox.Yes)
    except KeyError:
        QMessageBox.critical(MainWindow, "错误", "请设置合适的计算参数", QMessageBox.Yes, QMessageBox.Yes)
    except ValueError as Arguement:
        QMessageBox.critical(MainWindow, "错误", str(Arguement).split(":")[0], QMessageBox.Yes, QMessageBox.Yes)
    except RuntimeError:
        QMessageBox.critical(MainWindow, "错误", "程序将陷入死循环，请重新选择参数", QMessageBox.Yes, QMessageBox.Yes)



#  to calculate melting line T(p) 储存压力 Cannot find good Tmin #Inputs in Brent [0.000000,46524.000000] do not bracket the root.

#Failed to converge after 50 iterations, 气体温度 argument index out of range









def save_txt():
    try:
        file_name = f"""{ui.NN_gas_name.currentText()}_{ui.NN_gasT_in.text()}{ui.NN_ambT.currentText().split('（')[1].split('）')[0]}_{ui.NN_dia_in.text()}{ui.NN_dia.currentText().split('（')[1].split('）')[0]}_{ui.NN_P_in.text()}{ui.NN_P.currentText().split('（')[1].split('）')[0]}.txt"""
        file_path=os.getcwd()
        if len(string)<1:
            raise ValueError
        with open (file_name,'w') as f:
            f.write(string)
        printf(f'文件 "{file_name}" 已成功保存于 "{file_path}".')
        QMessageBox.information(MainWindow, "提示", "文件保存成功！", QMessageBox.Yes, QMessageBox.Yes)
    except:
        QMessageBox.critical(MainWindow, "错误", '无计算结果', QMessageBox.Yes, QMessageBox.Yes)
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
    ui.NN_dia.currentIndexChanged.connect(lambda: conv.length_refresh(ui.NN_dia,ui.NN_dia_in,g_length_output))
    ui.NN_dia_in.editingFinished.connect(lambda :conv.length_Unit_conversion(ui.NN_dia,ui.NN_dia_in,g_length_output))
    
    """压力单位换算及数值更新"""
    ui.NN_P.currentIndexChanged.connect(lambda :conv.pressure_refresh(ui.NN_P, ui.NN_P_in, NN_pressure_output))
    ui.NN_P_in.editingFinished.connect(lambda :conv.pressure_Unit_conversion(ui.NN_P, ui.NN_P_in, NN_pressure_output))
    ui.PA_P.currentIndexChanged.connect(lambda :conv.pressure_refresh(ui.PA_P, ui.PA_P_in, PA_pressure_output))
    ui.PA_P_in.editingFinished.connect(lambda :conv.pressure_Unit_conversion(ui.PA_P, ui.PA_P_in, PA_pressure_output))
    
    """温度单位换算及数值更新"""
    ui.PA_T.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.PA_T, ui.PA_T_in, PA_temperature_output))
    ui.PA_T_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.PA_T, ui.PA_T_in, PA_temperature_output))
    ui.NN_ambT.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.NN_ambT, ui.NN_ambT_in, Tamb_temperature_output))
    ui.NN_ambT_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.NN_ambT, ui.NN_ambT_in, Tamb_temperature_output))
    ui.NN_gasT.currentIndexChanged.connect(lambda :conv.temperature_refresh(ui.NN_gasT, ui.NN_gasT_in, Tgas_temperature_output))
    ui.NN_gasT_in.editingFinished.connect(lambda :conv.temperature_unit_conversion(ui.NN_gasT, ui.NN_gasT_in, Tgas_temperature_output))
    
    """密度单位换算及数值更新"""
    ui.PA_rho.currentIndexChanged.connect(lambda :conv.rho_refresh(ui.PA_rho, ui.PA_rho_in, g_rho_output))
    ui.PA_rho_in.editingFinished.connect(lambda :conv.rho_unit_conversion(ui.PA_rho, ui.PA_rho_in, g_rho_output))
    
    """开始计算流体参数"""
    ui.pushButton.clicked.connect(lambda :EOS_model())

    ui.NN_start_calc.clicked.connect(lambda :NN_model())

    ui.NN_output.clicked.connect(lambda :save_txt())

    MainWindow.show()
    sys.exit(app.exec_())