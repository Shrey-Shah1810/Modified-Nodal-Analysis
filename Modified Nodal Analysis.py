import sys
import traceback  # for command line arguments and exiting
import numpy as np  # for matrix mathematics for MNA
import cmath # complex mathematics for AC analysis
import math

ckt_start = ".circuit"
ckt_end = ".end"
AC=".ac"
ac_checker=False
ground_node="GND"

dict_elem=dict()
node_list=[]

# Decodes the value of string
def value_decoder(value):
    num_val=0
    str_val=""
    str_ticker=0
    for ind,i in enumerate(value):
        if i.isalpha():
            num_val=float(value[:ind])
            str_val=value[ind:]
            str_ticker=1
        else:
            if ind==(len(value)-1):
                return float(value)

    # split = re.split("([a-z]+)", value.lower())
    if len(str_val)>0:
        if str_val== "e":  # no suffix or exponent notation
            return float(value)
        str_val=str_val.lower()

        if str_val[0]=="k":  # kilo(*1000)
            return num_val * 1e3
        elif str_val[:3]=="meg":  # mega(*1000000)
            return num_val * 1e6
        elif str_val[0]=="g":  # giga(*1000000000)
            return num_val * 1e9
        elif str_val[0]=="t":  # tera(*1000000000000)
            return num_val * 1e12
        elif str_val[0]=="f":  # femto(*0.000000000000001)
            return num_val * 1e-15
        elif str_val[0]=="p":  # pico(*0.000000000001)
            return num_val * 1e-12
        elif str_val[0]=="n":  # nano(*0.000000001)
            return num_val * 1e-9
        elif str_val[0]=="u" or str_val[-5:]=="micro": #micro(*0.000001)
            return num_val * 1e-6
        elif str_val[0]==("m"):  # milli(*0.001)
            return num_val * 1e-3

class parent_Element():
    """
    Adding the basic elements common to every element in a common parent class
    """
    def __init__(self, line):
        self.name = line[0]
        self.type = self.name[0]

        self.n1 = line[1]
        self.n2 = line[2]
        if self.n1 not in node_list:
            node_list.append(self.n1)
        if self.n2 not in node_list:
            node_list.append(self.n2)
        self.value = value_decoder(line[-1])

class Passive_Element(parent_Element):

    def __init__(self,line):
        super().__init__(line)
        if self.name not in dict_elem:
            dict_elem[self.name]=self.value        

class Ind_Src(parent_Element):
    def __init__(self, line, ac_checker):
        super().__init__(line)
        if ac_checker:
            self.pp=line[4]
            self.phase=line[5]
        if self.name not in dict_elem:
            dict_elem[self.name]=self.value

class Dependent_Src(parent_Element):
    def __init__(self, line):
        super().__init__(line)
        if ac_checker:
            self.current = True
            if(len(line) == 6):
                self.node3 = line[3]
                self.node4 = line[4]
                self.current = False
            else:
                self.cVoltage = line[3]
        else:
            self.current=True
            if len(line)==6:
                self.n3=line[3]
                self.n4=line[4]
                self.current = False
            else:
                self.cVoltage = line[3]
                
        if self.name not in dict_elem:
            dict_elem[self.name]=self.value

num_variables = 0
voltage_dict = dict()
voltages = []
elements = []

try:
    with open(sys.argv[1], "r") as f:
        text = f.read().splitlines()
except FileNotFoundError as FNFE:
    print(FNFE)
    sys.exit(1)

try:
    start_idx = text.index(ckt_start) + 1
    end_idx = text.index(ckt_end)

    assert start_idx > 0 and start_idx < end_idx

    if end_idx + 1 != len(text):
        if text[end_idx + 1].split()[0] == AC:
            ac_checker = True  # ac_checker is enabled
            W = 2 * math.pi * float(str(value_decoder(text[end_idx + 1].split()[-1])))

    elements = []

    for line in text[start_idx:end_idx]:
        tokens = line.split("#")[0].split()
        if tokens[0][0] in ["R", "L", "C"]:
            element = Passive_Element(tokens)

        elif tokens[0][0] in ["V", "I"]:
            element = Ind_Src(tokens, ac_checker)
            if element.name[0] == "V":
                if ac_checker:
                    voltage_dict[element.name] = [element.n1,element.n2, element.value, element.pp, element.phase,]
                else:
                    voltage_dict[element.name] = [element.n1, element.n2, element.value]
                num_variables += 1
                voltages.append(element.name)

        else:
            element = Dependent_Src(tokens)
            if element.current:
                if element.name[0] == "H":
                    voltages.append(element.name)
                    num_variables += 1
            else:
                pass

                if element.name[0] == "E":
                    voltages.append(element.name)
                    num_variables += 1

        elements.append(element)
except (ValueError, TypeError) as e:
    sys.exit(f"Error: Invalid value in circuit file {traceback.print_exc()})")
except AssertionError:
    sys.exit("Error: Invalid circuit file")

num_variables += len(node_list)
if ac_checker:
    S = np.zeros((num_variables + 1, 1), dtype=complex)
    M = np.zeros((num_variables + 1, num_variables + 1), dtype=complex)

else:
    S = np.zeros((num_variables + 1, 1))
    M = np.zeros((num_variables + 1, num_variables + 1))

node_tracer = {"GND": 0} # Assuming GND as ground
GND_Flag = False
for node in node_list:
    if node not in node_tracer.keys():
        node_tracer[node] = len(node_tracer)
    if node == "GND":
        GND_Flag = True
        
assert GND_Flag, " Please label ground as GND"


"""   
CHANGE CONTENT BELOW
"""
# for the ground node, we add only 1 without subtracting -1 as GND voltage is taken as zero
M[-1][0] += 1
M[0][-1] += 1
try:
    for element in elements:
        if element.name[0] == "R":
            conductance = 1 / float(element.value)
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            M[n1][n1] += conductance
            M[n2][n2] += conductance
            M[n1][n2] -= conductance
            M[n2][n1] -= conductance

        elif element.name[0] == "L":
            conductance = complex(0, -1 / (W * float(element.value)))
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            M[n1][n1] += conductance
            M[n2][n2] += conductance
            M[n1][n2] -= conductance
            M[n2][n1] -= conductance

        elif element.name[0] == "C":
            conductance = complex(0, W * float(element.value))
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            M[n1][n1] += conductance
            M[n2][n2] += conductance
            M[n1][n2] -= conductance
            M[n2][n1] -= conductance

        elif element.name[0] == "V":
            if ac_checker:
                n1 = node_tracer[element.n1]
                n2 = node_tracer[element.n2]
                In1n2 = voltages.index(element.name) + len(node_list)
                phase = math.radians(float(element.phase))
                A = float(element.pp) / 2
                M[n1][In1n2] += 1
                M[n2][In1n2] -= 1
                M[In1n2][n1] += 1
                M[In1n2][n2] -= 1
                z = complex(A * math.cos(phase), A * math.sin(phase))
                S[In1n2] += z

            else:
                n1 = node_tracer[element.n1]
                n2 = node_tracer[element.n2]
                In1n2 = voltages.index(element.name) + len(node_list)
                M[n1][In1n2] += 1
                M[n2][In1n2] -= 1
                M[In1n2][n1] += 1
                M[In1n2][n2] -= 1

                S[In1n2] += float(element.value)

        elif element.name[0] == "I":
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            S[n1] -= float(element.value)
            S[n2] += float(element.value)

        elif element.name[0] == "G":
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            n3 = node_tracer[element.n3]
            n4 = node_tracer[element.n4]
            conductance = float(element.value)
            M[n1][n3] += conductance
            M[n1][n4] += -conductance
            M[n2][n4] += conductance
            M[n2][n3] += -conductance

        elif element.name[0] == "E":
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            n3 = node_tracer[element.n3]
            n4 = node_tracer[element.n4]
            In1n2 = voltages.index(element.name) + len(node_list)

            M[n1][In1n2] += 1
            M[n2][In1n2] -= 1
            M[In1n2][n1] += 1
            M[In1n2][n2] -= 1
            M[In1n2][n3] -= float(element.value)
            M[In1n2][n4] += float(element.value)

        elif element.name[0] == "H":
            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            value = float(element.value)
            source = voltage_dict[element.cVoltage]
            n3 = node_tracer[source[0]]
            n4 = node_tracer[source[1]]
            source_value = float(source[2])
            In1n2 = voltages.index(element.name) + len(node_list)
            In3n4 = voltages.index(element.cVoltage) + len(node_list)
            M[n1][In1n2] += 1
            M[n2][In1n2] -= 1
            M[In1n2][n1] += 1
            M[In1n2][n2] -= 1
            M[In1n2][In3n4] -= value

        elif element.name[0] == "F":

            n1 = node_tracer[element.n1]
            n2 = node_tracer[element.n2]
            value = float(element.value)
            source = voltage_dict[element.cVoltage]
            n3 = node_tracer[source[0]]
            n4 = node_tracer[source[1]]
            source_value = float(source[2])
            In3n4 = voltages.index(element.cVoltage) + len(node_list)
            M[n1][In3n4] += value
            M[n2][In3n4] -= value
    x = np.matmul(np.linalg.inv(M), S)
except np.linalg.LinAlgError:
    sys.exit("The elements of circuit defined are not solvable")
except ValueError:
    print("The values provided are in an invalid format")
    sys.exit()

if ac_checker:
    for index in range(num_variables + 1):
        if index < len(node_list):
            property_elem = cmath.polar(x[node_tracer[node_list[index]]][0])
            print(
                "Voltage at the node {0} is: {1:0.5f} ∠ ({2:0.4f}°)".format(
                    node_list[index], property_elem[0], property_elem[1]*180/math.pi))
        elif index < len(node_list) + len(voltages):
            property_elem = cmath.polar(x[index][0])
            print(f"Current through the Voltage source in {voltages[index - len(node_list)]}",
                "is: {0:0.5f} ∠({1:0.4f}°)".format(property_elem[0], property_elem[1]*180/math.pi),)

else:
    for index in range(num_variables + 1):
        if index < len(node_list):
            print("Voltage at the node {0} is: {1:0.5f}".format(node_list[index], x[index][0]))
        elif index < len(node_list) + len(voltages):
            print(f"Current through the Voltage source is {voltages[index - len(node_list)]}",
                "is: {0:0.5f}".format(x[index][0]),)