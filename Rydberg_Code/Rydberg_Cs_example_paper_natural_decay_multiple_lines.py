# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 10:32:42 2021
last edited on 14/12/2021
@Author: Liam W Bussey
Affiliation@ BT and University of Birmingham

This code will model the steady state and transient 
evolution of a Rydberg reciever using ladder 
configuration for atomic excitation. 

The aim of this program is the give experimental
prediction based on user inputs from the state the user
wishes to study and the interacting field parameter to
make these tranistion. 

when considering the time evolution of this system we 
will study what happend when the RF/laser field is switch
on and off. the hope is to be able to effectively model a
realisitic receiver using the principle of quantum mechanics

In this code we ust the python module delevoped by durham 
know as arc -Alkali-Rydberg-Calculator, which has a large
library of information relvent to the group 1 alkali metals
more information can be found at 
http://pypi.org/project/ARC-Alkali-Rydberg-Calculator

To preform the steady state and transient analysis we use
the Qutip library. More information from this library can
be found at https://qutip.org/docs/latest/index.html
"""
from calendar import c
from cmath import pi
from re import I
import sys
from typing import Type
from numpy.lib.function_base import append
from qutip.qobj import Qobj
from brokenaxes import brokenaxes
from typing import Union

sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory') 
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory/src')
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory/src/python')

import numpy as np
    #pip install numpy
import pandas as pd
    #pip install pandas
import qutip as qp
    #pip install qutip
    #https://qutip.org/docs/latest/index.html
import scipy 
from scipy import signal as sn
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit, root_scalar

from scipy.constants import k as k_b
from scipy.constants import c as c_c ## speed of light 
from scipy.constants import h as c_h # planck constant
from scipy.constants import hbar as hbar #reduced planks

import Transition_search as transistion
import atom_transition_wavelen as atom_wave
import EIT_to_main as EITSS  #EIT_to_main_doppler

import arc 
    #pip install ARC-Alkali-Rydberg-Calculator
    #https://pypi.org/project/ARC-Alkali-Rydberg-Calculator/

from arc import Potassium41 as K41
from arc import Rubidium as Rb
from arc import Rubidium85 as Rb85
from arc import Rubidium87 as Rb87
from arc import Caesium as Cs
from arc import Strontium88 as Sr



import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes, inset_axes
from matplotlib.transforms import IdentityTransform
#from matplotlib.gridspec import GridSpec
    #https://matplotlib.org/stable/tutorials/introductory/pyplot.html
    #conda install matplotlib
import liams_windows_config

#import fsv_data_extract_moku as expdata

from scipy.constants import k as k_b # boltzman constant
from scipy.constants import c as c_c ## speed of light 
from scipy.constants import h as c_h # planck constant
from scipy.constants import hbar as hbar #reduced planks

from scipy.signal import chirp, find_peaks, peak_widths
##### figure settings #########
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size']= 20
plt.rcParams['figure.figsize'] = [14,12]
plt.rcParams['figure.autolayout']= True
plt.rcParams['axes.labelsize']= 18
plt.rcParams['axes.titlesize']= 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.titlesize'] = 20

class laser_parameter:
    """The laser and radio frequency wave parameters. These are the key experimental parameters, such as power, polarisation, wavelength etc. 
    """
    def set_feild(self, f: str) -> None:
        """The type of EM field to be described. 

        Args:
            f (str): If the field is a laser field, f will be set to 'L', and if f is a radio wave field, then 'R'.
        """
        OAM_list = ['L','R']
        while True:
            if f != 'Unset':
                self.feild = str.upper(f)
                break
            else:   
                my_input = input('please enter L for laser or R for Radio field: ')
                try:
                    my_input = str.upper(my_input)
                    if (my_input not in OAM_list):
                        print("That makes no sense... it is not in "
                              + OAM_list)
                        continue
                    else:
                        self.feild = my_input
                        break
                except: 
                    print("That makes no sense...")
                    continue
    
    def set_temp(self, temp: float) -> None:
        """Sets the operating temperature of the atoms under investigation.

        Args:
            temp (float): The temperature in K.
        """
        while True:
            if temp != 'Unset':
                self.temp = float(temp)
                break
            else:   
                my_input = input('please enter the temperature in K')
                
                self.temp = float(my_input)
                break

    def set_MW_feild(self,feild: str, vm: float) -> None:
        """Sets the radio wave electric field strength. This function will return np.Nan if the laser field type is 'L'.

        Args:
            feild (str): Field type. This should be 'R' as it will relate to the radio wave field strength.
            vm (float): The radio wave field strength in V m^(-1).
        """
        if (feild == 'R'):
            if vm != 'Unset':
                self.p_rf = float(vm)
            else:           
                my_input = float(input('enter maxium field strength in V/m?: ')) 

                self.p_rf = my_input       
        else:
            self.p_rf = np.nan

    def set_power(self, feild: str, laser_p: float) -> None:
        """Sets the laser power. This function will return np.Nan if the radio field type is 'R'.

        Args:
            feild (str): Field type. This should be 'L' as it will relate to the laser power.
            laser_p (float): The laser power in mW.
        """
        if (feild == 'L'):
            if laser_p != 'Unset':
                self.p = float(laser_p/1000) # laser power in W
            else: 
                my_input=float(input('Please enter laser power in mW: '))
                self.p = my_input/1000 # laser power in W
        else:
            self.p = np.nan

    def set_waist(self, feild: str, waist: float) -> None:
        """Sets the beam waist for the laser

        Args:
            feild (str): Field type. This should be 'L' as it will relate to the laser power.
            waist (float): The beam waist diameter in um.
        """
        if (feild == 'L'):
            if waist != 'Unset':
                self.w = float(waist*10**(-6))# Micrometers
            else: 
                my_input=float(input(r'Please enter laser waist in $\mu$ M: '))
                self.w  = my_input*10**(-6) # micrometers
        else:
            self.w = np.nan
    
    def set_q(self, q_set: int) -> None:
        """Sets the polarisation of the EM field, either radio or laser.

        Args:
            q_set (int): The polarisation state of the interacting field. If '-1', the field in clockwise circulary polarised (-sigma transition), if '0' the field is linearly polarised (pi transition), and if '1' the field is anti-clockwise polarised (+sigma transition)
        """
        OAM_list =['-1', '0', '1']
        while True:
            if q_set != 'Unset':
                self.q = int(q_set)
                break
            else:
                my_input = input(r'please enter polaraisation q = -1,0,1(-$\sigma$ ,$\pi$, +$\sigma$): ' )
                try:
                    if (my_input not in OAM_list):
                        print('not possible polaristion. enter value -1, 0 or 1')
                        continue
                    else:
                        self.q = int(my_input)
                        break
                except: 
                    print("That makes no sense...")
                    continue    

    def set_onress(self, ress: str) -> None:
        """Sets the interacting field to be on resonance with the atomic transistion it is probing. This is a 'Y/N' option.

        Args:
            ress (str): 'Y' if the transition is on resonance, and 'N' if the transition is not on resonance.
        """
        OAM_list =['Y','N']
        while True: 
            if ress != 'Unset':
                if (ress not in OAM_list):
                    print("That makes no sense... it is not in " + OAM_list)
                    break
                else:
                    self.onress = str.upper(ress)
                    break

            else:  
                my_input = input('Is perfecttly on resonance [Y/N]: ')
                try:
                    my_input = str.upper(my_input)
                    if (my_input not in OAM_list):
                        print("That makes no sense... it is not in " + OAM_list)
                        continue
                    else:
                        self.onress = my_input
                        break
                except: 
                    print("That makes no sense...")
                    continue        
    
    def set_scanning_source(self, scan: str ,onress: str) -> None:
        """Sets which field to scan over when inspecting the atomic transition. Note, a field that is being scanned over will not be able to be set to on resonance. If the field is scanned over, this option overrides the onress option. If the field is not on resonance and not scanning, the upper and lower dephasing limits must be set to be equal. 

        Args:
            scan (str): 'Y' if the field is going to be changed/inspected. If 'N' the field is assumed to bea fixed wavelength.
            onress (str): 'Y' if the field is on resonance with the atomic transition, else 'N'.
        """
        OAM_list =['Y','N']
        while True:
            if (onress == 'Y'):
                self.scanning_source = 'N'
                break
            else:
                if scan != 'Unset' and scan != 'N':
                    #print(scanning_source)
                    self.scanning_source = 'Y'
                    break
                elif scan != 'Unset' and scan != 'Y':
                    self.scanning_source = 'N'
                    break
                else:
                    my_input = input('Is this the source your scanning with? [Y/N]:')
                    try:
                        my_input = str.upper(my_input) 
                        
                        if (my_input not in OAM_list):
                            print("That makes no sense... it is not in " +
                                  OAM_list)
                            continue
                        else:
                            self.scanning_source = my_input
                            
                            break
                    except: 
                        print("That makes no sense...")
                        continue 

    def set_dphase_l(self, onress: str, deph_l: float) -> None:
        """The lower bound of the dephasing (distance from on resonance). If the transition is on resonace, the dephasing is overwritten to equal 0.

        Args:
            onress (str): The option of whether the transition is fixed to the atomic transition.
            deph_l (float): The lower distance from the resonant value, in MHz.
        """
        if (onress == 'Y'):
            self.dephase_lower = 0
        else:
            if deph_l != 'Unset':
                self.dephase_lower = float(deph_l*(10**(6))) ## MHz
                #print(self.dephase_lower)
            else:
                myinput = float(input('please enter lower dephasing value in MHz:'))
                self.dephase_lower = myinput*(10**(6))## MHz

    def set_dphase_u(self, onress: str, deph_u: float) -> None:
        """The upper bound of the dephasing (distance from on resonance). If the transition if on resonace, the dephasing is overwritten to equal 0.

        Args:
            onress (str): The option of whether the transition is fixed to the atomic transition.
            deph_u (float): The upper distance from the resoant value, in MHz.
        """
        if (onress == 'Y'):
            self.dephase_upper = 0.0
        else:
            if deph_u != 'Unset':
                self.dephase_upper = float(deph_u*(10**(6)))## MHz
            else:
                myinput = float(input('please enter upper dephasing value in MHz:'))
                self.dephase_upper = myinput*(10**(6))
    
    def set_laser_line(self, feild: str, laser_line: float) -> None:
        """Sets the linewidth of the laser, this only impacts the linewidth if the field value is set to 'L'.

        Args:
            feild (str): The type of field, 'L'
            laser_line (float): The linewidth of the laser in kHz.
        """
        if (feild == 'L'):
            if laser_line != 'Unset':
                self.laser_line = float((laser_line*10**(3)))#kHz
            else: 
                my_input=float(input('Please enter laser linewidth in kHz: '))
                self.laser_line  = (my_input*10**(3))#kHz
        else:
            self.laser_line = 0.0
    
    def set_beam_prop(self, feild: str, beam_prop: int) -> None:
        """Sets the propagation direction of the laser beam. This only impacts the field if the field is an 'L' field.

        Args:
            feild (str): The type of field, 'L'
            beam_prop (int): The propogation direction along the optical axis. If 1, it propgates in the positive direction, and if -1, the negative direction.
        """
        OAM_list = ['1', '-1']
        if (feild == 'L'):
            while True:
                if beam_prop !='Unset':
                    self.beam_prop = int(beam_prop)
                    break
                else:
                    my_input = input(r'if beam travels in the postive direction please input 1 if beam travels in counter direction please put -1 ' )
                    try:
                        if (my_input not in OAM_list):
                            print('popergation direction must be enter value -1 or 1')
                            continue
                        else:
                            self.beam_prop = int(my_input)
                            break
                    except:
                        print("That makes no sense...")
                        continue
        else:
            self.beam_prop = 0
    
    def __str__(self):

        print_string = str(self.q) + str(self.p) + str(self.w) +  str(self.feild) +str(self.p_rf) + str(self.ress)+str(self.scan)+str(self.deph_l)+str(self.deph_u)+ str(self.laser_line) + str(self.temp)

        return(print_string)
    
   
    def __init__(self, f = 'Unset', vm ='Unset', laser_p ='Unset', waist = 'Unset', q_set= 'Unset', ress= 'Unset', scan = 'Unset', deph_l ='Unset', deph_u = 'Unset', laser_line='unset', temp = 'Unset',beam_prop = 'Unset'):
        """To create a field parameter object, all of the key experimental parameters are used here. These will also set the numerical 'experiment' parameters, e.g. which wavelengths are being scanned over.

        Args:
            f (str, optional): The type of field 'L' for a laser, 'R' for a radio wave. Defaults to 'Unset'.
            vm (float, optional): Field strength of an RF wave, in V m^(-1). Defaults to 'Unset'.
            laser_p (float, optional): The power of the laser beam in mW. Defaults to 'Unset'.
            waist (float, optional): The width of a laser beam in um. Defaults to 'Unset'.
            q_set (int, optional): The polarisation state of the field, from the set {-1, 0, 1}. Defaults to 'Unset'.
            ress (str, optional): Whether the interacting field is on ('Y') or off ('N') resonance with the associated atomic transtion. Defaults to 'Unset'.
            scan (str, optional): Whether the field will be scanned over the atomic transition. Defaults to 'Unset'.
            deph_l (float, optional): The lower limit of the scan. Defaults to 'Unset'.
            deph_u (float, optional): The upper limit of the scan. Defaults to 'Unset'.
            laser_line (float, optional): The linewidth of the laser in kHz. Defaults to 'unset'.
            temp (float, optional): The temperature of the system in K. Defaults to 'Unset'.
            beam_prop (int, optional): The propagation direction. Defaults to 'Unset'.
        """
        self.feild = "unset"
        self.onress = "unset"
        self.dephase_lower = np.nan
        self.dephase_upper = np.nan
        self.dephase_step =np.nan
        self.p_rf = np.nan
        self.q = np.nan
        self.p = np.nan
        self.w = np.nan
        self.scan = "unset"
        self.laser_line = np.nan
        self.temp = np.nan
        self.beam_prop= np.nan

        self.set_feild(f)
        self.set_MW_feild(self.feild, vm)
        self.set_power(self.feild, laser_p)
        self.set_waist(self.feild, waist)
        self.set_q(q_set)
        self.set_onress(ress)
        self.set_dphase_l(self.onress, deph_l)
        self.set_dphase_u(self.onress, deph_u)
        self.set_scanning_source(scan, self.onress)
        self.set_laser_line(self.feild, laser_line)
        self.set_temp(temp)
        self.set_beam_prop(self.feild, beam_prop)       

class state_describing_object:
    """ The is the Qunautm numbers assiociated with the N level system. 
    """
    def set_N(self, N: int) -> None:   
        """Sets tthe numerical value for the principal quantum number

        Args:
            N (int): The principal quantum number
        """
        while True:
             
            if (N != 'Unset'):
                self.energy_level = int(N)
                break
            else:      
                my_input = input("Atomic Energy Level (N) [Integer value]:")

                try:
                    
                    my_input = int(my_input)

                    if (type(my_input) != type(1) or my_input <=0 ):
                        print("That makes no sense...")
                        continue
                    else: 
                        self.energy_level = my_input
                        break
                except: 
                    print("That makes no sense...")
                    continue

    def get_L(self, l: str):
        """Sets the string representation of the orbital angular momentum symbol.

        Args:
            l (str): The string representing the orbital angular momentum.
        """
        
        OAM_list = ['S','P','D','F','G']

        while True: 
            if l !='Unset':
                self.L = str.upper(l)
                break
            else:
                my_input = input("What is L [String value]:")

                try:
                    
                    my_input = str.upper(my_input)

                    if (my_input not in OAM_list):
                        print("That makes no sense... it is not in " +
                              OAM_list)
                        continue
                    else: 
                        self.L = my_input
                    break
                except: 
                    print("That makes no sense...")
                    continue

    def get_num_L(self, L: str) -> None:
        """Sets the numerical value for the orbital angular momentum.

        Args:
            L (str): The symbol, as part of a term symbol, for the orbital angular momentum, e.g. 'S', 'P', etc.
        """
        if L =="S":
            num_L = 0
        elif L == "P":
            num_L= 1
        elif L == "D":
            num_L= 2
        elif L == "F":
            num_L= 3
        elif L == "G":
            num_L= 4
        else: 
            num_L = np.nan
            print("Invalid term symbol")
        self.numerical_L = num_L

    def get_J(self, L: int, j: str) -> None:
        """Sets the total angular momentum to a value of L+0.5 or L-0.5.

        Args:
            L (int): The integer value of the orbital angular momentum (the quantum number).
            j (str): A string denoting whether the total angular momentum is positive or negative. This relates to the spin of the promoted electron, i.e. is it co- or counter-propogating with respect to the orbit.
        """
        while True: 
            if j != 'Unset':
                j = str.upper(j)
                if j == "Y": 
                    self.J = L+0.5
                    break
                else: 
                    self.J = L-0.5
                    break        
            else:
                my_input = input("High Total Angular Momentum [Y/N]?:")

                try:
                    
                    my_input = str.upper(my_input)

                    if ((my_input != "Y" and my_input != "N") or
                        (my_input == "N" and L == 0)):
                        print(L)
                        print(my_input)
                        print("That makes no sense...")
                        continue
                    else: 
                        if my_input == "Y": 
                            self.J = L+0.5
                            break
                        else: 
                            self.J = L-0.5
                            break
                except: 
                    print("That makes no sense (crashed)...")
                    continue 
    
    def get_mj(self, J: float, MJ: Union[float, str]) -> None:
        """Sets the mj value. mj = -j, -j+1 ...j-1, j. 

        Args:
            J (float): Total angualer momentum
            MJ (Union[float, str]): The second angular momentum number.
        """
        while True:
            if MJ != 'Unset':
                self.mj = float(MJ)
                break
            else:
                myinput = input(
                    'enter value for Secondary Angular Momentum, mj:')

                try:
                    myinput = float(myinput)

                    if ((myinput > J) or (myinput < -J) or
                        ((myinput % 0.5) != 0)):
                        print('that makes no sense. mj = -j, -j+1 ...j-1, j:')
                        print(J)
                    else:
                        self.mj = myinput
                        break   
                except ValueError: 
                    print("That makes no sense (crashed)...")
                    continue  

    def __str__(self) -> str:

        print_string = str(self.energy_level) + str(self.L) + str(self.numerical_L) +  str(self.J) +str(self.mj)

        return print_string
    
    def get_name(self, n: int, l: str, j: float) -> None:
        """Sets the name property of the state describing object. This is a representation of the term symbol, i.e. 5S_{0.5}.

        Args:
            n (int): The principal quantum number.
            l (str): The orbital angular momentum symbol.
            j (float): The total angular momentum.
        """

        self.name =str(n)+str(l)+str('_')+str('{')+str(j)+str('}')
        

    def __init__(self, n = 'Unset' , l = 'Unset', j ='Unset', MJ = 'Unset') -> None:
        """When creating a state describing object, any 'Unset' parameters will be requested from the user via the commmand line. This allows for its use in both command line applications, and as part of larger simulations.

        Args:
            n (int, optional): The principal quantum number. Defaults to 'Unset'.
            l (str, optional): A single character that describes the orbital angualar momentum, which will be: 'S' if L=0, 'P' if L=1 etc. This is not case sensitive. Defaults to 'Unset'.
            j (str, optional): This ask if the total angular moment is positive giving the higher value j = l+1/2. Chose from 'Y' is postive j = l+1/2 or 'N' for negative j = l-1/2. Defaults to 'Unset'.
            MJ (float, optional): This is the secondary angular momentum. This should be mj = -j, -j+1 ...j-1, j. Defaults to 'Unset'.
        """
        self.energy_level = "Unset"
        self.L = "unset"
        self.numerical_L = np.nan
        self.J = np.nan
        self.mj = np.nan
        self.name = "unset"
        
        self.set_N(n)
        self.get_L(l)
        self.get_num_L(self.L)
        self.get_J(self.numerical_L,j)
        self.get_mj(self.J, MJ)
        self.get_name(self.energy_level,self.L, self.J)
        pass 
    
def create_subplot(x, y, title, ax=None, color=None, string_data=None, labels=None, linestyles=None):
    if ax is None:
        ax = plt.axes()    

    ax.grid(color='black', linestyle='--', alpha=0.5)
    ax.grid(which='minor', color='grey', linestyle='--', alpha=0.5)
    ax.set_title(title)
       
    if labels is not None and string_data is not None:
        label_update = labels.format(*string_data)
        plots = ax.plot(x, y, color=color, label=label_update, linestyle=linestyles)
        ax.legend(loc='upper right', ncols=1, prop={'size': 18})
    else:
        plots = ax.plot(x, y, color=color)
    
    ax.set_xlim(-25, 25)
    
    return plots

def create_inset(ax, x, y, color=None, linestyles=None, x1=None, x2=None, y1=None, y2=None, Z=None):
    if x1 is not None:
        #axins = zoomed_inset_axes(ax, zoom=Z, loc=2)
        axins = inset_axes(ax, width="20%", height="80%",  loc='upper left', bbox_to_anchor=(0.04, 0, 1, 0.95), bbox_transform= ax.transAxes)#, loc='upper left' 

        for i in range(len(y)):
            axins.plot(x[i], y[i], color=color[i], linestyle=linestyles[i])
        
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        
        axins.yaxis.get_major_locator().set_params(nbins=3)
        axins.xaxis.get_major_locator().set_params(nbins=3)
        #axins.tick_params(labelleft=False, labelbottom=False)
        axins.grid(color='black', linestyle='--', alpha=0.5)
        axins.grid(which='minor', color='grey', linestyle='--', alpha=0.5)
    
        mark_inset(ax, axins, loc1=3, loc2=1, fc="none", ec="black", color = 'black')


if __name__ =='__main__':
   
    
   ########## data for main program theory#######
   
    config = liams_windows_config.liams_windows_config() # this is my windows save location 
    
    atom = [Rb85() ] # this is the atom of interest Rb85(), Cs() 

    #this asks the user how may state of the atom they wish to search
    natline = [[],[],[],[],[]]
    EITcon = [[],[],[],[],[]]
    laserlin =[]
     
    # ax1 ,ax2,ax1 ax2,
    fig, ( ax3) = plt.subplots(1,1, sharex =True, figsize= (14,12))
    
    
    title =['A) Probe', 'A) EIT Contrast','A) RF Sensing']
    Supertitle = '$^{133} Cs$ Absorption Spectra Changes for off Resonance RF'
    
    fig.supylabel('Transmission')
    fig.supxlabel('$\Delta$p (MHz)')
    fig.suptitle(Supertitle)
    
    test_chi_f_list_nat = [[],[],[],[],[]]
    test_chi_f_list_nat_lin = [[],[],[],[],[]]    
    delta_p =[[],[],[],[],[]]
    trans_nat=[[],[],[],[],[]]


    Linewidth = [0.1, 1, 10, 100, 1000] # laser linewidths
    Rf_ress= [-20, -10, 0, 10, 20] # off resonances for RF
    couplingpower = [32, 64, 96, 128, 160 ]  #Cs coupling [17.3, 50, 100, 150, 200] #Rb coupling [32, 64, 96, 128, 160 ]
    rfpower = [0.5, 1, 1.5, 2, 3.57]#rb RF poweer[0.5, 1, 1.5, 2, 3.57] #Cs RF power [0.05, 0.25, 0.5, 0.75, 1] 
    probe_power = Rf_ress #rb probe powers [ 0.005, 0.01, 0.017, 0.025, 0.05 ]# Cs probe powers [0.0032, 0.0064, 0.0096,0.0128, 0.016]
    
    ###### variables ######    
    blank = [] # Full half maxium store for probe only transition
    EITconl = []# Full half maxium store for probe only transition
    rfspilt =[]# RF splitting measurement
    
    for p in range(len(probe_power)): # this loops through the probe powers
        print(probe_power[p])
        for number_states in range (2, 5): ## this will display rf senisng when (2, 5)
        #number_states = 2 #int(input('please enter the number of state you are interested in: '))
            print(number_states)
            
            waist_p =160 #hree L2000  # Cs paper waist 708, Rb paper wait 160
            #print(waist_p )
            waist_c =244 #three L 2000    #Cs paper waist 708, Rb paper waist 244
        
            
            chi_list_p =[]
            chi_list_p_basic= []
            chi_list_c = []
            chi_list_rf =[]
            
            trans_doppler =[]
            
            chi_max =[] 
            chi_p = [] 
            
            a = [] #used to store the probe power as a list
        
            temp = 298.15 #20.85 #298.15 # 25 C
    
    
            for i in range(0,   1):
                
                state_array = [] # this is a list generated by the user for state information 
                laser_trans = [] # this is a list to store laser information
                #a.append(probe_p)
               # print('temp= '+ str(temp[i]))
               
                for j in range(number_states):
                    #this is used to build the state information array from user input 
                    #print('please enter information for ' + str(i) + 'th tranistion-----')
                    if j == 0:
                        this_state = state_describing_object(5 , 's',  'y', 0.5) #this is the object built from user input
                    elif j == 1:
                        this_state = state_describing_object(5 , 'p', 'y', 1.5) #this is the object built from user input
                    elif j == 2:
                        this_state = state_describing_object(26, 'd', 'Y', 2.5)#83, 'D',  'y', 2.5) #this is the object built from user input
                    elif j == 3:
                        this_state = state_describing_object(27, 'p', 'y', 1.5)
                    #else:
                    #    this_state = state_describing_object(44,'g','y', 4.5)# 84, 'd',  'y', 2.5) #this is the object built from user input    
                                                
                    state_array.append(this_state) # this is n number of state informatin the user asked for
                
                for b in range(number_states-1):
                    #print('please enter laser information for ' + str(i) + 'th tranistion-----')
                    if b == 0:  # 0.000589, 0.00101, 0.001313, 0.002095,0.002589, 0.00334, 0.00354
                        laser_param = laser_parameter('L','Unset',0.0175 , waist_p, 1,'N','Y',-250, 250,  1, temp, 1) #112 Isat = 0.0005208probe beam has 100 micron beam waist
                    elif b == 1: #44.15
                        laser_param = laser_parameter('L','Unset', 32, waist_c, 1,'Y', 'n','Unset','Unset', 1, temp, 1) #0.00512 3 laser system 
                    
                    elif b == 2: #65
                        laser_param = laser_parameter('r', 0.5, 'Unset', 'Unset', -1,'N','N',probe_power[p],probe_power[p], 100, temp, -1)    
                    #else:
                    #    laser_param= laser_parameter('R',5.9, 'Unset', 'Unset', 1, 'Y', 'N', 'Unset','Unset', 'Unset', temp)
                    #laser_parameter()
                    laser_trans.append(laser_param)  
                
                dephasing_step = 9003#1000001 #int(input('please enter number of steps in the dephasing:'))
                #print(laser_trans)   
                
        
                # This will search and merge the state of interest with the data frame of state generated in CSV 
                # Rb_state_lifetimes_test.csv (this file is generated using the Rb_lifetimes.py program) 
                intrested_states = transistion.transition_picker(state_array, atom)
                
                config = liams_windows_config.liams_windows_config() # This is the config route and file location on my work laptop
        
                intrested_states.to_csv(config.project_loc+"/Data/CSV/extracted.CSV") # test CSV for pulled states
                
                
        
                Trans_wave = atom_wave.transition_generator_excitation(atom, state_array, laser_trans) # This calls the function to get the tranistion rabi frequency 
                decay_wave = atom_wave.transition_generator_decay(atom, state_array, laser_trans)# this calls the function to get the decay rates
        
                
                Trans_wave.to_csv(config.project_loc+"/Data/CSV/trans_wave.CSV") # test CSV for pulled states
                decay_wave.to_csv(config.project_loc+"/Data/CSV/decay_wave.CSV") # test CSV for pulled states
        
        
                ### steady state builder for EIT for an nxn steady state system
                chi, delta_p_1 = EITSS.EIT_builder(number_states,Trans_wave, decay_wave, laser_trans, dephasing_step)
            
                
                
                
                chi_list_p_basic.append(chi[:,0]) ## non doppler profile 
                delta_p[p].append(delta_p_1)
                
                
                #### Calculate the real Transmission Units chek
                ####constanst values
                
                e_0 = np.longdouble(8.85418782*(10**(-12)))
                Bolt = np.longdouble(1.3806*10**(-23))## m^2 kg s^-2 K^-1
                a0 = np.longdouble(5.29*10**(-11)) # bohr radius m
                e = np.longdouble(1.602*(10**(-19)))# charge of an electron C
                wave_num = np.longdouble(2*np.pi/(780*10**(-9)))
                
                T = laser_param.temp #
                pres2 = np.longdouble(atom[0].getPressure(T)) #Pa
                n2 = np.longdouble(atom[0].getNumberDensity(temp))*((np.pi*((waist_p*10**(-6))**2)*0.072)/(np.pi*((0.0125)**2)*0.072))
                
                
        
                
                rabi_p = np.longdouble((Trans_wave.iloc[0]['rabi Freq rads/s']))
                
                dipole_12 = np.longdouble(Trans_wave.iloc[0]['Dipole_moment']*(a0*e))
                
                E = np.longdouble((hbar*rabi_p)/dipole_12)
                
                L = 0.0718
        
                Alpha = ((2*dipole_12)/(E*e_0))
                
                ## non doppler Alpha*n* Alpha*n2*
                Chi_val_nat = Alpha*n2*(np.imag(chi_list_p_basic[0]))
                
                test_chi_f_list_nat[p].append(Chi_val_nat)
               ### peak and splitting for spectra. 
                if number_states == 2:            
                    peaks, _ = sn.find_peaks(Chi_val_nat) #Chi_val_nat peak finder
                    
                    result_half_n, _, _,_ = sn.peak_widths((Chi_val_nat), peaks, rel_height = 0.5)
                    
                    
                    coarse_upper =(laser_trans[0].dephase_upper - 25*10**6)/(dephasing_step/3)
                    coarse_lower = (-25*10**6-laser_trans[0].dephase_lower)/(dephasing_step/3)
                    high_res_mid =(25*10**6 + 25*10**6)/(dephasing_step/3)
                    
                    
                    natline[p].extend(result_half_n*-1*(coarse_upper-high_res_mid-coarse_lower)/(10**6))
                    blank.append(natline[p][0]/(2*np.pi))
                    
                if number_states == 3:            
                    peaks, _ = sn.find_peaks(test_chi_f_list_nat[p][0]-test_chi_f_list_nat[p][1]) #Chi_val_nat
                    result_half_n_EIT, _, _,_ = sn.peak_widths((test_chi_f_list_nat[p][0]-test_chi_f_list_nat[p][1]), peaks, rel_height = 0.5)
                                   
                    coarse_upper =(laser_trans[0].dephase_upper - 25*10**6)/(dephasing_step/3)
                    coarse_lower = (-25*10**6-laser_trans[0].dephase_lower)/(dephasing_step/3)
                    high_res_mid =(25*10**6 + 25*10**6)/(dephasing_step/3)
                    
                    
                    EITcon[p].extend(result_half_n_EIT*-1*(coarse_upper-high_res_mid-coarse_lower)/(10**6))
                    
                    EITconl.append(EITcon[p][0]/(2*np.pi))
                    
                if number_states == 4:            
                    peaks, _ = sn.find_peaks(test_chi_f_list_nat[p][0]-test_chi_f_list_nat[p][2]) #Chi_val_nat
                    if len(peaks) > 1:
                        result_half_n_rf, _, _,_ = sn.peak_widths((test_chi_f_list_nat[p][0]-test_chi_f_list_nat[p][2]), peaks, rel_height = 0.5)
                    
                        peak_0 = int(peaks[0])
                        peak_1 = int(peaks[1])
                        peaks_sp = peak_1-peak_0
                    
                        coarse_upper =(laser_trans[0].dephase_upper - 25*10**6)/(dephasing_step/3)
                        coarse_lower = (-25*10**6-laser_trans[0].dephase_lower)/(dephasing_step/3)
                        high_res_mid =(25*10**6 + 25*10**6)/(dephasing_step/3)
                    
                    
                        x = ((peaks_sp)*-1*(coarse_upper-high_res_mid-coarse_lower)/(10**6))
                    
                        rfspilt.append(x)
                    else:
                        rfspilt.append(0)
                    
                    

        
                exp_dat_1 = (-wave_num*L*(Chi_val_nat))
                trans_nat[p].append(np.exp(exp_dat_1)) ### Transmission data 
       
    #Functionalise plot data for multiple runs
    color = ['black', 'blue', 'red', 'purple', 'green'] #list colors
    linestyles  =['solid', 'dashed','dotted', 'dashdot', (0,(1,1))]
    labelone = r'2$\pi\cdot${1:.3f} MHz' #2$\pi$ $\mu$W,
    labeltwo = r'2$\pi\cdot${1:.3f} MHz' #2$\pi$ {0:.2f} mW,{0:.2f} mW,
    labelthree = r'{0:.0f} MHz, {1:.3f} MHz' #{0:.2f} V/m,
    
    # Loop over probe_power to create plots
    for p in range(len(probe_power)):
        #create_subplot(delta_p[p][0]/(10**6), trans_nat[p][0], title[0], ax1, color[p], (probe_power[p], blank[p],), labelone, linestyles[p])
        #create_subplot(delta_p[p][1]/(10**6), trans_nat[p][1], title[1], ax2, color[p], (probe_power[p], EITconl[p],), labeltwo, linestyles[p])
        create_subplot(delta_p[p][2]/(10**6), trans_nat[p][2], title[2], ax3, color[p], (probe_power[p], rfspilt[p],), labelthree, linestyles[p])
    
    #Create inset
    #create_inset(ax1, [delta_p[i][0]/(10**6) for i in range(len(probe_power))], [trans_nat[i][0] for i in range(len(probe_power))], color, linestyles, -1, 1, 0.505, 0.525, 10)

    #create_inset(ax2, [delta_p[i][1]/(10**6) for i in range(len(probe_power))], [trans_nat[i][1] for i in range(len(probe_power))], color, linestyles, -4, 4, 0.5, 0.9, 1.1)

    #create_inset(ax3, [delta_p[i][2]/(10**6) for i in range(len(probe_power))], [trans_nat[i][2] for i in range(len(probe_power))], color, linestyles, -6, -3, 0.5, 0.9, 1.1)
    
    plt.show()    
    plt.savefig(config.project_loc + '/Data/Graphs/'+ Supertitle +'.png')
    plt.close()
    print('computation complete')
    pass