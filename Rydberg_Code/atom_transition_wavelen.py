# -*- coding: utf-8 -*-
"""

 test file wavelength search
     
 
"""
import sys
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory') 
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory/src')
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/Github/rydberg_theory/src/python')

import pandas as pd 
import numpy as np
import liams_windows_config

import arc

def state_array_dataframe_generator(state_array, laser):
  #this created the data frame for user input data from the main program which is then ready to be amended 
    state_df = pd.DataFrame(columns = ['transition','feild','n1','l1','j1',
                                       'mj1','n2','l2','j2','mj2','q','p','w',
                                       'p_rf', 'linewidth'])
    #print(len(state_array)-1)
    for i in range(len(state_array)-1):
        single_record = pd.DataFrame({"transition": state_array[i].name +str(' -> ')+state_array[i+1].name,
                         "n1": [state_array[i].energy_level],
                         "l1": [state_array[i].numerical_L],
                         "j1": [state_array[i].J],
                         "mj1": [state_array[i].mj],
                         "n2": [state_array[i+1].energy_level],
                         "l2": [state_array[i+1].numerical_L],
                         "j2": [state_array[i+1].J],
                         "mj2": [state_array[i+1].mj],
                         "feild": [laser[i].feild],
                         "q": [laser[i].q],
                         "p": [laser[i].p],
                         "w": [laser[i].w],
                         "p_rf": [laser[i].p_rf],
                         'temp': [laser[i].temp],
                         'linewidth': [laser[i].laser_line]})
        state_df = pd.concat([state_df, single_record])
    return(state_df)

def state_array_dataframe_generator_decay(state_array, laser):
    #this created the data frame for user input data from the main program which is then ready to be amended 
    state_df = pd.DataFrame(columns = ['transition','n2','l2','j2',
                                       'mj2','n1','l1','j1','mj1'])

    for i in range(len(state_array)-1):
        single_record = pd.DataFrame({"transition": state_array[i+1].name +str(' -> ')+state_array[i].name,
                         "n2": [state_array[i].energy_level],
                         "l2": [state_array[i].numerical_L],
                         "j2": [state_array[i].J],
                         "mj2": [state_array[i].mj],
                         "n1": [state_array[i+1].energy_level],
                         "l1": [state_array[i+1].numerical_L],
                         "j1": [state_array[i+1].J],
                         "mj1": [state_array[i+1].mj],
                         'temp': [laser[i].temp]})
        state_df = pd.concat([state_df, single_record])
    return(state_df)

def getwavelength(atom,n1,l1,j1,n2,l2,j2):
    # this returns the wavelength associated with the transition
    val_list =[(atom[i].getTransitionWavelength(n1,l1,j1,n2,l2,j2)*1e9) for i in range(len(atom))]
    return(val_list)

def getfrequency(atom,n1,l1,j1,n2,l2,j2):
    # this returns the frequency associated with the transition
    val_list =[(atom[i].getTransitionFrequency(n1,l1,j1,n2,l2,j2)) for i in range(len(atom))]
    return(val_list)

def getdecay(atom, n1, l1, j1, n2, l2, j2, T):
    ### this finds the spontanous decay from state 1 to state 2
    val_list = [(atom[i].getTransitionRate(n1, l1, j1, n2, l2, j2, temperature=T)) for i in range(len(atom))]
    ### return tranisition rante in s-1 (Hz)
    return(val_list)

def getmass(atom):
    val_list =[atom[i].mass for i in range(len(atom))]
    return(val_list)

def getdens(atom, temp):
    val_list =[atom[i].getNumberDensity(temp) for i in range(len(atom))]
    return(val_list)

def transition_generator_excitation(atom_list, interested_state,laser):
    #atom = atomic species and natrually occuring or isotopic
    #interested states provided the user defined states
    #This will find the relevant tranistion wavelength between the given states
    #that are of interest.
    # i.e Rb 5S_{1/2} -> 5P_{3/2} == 780nm
    state_df = state_array_dataframe_generator(interested_state, laser)
    name_list = ['wavelength in nm for Rb']
    name_list2 = ['freq in Hz for Rb']
    name_list3 = ['Transition Rate(Hz) with BBR for Rb']
    name_list4 = ['Transition Rate(Hz) at 0.1K for Rb']
    name_list5= ['rabi Freq rads/s']
    name_list6 = ['Dipole_moment']
    name_list7 = ['mean_speed']
    name_list8 =['mass']
    name_list9 = ['Atomic_number_dens']
    trans_wave = pd.DataFrame (columns = ['transition','feild','n1','l1','j1',
                                'mj1','n2','l2','j2','mj2', 'q','p','w','p_rf', 'temp'])

    trans_wave = pd.concat([trans_wave, state_df])

    trans_wave[name_list] = trans_wave.apply(lambda x: getwavelength(atom_list,
                                                    x.n1, x.l1, x.j1, x.n2, x.l2, x.j2),
                                                    axis=1, result_type='expand')
    trans_wave[name_list2] = trans_wave.apply(lambda x: getfrequency(atom_list,
                                                    x.n1, x.l1, x.j1, x.n2, x.l2, x.j2),
                                                    axis=1, result_type='expand')
    trans_wave[name_list3] = trans_wave.apply(lambda x: getdecay(atom_list,
                                                        x.n1, x.l1, x.j1, x.n2, x.l2, x.j2, T=x.temp),
                                                        axis=1, result_type='expand')
    trans_wave[name_list4] = trans_wave.apply(lambda x: getdecay(atom_list,
                                                        x.n1, x.l1, x.j1, x.n2, x.l2, x.j2, T=0.1),
                                                        axis=1, result_type='expand')
    trans_wave[name_list5] = trans_wave.apply(lambda x: get_rabi(atom_list, x.feild,
                                                       x.n1, x.l1, x.j1, x.mj1, x.n2, x.l2, x.j2,
                                                        x.q , x.p, x.w, x.p_rf, x.linewidth),
                                                        axis=1, result_type='expand')
    
    trans_wave[name_list6] = trans_wave.apply(lambda x: get_dipole(atom_list, x.feild,
                                                       x.n1, x.l1, x.j1, x.mj1, x.n2, x.l2, x.j2, x.mj2,
                                                        x.q , x.p, x.w, x.p_rf),
                                                        axis=1, result_type='expand')
    
    trans_wave[name_list7] = trans_wave.apply(lambda x: get_speed(atom_list, x.feild,
                                                       x.n1, x.l1, x.j1, x.mj1, x.n2, x.l2, x.j2, x.mj2,
                                                        x.q , x.p, x.w, x.p_rf, x.temp),
                                                        axis=1, result_type='expand')
    trans_wave[name_list8] = trans_wave.apply(lambda x: getmass(atom_list),
                                                        axis=1, result_type='expand')
    trans_wave[name_list9] = trans_wave.apply(lambda x: getdens(atom_list, x.temp),
                                                        axis=1, result_type='expand')
    return (trans_wave)

def transition_generator_decay(atom_list, interested_state, laser):
    #atom = atomic species and natrually occuring or isotopic
    #interested states provided the user defined states
    #This will find the relevant tranistion wavelength between the given states
    #that are of interest.
    # i.e Rb 5S_{1/2} -> 5P_{3/2} == 780nm
    state_df = state_array_dataframe_generator_decay(interested_state, laser)
    name_list = ['wavelength in nm for Rb']
    name_list2 = ['freq in Hz for Rb']
    name_list3 = ['Transition Rate (Hz) with BBR for Rb']
    name_list4 = ['Transition Rate at (Hz) 0.1K for Rb']
    
    trans_wave = pd.DataFrame (columns = ['transition','n1','l1','j1',
                                'mj1','n2','l2','j2','mj2', 'temp'])

    trans_wave = pd.concat([trans_wave, state_df])

    trans_wave[name_list] = trans_wave.apply(lambda x: getwavelength(atom_list,
                                                    x.n1, x.l1, x.j1, x.n2, x.l2, x.j2),
                                                    axis=1, result_type='expand')

    trans_wave[name_list2] = trans_wave.apply(lambda x: getfrequency(atom_list,
                                                    x.n1, x.l1, x.j1, x.n2, x.l2, x.j2),
                                                    axis=1, result_type='expand')

    trans_wave[name_list3] = trans_wave.apply(lambda x: getdecay(atom_list,
                                                        x.n1, x.l1, x.j1, x.n2, x.l2, x.j2, x.temp),
                                                        axis=1, result_type='expand')

    trans_wave[name_list4] = trans_wave.apply(lambda x: getdecay(atom_list,
                                                        x.n1, x.l1, x.j1, x.n2, x.l2, x.j2, 0.0),
                                                        axis=1, result_type='expand')
    return (trans_wave)
#### Using the FWHM to calculate laser waist (must check units as linewidth is giveing in KHz) can be changed back if needed using   w
                                    
def get_rabi(atom_list, feild ,n1, l1, j1, mj1, n2, l2, j2, q, p, w, p_rf, linewidth):
    
    if (feild == 'L'):
    
        rabiFreq = [(atom_list[i].getRabiFrequency(n1, l1, j1, mj1,
                                                    n2, l2, j2,
                                                    q, p, w)) for i in range(len(atom_list))]
        rabiFreq = rabiFreq ### return in rad-1 divide by 2Pi for Hz
        return(rabiFreq)
        
    else:

        rabiFreq = [(atom_list[i].getRabiFrequency2(n1, l1, j1, mj1,
                                                    n2, l2, j2,q,
                                                    p_rf)) for i in range(len(atom_list))]
        rabiFreq =rabiFreq  ### return in rad-1 divide by 2Pi for Hz
        return(rabiFreq)
    
def get_dipole(atom_list, feild ,n1, l1, j1, mj1, n2, l2, j2, mj2, q, p, w, p_rf):
    
    dipolemoment = [(atom_list[i].getDipoleMatrixElement(n1, l1, j1, mj1,
                                                         n2, l2, j2, mj2, q)) 
                                                        for i in range(len(atom_list))]
    return(dipolemoment)

def get_speed(atom_list, feild ,n1, l1, j1, mj1, n2, l2, j2, mj2, q, p, w, p_rf, temp):
    
    meanspeed = [atom_list[i].getAverageSpeed(temp) #returns ms-1
                 for i in range(len(atom_list))]
    return(meanspeed)

if __name__ =='__main__':
    print('this should be inported into main file')
    
    pass 

