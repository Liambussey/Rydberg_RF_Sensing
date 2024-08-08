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


#from Rydberg_Cs_example_paper_natural_decay_multiple_lines import state_describing_object, laser_parameter

def state_array_dataframe_generator(state_array: list["state_describing_object"], laser: list["laser_parameter"]) -> pd.DataFrame:
    """This collates the information require for excitation tranisition between neighbouring states in the State_array, these transition must obey Hunds selection rules.

    Args:
        state_array (list[state_describing_object]): An object that will contain the specific atomic infomation.
        laser (list[laser_parameter]): A list containing the interation feilds' infomation.

    Returns:
        pd.DataFrame: This return the dataframe contain the atomic excitation transition data only neighbouring states in the state_array will have transitions occuring. 
    """
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
    return state_df

def state_array_dataframe_generator_decay(state_array: list["state_describing_object"], laser: list["laser_parameter"]) -> pd.DataFrame:
    """This collates the information require for decay tranisition between neighbouring states in the State_array, these transition must obey Hunds selection rules.
        To create this decay dataframe the state and laser arrays need to stay in the same order.

    Args:
        state_array (list[state_describing_object]): An object that will contain the specific atomic infomation.
        laser (list[laser_parameter]): A list containing the interation feilds' infomation.

    Returns:
        pd.DataFrame: This return the dataframe contain the atomic decay transition data only neighbouring states in the state_array will have transitions occuring. 
    """
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

def getwavelength(atom: list[arc.Rubidium85], n1: int, l1: int, j1: float, n2: int, l2: int, j2: float) -> list[float]:
    """Gets the wavelength for the defined tranision using ARC `getTransitionWavelength` 

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        n1 (int): Lower principal quantum number.
        l1 (int): lower orbital angular momentum quantum number.
        j1 (float): lower total angular momentum.
        n2 (int): Upper principal quantum number.
        l2 (int): Upper orbital angular momentum quantum number.
        j2 (float): Upper total angular momentum.

    Returns:
        list[float]: return the wavelength assiociate with the transition (nm).
    """
    
    return [(atom[i].getTransitionWavelength(n1,l1,j1,n2,l2,j2)*1e9) for i in range(len(atom))]
    

def getfrequency(atom: list[arc.Rubidium85], n1: int, l1: int, j1: float, n2: int, l2: int ,j2: float) -> list[float]:
    """Gets the Frequency for the defined tranision using ARC `getTransitionFrequency` 

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        n1 (int): Lower principal quantum number.
        l1 (int): lower orbital angular momentum quantum number.
        j1 (float): lower total angular momentum.
        n2 (int): Upper principal quantum number.
        l2 (int): Upper orbital angular momentum quantum number.
        j2 (float): Upper total angular momentum.

    Returns:
        list[float]: return the Frequency assiociate with the transition (Hz).
    """
    
    return [(atom[i].getTransitionFrequency(n1,l1,j1,n2,l2,j2)) for i in range(len(atom))]
    
def getdecay(atom: list[arc.Rubidium85], n1: int, l1: int, j1: float, n2: int, l2: int, j2: float, T: float) -> list[float]:
    """this finds the spontanous decay from state 1 to state 2 using ARC `getTransitionRate`

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        n1 (int): Lower principal quantum number.
        l1 (int): lower orbital angular momentum quantum number.
        j1 (float): lower total angular momentum.
        n2 (int): Upper principal quantum number.
        l2 (int): Upper orbital angular momentum quantum number.
        j2 (float): Upper total angular momentum.
        T (float): Temperature in K.

    Returns:
        list[float]: Transition rates in Hz
    """
    return [(atom[i].getTransitionRate(n1, l1, j1, n2, l2, j2, temperature=T)) for i in range(len(atom))]
 
    

def getmass(atom: list[arc.Rubidium85]) -> list[float]:
    """Gets the average mass of the atomic spiecies using ARC `mass` 

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.

    Returns:
        list[float]: Return the mass of the atom in kg
    """
    return [atom[i].mass for i in range(len(atom))]
    

def getdens(atom: list[arc.Rubidium85], temp: float) -> list[float]:
    """Atom number density at given temperature

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        temp (float): Temperature in K.

    Returns:
        list[float]: Returns the number desnity of atoms given a particular temperature in m^(-3)
    """
    return [atom[i].getNumberDensity(temp) for i in range(len(atom))]
    

def transition_generator_excitation(atom_list: list[arc.Rubidium85], interested_state: pd.DataFrame, laser: list["laser_parameter"])-> pd.DataFrame:
    """This calculate the physical parameter which impact the excitation pathways in the defined atomic system

    Args:
        atom_list (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        interested_state (pd.DataFrame): This is a Dataframe of atomic state the user has entered and find the state lifetime assoicated with those states.  
        laser (list[laser_parameter]): A list containing the interation feilds' infomation.

    Returns:
        pd.DataFrame: This returns a Dataframe with the expanded information on the excitation pathway gathered using ARC
    """
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

def transition_generator_decay(atom_list: list[arc.Rubidium85], interested_state: pd.DataFrame, laser: list["laser_parameter"]) -> pd.DataFrame:
    """This calculate the physical parameter which impact the decay pathways in the defined atomic system

    Args:
        atom_list (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        interested_state (pd.DataFrame): This is a Dataframe of atomic state the user has entered and find the state lifetime assoicated with those states.  
        laser (list[laser_parameter]): A list containing the interation feilds' infomation.

    Returns:
        pd.DataFrame: This returns a Dataframe with the expanded information on the decay pathway gathered using ARC
    """
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
                                    
def get_rabi(atom_list: list[arc.Rubidium85], feild: str, n1: int, l1: int, j1: float, mj1: float, n2: int, l2: int, j2: float, q: int, p: float, w: float, p_rf: float, linewidth: float) ->list[float]:
    """Gets the Rabi Frequency for the defined tranision using ARC `getRabiFrequency` or `getRabiFrequency2`

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        feild (str):  If the field is a laser field 'L' this will use `getRabiFrequency` else it will use `getRabiFrequency2` for Radio frequency 'R'
        n1 (int): Lower principal quantum number.
        l1 (int): lower orbital angular momentum quantum number.
        j1 (float): lower total angular momentum.
        mj1 (float): lower second angular momentum number.
        n2 (int): Upper principal quantum number.
        l2 (int): Upper orbital angular momentum quantum number.
        j2 (float): Upper total angular momentum.
        q (int): The polarisation state of the field, from the set {-1, 0, 1}.
        p (float): The power of the laser beam in mW.
        w (float): The width of a laser beam in um.
        p_rf (float):Field strength of an RF wave, in V m^(-1).
        linewidth (float): laser linewidth. not used in this function 

    Returns:
        list[float]: Returns the rabi frequency for a given transition in rad s^(-1)
    """
    if (feild == 'L'):
    
        return [(atom_list[i].getRabiFrequency(n1, l1, j1, mj1,
                                                    n2, l2, j2,
                                                    q, p, w)) for i in range(len(atom_list))]
        
        
    else:

        return [(atom_list[i].getRabiFrequency2(n1, l1, j1, mj1,
                                                    n2, l2, j2,q,
                                                    p_rf)) for i in range(len(atom_list))]

    
def get_dipole(atom_list, feild: str ,n1: int, l1: int, j1: float, mj1: float, n2: int, l2: int, j2: float, mj2: float, q: int, p: float, w: float, p_rf: float) -> list[float]:
    """Gets the Dipole matrix element for the defined tranision using ARC `getDipoleMatrixElement`

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        feild (str):  If the field is a laser field 'L' this will use `getRabiFrequency` else it will use `getRabiFrequency2` for Radio frequency 'R'
        n1 (int): Lower principal quantum number.
        l1 (int): lower orbital angular momentum quantum number.
        j1 (float): lower total angular momentum.
        mj1 (float): lower second angular momentum number.
        n2 (int): Upper principal quantum number.
        l2 (int): Upper orbital angular momentum quantum number.
        j2 (float): Upper total angular momentum.
        mj2 (float): upper second angular momentum number. 
        q (int): The polarisation state of the field, from the set {-1, 0, 1}.
        p (float): The power of the laser beam in mW.
        w (float): The width of a laser beam in um.
        p_rf (float):Field strength of an RF wave, in V m^(-1).


    Returns:
        list[float]: This calculated the Dipole matrix element units (a_0 e)
    """
    return [(atom_list[i].getDipoleMatrixElement(n1, l1, j1, mj1,
                                                         n2, l2, j2, mj2, q)) 
                                                        for i in range(len(atom_list))]

def get_speed(atom_list: list[arc.Rubidium85], feild: str ,n1: int, l1: int, j1: float, mj1: float, n2: int, l2: int, j2: float, mj2: float, q: int, p: float, w: float, p_rf: float, temp: float) -> list[float]:
    """Gets the average speed of the atoms using ARC `getAverageSpeed`

    Args:
        atom (list[arc.Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        feild (str):  If the field is a laser field 'L' this will use `getRabiFrequency` else it will use `getRabiFrequency2` for Radio frequency 'R'. not used in this function 
        n1 (int): Lower principal quantum number. not used in this function 
        l1 (int): lower orbital angular momentum quantum number. not used in this function 
        j1 (float): lower total angular momentum. not used in this function 
        mj1 (float): lower second angular momentum number. not used in this function 
        n2 (int): Upper principal quantum number. not used in this function not used in this function 
        l2 (int): Upper orbital angular momentum quantum number.not used in this function 
        j2 (float): Upper total angular momentum.not used in this function 
        mj2 (float): upper second angular momentum number. not used in this function 
        q (int): The polarisation state of the field, from the set {-1, 0, 1}.not used in this function 
        p (float): The power of the laser beam in mW.not used in this function 
        w (float): The width of a laser beam in um.not used in this function 
        p_rf (float):Field strength of an RF wave, in V m^(-1).not used in this function  
        temp (float): Temperature in K

    Returns:
        list[float]: This return the average velocity of the atoms
    """
    return [atom_list[i].getAverageSpeed(temp) #returns ms-1
                 for i in range(len(atom_list))]
    

if __name__ =='__main__':
    print('this should be inported into main file')
    
    pass 

