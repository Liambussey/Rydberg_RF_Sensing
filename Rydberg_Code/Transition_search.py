# -*- coding: utf-8 -*-
"""
   '''
    @author: Liam Bussey 
    @return: Panda Dataframe, table of....
    This program searches the data frame created by Rydberg lifetime to pull two indivdual tranisition
    @param: n is the lower state, principle qm
    @param: l is the orbital angular momentum l = 0,1,...n-1
    @param: j is the total angular momentum j = S+l
    @param: mj: secondary total angular momentum quantum number mj = -j, -j+1...j-1, j
    @param: verbose = false wont print, if true willl print
 
    '''  
 
"""
import sys
sys.path.append('C:/Users/liamw/Desktop/RYDBERG_GIT-1') 
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/src')
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/src/PhD_Rydberg_code_version_1')
sys.path.append('C:/Users/liamw/Desktop/Rydberg_Git/src/PhD_Rydberg_code_version_1/Springer/Rydberg_Code')

import pandas as pd 
import numpy as np
import liams_windows_config
import csv

from arc import Rubidium85





def state_array_dataframe_generator(state_array: list["state_describing_object"]) -> pd.DataFrame:
    """Creates a dataframe containing the information for the relevant atomic transitions for a given system.

    Args:
        state_array (list[state_describing_object]): An object that will contain the specific atomic infomation.

    Returns:
        pd.DataFrame: A dataframe containing the quantum numbers of all of the transitions of the system.
    """
    #from Rydberg_Cs_example_paper_natural_decay_multiple_lines import state_describing_object
    state_df = pd.DataFrame(columns = ['name','n','l','j','mj'])

    for i in range(len(state_array)):
        single_record = pd.DataFrame({"name": [state_array[i].name],
                         "n": [state_array[i].energy_level],
                         "l": [state_array[i].numerical_L],
                         "j": [state_array[i].J],
                         "mj": [state_array[i].mj],
                         })

        state_df = pd.concat([state_df,single_record])

    return(state_df)

def get_state_lifetime(atom: list[Rubidium85], N :int, L: int ,J: float ,T=300) -> float:
    """_summary_

    Args:
        atom (list[Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        N (int): Principal quantum number.
        L (int): The orbital angular momentum quantum number.
        J (float): The total angular momentum.
        T (int, optional): The temperature of the system in K. Defaults to 300.

    Returns:
        float: The state (natural) lifetime of the atomic transtion, in s.
    """

    #atom = atomic species and natrually occuring or isotopic
    #N = Princple QN, L = Orbital Angular Momentum QN,
    #J= Total angular momentum,
    #T = Temperature in Kelvin
    
    lifetime = [atom[i].getStateLifetime(int(N), 
                                int(L), 
                                J, 
                                temperature = T, 
                                includeLevelsUpTo=int(N+10)
                                )for i in range(len(atom))]
    return lifetime

def transition_picker(state_array: list["state_describing_object"], atom: list[Rubidium85], vbose=False) -> pd.DataFrame:
    """ This is used to create a Dataframe of atomic state the user has entered and find the state lifetime assoicated with those states.  

    Args:
        state_array (list[state_describing_object]): An object that will contain the specific atomic infomation.
        atom (list[Rubidium85]): This is a list containing the ARC object of Rb85, but may contain any relevant atomic species. This must be an ARC atom class.
        vbose (bool, optional): A flag set for writing CSVs. Do not change this value. Defaults to False.

    Returns:
        pd.DataFrame: Dataframe containing information on the atomic transitions of interest and calculate the state lifetime up N+10 levels
    """

    config = liams_windows_config.liams_windows_config()  

    state_df = state_array_dataframe_generator(state_array)
    
    name_list = [atom[i].__str__().split(".")[2].split(" ")[0]+(' in s') 
                            for i in range(len(atom))]
    
    state_df1 = pd.DataFrame(columns=['name','n','l','j','mj'])
    
    state_df1=pd.concat([state_df1, state_df])
    
    state_df1[name_list] = state_df.apply(lambda x: get_state_lifetime(atom, x.n, x.l, x.j), axis= 1, result_type = 'expand')


    return state_df1
    
if __name__ =='__main__':
    
    print('this should be inported into main file')
    pass 