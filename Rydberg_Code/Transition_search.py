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

import pandas as pd 
import numpy as np
import liams_windows_config
import csv

def state_array_dataframe_generator(state_array):

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

def get_state_lifetime(atom,N,L,J,T=300):
    #atom = atomic species and natrually occuring or isotopic
    #N = Princple QN, L = Orbital Angular Momentum QN,
    #J= Total angular momentum,
    #T = Temperature in Kelvin
    
    val = [atom[i].getStateLifetime(int(N), 
                                int(L), 
                                J, 
                                temperature = T, 
                                includeLevelsUpTo=int(N+10)
                                )for i in range(len(atom))]
    return(val)

def transition_picker(state_array, atom, vbose = False):

    config = liams_windows_config.liams_windows_config()  

    state_df = state_array_dataframe_generator(state_array)
    
    name_list = [atom[i].__str__().split(".")[2].split(" ")[0]+(' in s') 
                            for i in range(len(atom))]
    
    state_df1 = pd.DataFrame(columns=['name','n','l','j','mj'])
    
    state_df1=pd.concat([state_df1, state_df])
    
    state_df1[name_list] = state_df.apply(lambda x: get_state_lifetime(atom, x.n, x.l, x.j), axis= 1, result_type = 'expand')


    return(state_df1)  
    
if __name__ =='__main__':
    
    print('this should be inported into main file')
    pass 