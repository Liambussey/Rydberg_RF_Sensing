import qutip as qp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

from scipy.constants import k as k_b
from scipy.constants import c as c_c ## speed of light 
from scipy.constants import h as c_h # planck constant
from scipy.constants import hbar as hbar #reduced planks

#from Rydberg_Cs_example_paper_natural_decay_multiple_lines import laser_parameter


from typing import Union

def basis_generator(number_states: int) -> list[qp.Qobj]:
    """This generates the number of Fock basis state

    Args:
        number_states (int): This is the Interger value for the total number of atomic state in the system.

    Returns:
        list[qp.Qobj]: List of N = number_states Fock states 
    """
    #### basis state for a N system
    basis_list =[]

    for i in range(number_states):
        basis = qp.basis(number_states, i)
        basis_list.append(basis)

    return(basis_list)

def atomic_operators(basis: list[qp.Qobj], number_states: int) -> list[qp.Qobj]:
    """This function uses the basis state and number of state to create the atomic jump operators

    Args:
        basis (list[qp.Qobj]): List of N = number_states Fock states 
        number_states (int): This is the interger value for the total number of atomic state in the system.

    Returns:
        list[qp.Qobj]: list of creation and annihilation operators 
    """
    #### sigma operator or atomic jump operator 
    sigma_operator = []#np.zeros(shape = ((number_state),(number_state)))

    for i in range(number_states):
        for b in range(number_states):
            sigma = basis[i]*basis[b].dag()
            
            sigma_operator.append(sigma)
            
    return(sigma_operator)

def Hamil_builder(number_states: int, Trans_wave: pd.DataFrame, dephasing: list[float]) -> np.array:
    """_summary_

    Args:
        number_states (int): This is the interger value for the total number of atomic state in the system.
        Trans_wave (pd.DataFrame): Dataframe with the expanded information on the excitation pathway gathered using ARC
        dephasing (list[float]): This is a list of floats containing the dephasing scanning range and number of steps

    Returns:
        np.array: This 2D array contains the an nxn Hamiltionian 
    """
    Hamiltonian = np.zeros(shape=((number_states),(number_states)))                     ### NxN matrix of zeros to form the base Hamiltonian 
    

 #### how to build this a function builder 
    for i in range (number_states- 1):
        
        Hamiltonian[i][i+1] = -(1/(4*np.pi))*(Trans_wave.iloc[i]['rabi Freq rads/s'])   ### fills the upper off diagonal element for an NXN matrix 
        Hamiltonian[i+1][i] = -(1/(4*np.pi))*(Trans_wave.iloc[i]['rabi Freq rads/s'])   ### fills the lower off diagonal element for an NXN matrix 

    for i in range (number_states):
        if i == 0: 
            Hamiltonian[i][i] = 0.0
        else:
            Hamiltonian[i][i] = dephasing[i-1]                                            # forms the diagonal element of the system
    
    return(Hamiltonian)

def collapse_op(number_state: float, basis: list[qp.Qobj], decay_wave: pd.DataFrame) -> np.array :
    """Natrual Decay collapse operator 

    Args:
        number_state (float): This is the interger value for the total number of atomic state in the system.
        basis (list[qp.Qobj]): List of N = number_states Fock states 
        decay_wave (pd.DataFrame): This a Dataframe with the expanded information on the decay pathway gathered using ARC.

    Returns:
        np.array: This return an array of collapse operature generaterd using the natural decay taken from ARC
    """
    natural_line =[]
    for i in range(number_state - 1):
        #spontaneous decay in Hz
        natural_line.append((decay_wave.iloc[i]['Transition Rate (Hz) with BBR for Rb']))

    ls = []
    dim = number_state
    n_l = []
    
    for j in range(1, dim):
        op = np.zeros((dim, dim))
        natural_line_value = natural_line[j-1]

        op[j-1,j] = natural_line_value
        ls.append(op)
    
    ls = np.array(ls)

    return(ls)
    

def collapse_op_laser(number_state: float, laser_trans: list["laser_parameter"], basis: list[qp.Qobj], Trans_wave: pd.DataFrame) -> list:
    """Laser linewidth collapse operator, allow for the braodening of the transition based on the laser linewidth, this will set the RF linewidth to zero

    Args:
        number_state (float): This is the interger value for the total number of atomic state in the system.
        laser_trans (list[laser_parameter]): This is a list of field parameter object, all of the key experimental parameters are used here. These will also set the numerical 'experiment' parameters, e.g. which wavelengths are being scanned over.
        basis (list[qp.Qobj]): List of N = number_states Fock states 
        Trans_wave (pd.DataFrame): This a Dataframe with the expanded information on the excitation pathway gathered using ARC.

    Returns:
        list: returns a list of collapse operator for the laser linewidth
    """
 #laser linewidth in Hz
    laser_line =[]
    for i in range(number_state - 1):
        if  Trans_wave.iloc[i]['feild'] == 'R':
            laser_line.append(0)
        else:
            laser_line.append((Trans_wave.iloc[i]['linewidth']))
    ls = []
    dim = number_state
    q_l = []
    op = np.zeros((dim, dim)) 

    for j in range(1, dim):
        op = np.zeros((dim, dim))
        laser_line_value = laser_line[:j]
        laser_line_value.reverse()
        
        for i in range(0, j):
            op[i, j] = (np.sum(laser_line_value[:(j-i)]))
            
        ls.append(op)
    return ls
  
def dephase_dataframe_generator(laser_trans: list["laser_parameter"], Trans_wave: pd.DataFrame)->pd.DataFrame:
    """This function merges a set of list parameter and items store in a dataframe into a new datafram to allow for easier removal of data in other functions

    Args:
        laser_trans (list[laser_parameter]): This is a list of field parameter object, all of the key experimental parameters are used here. These will also set the numerical 'experiment' parameters, e.g. which wavelengths are being scanned over.
        Trans_wave (pd.DataFrame): This a Dataframe with the expanded information on the excitation pathway gathered using ARC.

    Returns:
        pd.DataFrame: DataFrame contain infromation on the resonance, scanning source, dephasing bounds, wavelength, speed, frequency, mass and beam propigation direction
    """
    dephase_df = pd.DataFrame(columns=['onress', 'dephase_lower',
                                       'dephase_upper',"scanning_source",
                                       'wavelength',"mean_speed",'beam_prop'])

    for i in range(len(laser_trans)):
        single_record = pd.DataFrame({"onress": [laser_trans[i].onress],
                         "dephase_lower": [laser_trans[i].dephase_lower],
                         "dephase_upper": [laser_trans[i].dephase_upper],
                         "scanning_source": [laser_trans[i].scanning_source],
                         "wavelength": [Trans_wave.iloc[i]['wavelength in nm for Rb']/(1*10**9)],
                         "mean_speed": [Trans_wave.iloc[i]['mean_speed']],
                         'transf': [Trans_wave.iloc[i]['freq in Hz for Rb']],
                         'mass': [Trans_wave.iloc[i]['mass']],
                         'beam_prop': [laser_trans[i].beam_prop]})
                          
        dephase_df = pd.concat([dephase_df,single_record])             
    
    return dephase_df

def dephase(laser_trans: list["laser_parameter"], number_states: float, dephasing_step: float, Trans_wave: pd.DataFrame) -> list:
    """creates the dephasing list based on the number of steps and boundaries.

    Args:
        laser_trans (list[laser_parameter]): This is a list of field parameter object, all of the key experimental parameters are used here. These will also set the numerical 'experiment' parameters, e.g. which wavelengths are being scanned over.
        number_states (float): This is the interger value for the total number of atomic state in the system.
        dephasing_step (float): This is how may step are used to creat the dephasing scan between the upper and lower scan freqeucy
        Trans_wave (pd.DataFrame): This a Dataframe with the expanded information on the excitation pathway gathered using ARC.

    Returns:
        list: list of dephasing values bound between the lower and upper bounds 
    """
    dephase_df = dephase_dataframe_generator(laser_trans, Trans_wave)
    #print(dephase_df.scanning_source)
    dephasing_range = np.zeros(shape=(dephasing_step, (number_states-1)))
    dephasing = np.zeros(shape=(dephasing_step, (number_states-1)))
    dephasing_range_concat =[]
    
    for i in range(len(dephase_df)):
        (dephasing_step/3)
        upper= dephase_df.iloc[i]['dephase_upper']
        lower = dephase_df.iloc[i]['dephase_lower']
       
        if upper == 0 and lower == 0:
            step = dephasing_step    
            Dephasing_range = np.zeros(shape=(step))
            dephasing_range[:, i] = Dephasing_range
        if dephase_df.iloc[i]['scanning_source'] == 'N' and dephase_df.iloc[i]['onress'] == 'N':
            step = dephasing_step  
            Dephasing_range = [upper]*step
            dephasing_range[:, i] = Dephasing_range    
        
        elif dephase_df.iloc[i]['scanning_source'] == 'Y' and dephase_df.iloc[i]['onress'] == 'N':
            step = int(dephasing_step/3)
            Dephasing_range = np.linspace([lower, -25*10**6,  25*10**6], [-25*10**6, 25*10**6, upper],step, axis = 1)
            Dephasing_range.reshape(dephasing_step)
            dephasing_range_concat = np.ravel(Dephasing_range) 
            dephasing_range[:, i] = dephasing_range_concat


        
        
    #print(dephasing_range)
        
    for j in range(number_states - 1):
        if j == 0:
            #dephasing[:,j] = Dephasing_range
            dephasing[:,j] = dephasing_range[:,j]
        else:
            dephasing[:,j] = (dephasing[:,j-1]+dephasing_range[:,j])#dephasing_range[:,j])

    dephasing = dephasing
    return dephasing

def EIT_builder(number_state: float, Trans_wave: pd.DataFrame, decay_wave: pd.DataFrame, laser_trans: list["laser_parameter"], dephasing_step: float) -> Union[pd.DataFrame, list]:
    """This function calls multiple subfunction required build the desire hamiltonian and collapse operators to solve the steady state solution of the density matrix using Qutip steadystate solve.

    Args:
        number_state (float): This is the interger value for the total number of atomic state in the system.
        Trans_wave (pd.DataFrame): This a Dataframe with the expanded information on the excitation pathway gathered using ARC.
        decay_wave (pd.DataFrame): This a Dataframe with the expanded information on the decay pathway gathered using ARC.
        laser_trans (list[laser_parameter]): This is a list of field parameter object, all of the key experimental parameters are used here. These will also set the numerical 'experiment' parameters, e.g. which wavelengths are being scanned over.

        dephasing_step (float): This is how may step are used to creat the dephasing scan between the upper and lower scan freqeucy.

    Returns:
        Union[pd.DataFrame, list]: This return a pd.Dataframe containing the steady solution to the density matrix rho_12, rho_13.. etc. and a list containing the dephasing scan to be plotted over. 
    """

    ### bais is the ket vector for the number of states in the atom 
    basis_list = basis_generator(number_state)

    #### atomic operator is the jump operator ie basis[0]*basis[1].dag() this will from an nxn matrix 
    atomic_operator = atomic_operators(basis_list, number_state)
    
    ## this calls the collapse function to form the collapse operator for the system
    collapse_2 = collapse_op(number_state, basis_list, decay_wave) # Natural decay 

    collapse_l = collapse_op_laser(number_state, laser_trans, basis_list, Trans_wave) # laser linewidth
    
    collapse_n_l_l =[]
    collapse_n_list =[]
    
    for i in range(number_state-1):
        
        collapse_n = collapse_2[i]
        
        collapse_n_list.append(-np.sqrt(collapse_n))#
        
        collapse_n_l =  -((np.sqrt(collapse_l[i]+collapse_2[i]))) #
        collapse_n_l_l.append((collapse_n_l))
        
    collapse_n_l_q = []
    collapse_n_q = []
    
    
    for i in range(number_state-1):
        collapse_n_q.append(qp.Qobj(collapse_n_list[i]))
        collapse_n_l_q.append(qp.Qobj(collapse_n_l_l[i]))
 

    
    chi_1 = pd.DataFrame()### first column electric susceptability for the |0> ---> |1>
                        ### second column is chi for the |1> ---> |2>
                        ### third column is chi for the |2>--> |3> and so on
    
    chi_2 = pd.DataFrame()
    
    
    chi_value = np.zeros(shape=(dephasing_step,(number_state-1)), dtype=complex)
    chi_value_2 = np.zeros(shape=(dephasing_step,(number_state-1)), dtype=complex)
    
    chi_test_1 = np.zeros(shape=(dephasing_step,(number_state-1)), dtype=complex)
    chi_test_2 = np.zeros(shape=(dephasing_step,(number_state-1)), dtype=complex)

    dephasing_2 = dephase(laser_trans,number_state, dephasing_step, Trans_wave)
     
    for i in range(dephasing_step):
        
       
        ##non doppler Hamiltonian 
        Hamiltonian = (Hamil_builder(number_state, Trans_wave, dephasing_2[i, :])) # this calls the hamiltonian builder function

        quantum_obj_2 = qp.Qobj(Hamiltonian)

        ##### steady state solution to non doppler hamiltonian
        rho_ss_1 = qp.steadystate(quantum_obj_2, collapse_n_q) #natural
        #rho_ss_2 = qp.steadystate(quantum_obj_2, collapse_n_l_q) # natrual +laser
        
        ## expectation value of systemm
        for b in range(1, number_state):

            a = b + (number_state*(b-1))
           ## non doppler expectation value 

            chi_test_1[i, b-1] = qp.expect(atomic_operator[a], rho_ss_1) #natural
            #chi_test_2[i, b-1] = qp.expect(atomic_operator[a], rho_ss_2) # natrual +laser
            
            
    # electric susceptability for the off diagonal elements 
    
    chi_1 = chi_test_1 # natrual 
    #chi_2 =chi_test_2 # natrual +laser
    
    dephasing_graph = [] # transferes the probe dephasing
    
    laser_state_df = dephase_dataframe_generator(laser_trans, Trans_wave)

    # dephasing list not 100% needed
    for i in range(number_state-1):
        if ((laser_state_df.iloc[i]['onress'] == 'N') & (laser_state_df.iloc[i]['scanning_source'] == 'Y')):
            dephasing_graph = dephasing_2[:,i]#*2*np.pi  # transferes the probe dephasing (need to work out how to transfere the relevent dephasing)

    return(chi_1, dephasing_graph)
   
if __name__ == '__main__': 
    print('this should be inported into main file')