import numpy as np
import pandas as pd

############################
##### Helper Functions #####
############################

def wenn(condition, option1, option2):
    """
    This function implements the ternary operator 'wenn' in Python. It takes three arguments: a condition, an option1 
    if the condition is true, and an option2 if the condition is false or not satisfied. It is used as helper function 
    simplyfing the implementation from excel.
    
    Parameters:
        condition (bool): The condition to check.
        option1 : The return value when the condition is True.
        option2 : The return value when the condition is False.
        
    Returns:
        Evaluated based on condition, returns either 'option1' or 'option2'. 
    """
    if condition:
        return option1
    else:
        return option2

def compute_QKL_steel_I(cls_type,N_Ed, M_yEd, f_y, M_zEd=0):
    """
    Computes the quality factor of a steel I.
    
    Parameters:
        cls_type (str): Type of I element
        N_Ed (float): Load
        M_yEd  (float): Load
        f_y (float): 
        M_zEd (float, optional): Load (default is 0)
        
    Returns:
        int : Quality factor.
    """
    
    #cls_type: class of the steel (e.g. "HEA 160")
    try:
        qs_info = get_steel_I_QS_info(cls_type)
        r = float(qs_info['r'])
        b = float(qs_info['b'])
        d = float(qs_info['d'])
        t_w = float(qs_info['tw'])
        t_f = float(qs_info['tf'])
        W_ely = float(qs_info['Wel,y'])
        h = float(qs_info['h'])
        A = float(qs_info['A'])
        c = (b - 2 * r - t_w) / 2

        # Case 1: H22 Querschnitt lediglich auf Zug beansprucht?
        if N_Ed < 0 and M_yEd == 0 and M_zEd == 0:
            return 1

        # Case 2: H27 Querschnitt lediglich auf Druck beansprucht?
        if N_Ed > 0 and M_yEd == 0 and M_zEd == 0:
            if d/t_w <= 33 and c/t_f <= 9:
                return 1
            elif d/t_w <= 38 and c/t_f <= 10:
                return 2
            elif d/t_w <= 42 and c/t_f <= 14:
                return 3
            else:
                return 4

        # Case 3: Q27 Querschnitt lediglich auf Biegung um die y-Achse beansprucht?
        if N_Ed == 0 and M_yEd > 0 and M_zEd == 0:
            if d/t_w <= 72 and c/t_f <= 9:
                return 1
            elif d/t_w <= 83 and c/t_f <= 10:
                return 2
            elif d/t_w <= 124 and c/t_f <= 14:
                return 3
            else:
                return 4

        # Case 4: Y27 Querschnitt lediglich auf Biegung um die z-Achse beansprucht?
        if N_Ed == 0 and M_yEd == 0 and M_zEd > 0:
            pass # not implemented as M_zEd = 0

        # Case 5: Q55 Querschnitt auf Biegung um die y-Achse und N beansprucht?
        if N_Ed != 0 and M_yEd > 0 and M_zEd == 0:
            L62 = np.abs(N_Ed * 1000/f_y)
            if N_Ed > 0:
                L70 = (d/2 + (L62/t_w)/2)/d
            else:
                L70 = (d/2 - (L62/t_w)/2)/d

            # Checking QKL 1
            if (np.abs(N_Ed*1000/f_y)/t_w) < d: # P65
                if L70 > 0.5:
                    Q76 = d/t_w <= 396*1/(13*L70-1)
                else:
                    Q76 = d/t_w <= 36*1/L70
                P79 = c/t_f <= 9

                if Q76 == True and P79 == True: # K80 = "Einordnung in QKL1 möglich"
                    return 1
                else:
                    pass
            else:
                if (d/t_w <= 33 and c/t_f <= 9) or N_Ed <= 0: # K99 = "Einordnung in QKL1 möglich"
                    return 1
                else:
                    pass
            
            # Checking QKL 2
            if (np.abs(N_Ed*1000/f_y)/t_w) < d: # P65
                if L70 > 0.5:
                    Q86 = d/t_w <= 456*1/(13*L70-1)
                else:
                    Q86 = d/t_w <= 41.5*1/L70
                P89 = c/t_f <= 10

                if Q86 and P89: # 90 = "Einordnung in QKL2 möglich"
                    return 2
                else:
                    pass
            else:
                if (d/t_w <= 38 and c/t_f <= 10) or N_Ed <= 0: # 106 = "Einordnung in QKL2 möglich"
                    return 2
                else:
                    pass
            
            # Checking QKL 3 
            M118 = M_yEd * 1000000 / W_ely * 1000 / h * d + N_Ed * 1000 / A * 100
            M121 = - M_yEd * 1000000 / W_ely * 1000 / h * d + N_Ed * 1000 / A * 100
            if np.max(M118, M121) > 0:
                M125, M126 = np.max(M118, M121), np.min([M118, M121])
                K128 = M126/M125
                N123 = M118 <= 0 and M121 <= 0

                if N123:
                    Q142 = True
                else:
                    if K128 > 1:
                        Q142 = (not N123) and ((d/t_w) <= (42/(0.62+0.33*K128)))
                    else:
                        Q142 = (not N123) and (K128 <= 0) and (c/t_f <= 14)

                if N123:
                    P145 = True
                else:
                    P145 = c/t_f <= 14

                if Q142 and P145:
                    return 3
                else:
                    pass
            else:
                return 3
            
            return 4

        # Case 6: Y55 Querschnitt auf Biegung um die z-Achse und N beansprucht?
        if N_Ed != 0 and M_yEd == 0 and M_zEd != 0:
            pass # not implemented as M_zEd = 0

        # Case 7: AG27 Querschnitt auf Biegung um die y- sowie die z-Achse beansprucht (ohne N)?
        if N_Ed == 0 and M_yEd > 0 and M_zEd > 0:
            pass # not implemented as M_zEd = 0

        # Case 8: AG55 Querschnitt auf Biegung um die y- sowie die z-Achse beansprucht (mit N)?
        if N_Ed != 0 and M_yEd > 0 and M_zEd > 0:
            pass # not implemented as M_zEd = 0
        
        raise NotImplementedError("Error occured while computing QKL. Load characteristics not known!")
    except Exception as err:
        raise ValueError(f"Error occured while computing QKL: {err}")
    
def compute_QKL_steel_O(cls_type,f_y):
    """
    Computes the quality factor of a steel O.
    
    Parameters:
        cls_type (str): Type of O element
        f_y (float): 
        
    Returns:
        int : Quality factor.
    """
    try:
        qs_info = get_steel_O_QS_info(cls_type)
        D = float(qs_info['D'])
        t = float(qs_info['T'])

        C7 = f_y
        D27 = np.sqrt(235/C7)
        C11 = t
        C10 = D
        C17 = C10/C11
        qkl = wenn(C17<=50*D27**2,1,wenn(C17<=70*D27**2,2,wenn(C17<=90*D27**2,3,4)))
    except Exception as err:
        raise ValueError(f"Error occured while computing QKL: {err}")
    return qkl

def compute_QKL_steel_RH(f_y, r_i, t, h, b):
    """
    Computes the quality factor of a steel RH.
    
    Parameters:
        f_y (float): 
        ri_i (float): 
        t (float): 
        h  (float): 
        b (float): 
        
    Returns:
        int : Quality factor.
    """
    epsilon = np.sqrt(235 / f_y)
    c_druck = h - 2 * t - 2 * r_i 
    t_druck = t
    c_biege = b - 2 * t - 2 * r_i
    t_biege = t

    if c_druck/t_druck <= 33 * epsilon and c_biege/t_biege <= 72 * epsilon:
        return 1
    elif c_druck/t_druck <= 38 * epsilon and c_biege/t_biege <= 83 * epsilon:
        return 2
    elif c_druck/t_druck <= 42 * epsilon and c_biege/t_biege <= 124 * epsilon:
        return 3
    else:
        return 4
    
def get_steel_I_QS_info(qs):
    """
    Retrieves information related to  a steel I element based on its type qs.
    
    Parameters:
        qs (str): Type
        
    Returns:
        qs_info (dict): A dictionary containing various properties of the steel element 
    """
    try:
        df = pd.read_csv("data/static_computation/steel/QS_Stahl_Profile_I.csv")
        qs_info = df[df['id'] == qs]
        if qs_info.empty:
            raise NotImplementedError(f"No QS information availabel for type {qs}!")
        else:
            return qs_info.iloc[0]
    except Exception as err:
        print(f"Error while retrieving QS data: {err}")
        raise   

def get_steel_O_QS_info(qs):
    """
    Retrieves information related to  a steel O element based on its type qs.
    
    Parameters:
        qs (str): Type
        
    Returns:
        qs_info (dict): A dictionary containing various properties of the steel element 
    """
    try:
        df = pd.read_csv("data/static_computation/steel/QS_Stahl_Profile_O.csv") 
        qs_info = df[df['section'] == qs]
        if qs_info.empty:
            raise NotImplementedError(f"No QS information availabel for type {qs}!")
        else:
            return qs_info.iloc[0]
    except Exception as err:
        print(f"Error while retrieving QS data: {err}")
        raise   

def get_steel_RH_QS_info(qs):
    """
    Retrieves information related to  a steel RH element based on its type qs.
    
    Parameters:
        qs (str): Type
        
    Returns:
        qs_info (dict): A dictionary containing various properties of the steel element 
    """
    try:
        df = pd.read_csv("data/static_computation/steel/QS_Stahl_Profile_Hohl.csv")
        qs_info = df[df['section'] == qs]
        if qs_info.empty:
            raise NotImplementedError(f"No QS information availabel for type {qs}!")
        else:
            return qs_info.iloc[0]
    except Exception as err:
        print(f"Error while retrieving QS data: {err}")
        raise   

########################
##### Steel Column #####
########################

def only_static_steel_I_column(cls_type, f_y, y_M1, L, E, v, L_T): 
    """
    This function computes the static loads for a steel column of type I.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        L (float): 
        E (float): 
        L_T (float): 
        v  (float):
    
    Returns:
        N_bRdy (float): Maximal load
        N_bRdz (float): Maximal load
        N_bRdT (float): Maximal load
    """

    try:
        qs_info = get_steel_I_QS_info(cls_type)
       
        t_w = float(qs_info['tw'])
        t_f = float(qs_info['tf'])
        h = float(qs_info['h'])
        b = float(qs_info['b'])
        r = float(qs_info['r'])
        c = (b - 2 * r - t_w) / 2

        z_p = -h/2
        # N_b,Rd,y
        J36 = h/b > 1.2 and t_f <= 40
        J38 = h/b <= 1.2 and (t_f <= 100 and t_f > 40)
        J40 = h/b <= 1.2 and t_f <= 100
        J42 = h/b <= 1.2 and t_f > 100
        H33 = f_y <= 420
        I33 = f_y >= 460 and f_y <= 700

        K36 = wenn(H33 and J36, "a", wenn(I33 and J36, "a_0", False))
        K38 = wenn(H33 and J38, "b", wenn(I33 and J38, "a", False))
        K40 = wenn(H33 and J40, "b", wenn(I33 and J40, "a", False))
        K42 = wenn(H33 and J42, "d", wenn(I33 and J42, "c", False))

        D46 = wenn(K36, K36, wenn(K38, K38, wenn(K40, K40, wenn(K42, K42, False))))

        F46 = wenn(D46 == "a_0", 0.13, wenn(D46 == "a", 0.21, wenn(D46 == "b", 0.34, wenn(D46 == "c", 0.49, 0.76))))

        F15 = ((1/1)*(c*t_f**3/3+ c* t_f*( h-t_f)**2+(b-2* c)*t_f**3/48+t_f/4*( h- t_f/2)**2*( b-2* c))+(2/2)*(( h-2* t_f)**3*t_w/12+( r**4/12- r**4/2*(np.pi/8-8/(9*np.pi))+r**2*( h/2- t_f- r/2)**2- r**2*(1/4*np.pi)*(h/2- t_f- r+4/3* r/np.pi)**2)*4)+(3/3)*(( b-2*c)* t_f**3/48+( b-2* c)* t_f/4*(h-3* t_f/2)**2))/(10**4)
        C15 = ((1/1)*(t_f*c*4+(b-2*c)*t_f)+(2/2)*((h-2*t_f)*t_w+r**2*(4-np.pi))+(3/3)*((b-2*c)*t_f))/10**2
  
        E100 = (E * F15 * 10000 * np.pi**2 / (L_T * 1000)**2) / 1000

        E99 = np.sqrt(((C15 * 100 * f_y)/1000)/E100)
        E98 = 0.5 * (1 + F46 * (E99 - 0.2) + E99**2)
        E97 = np.min([1 / (E98 + np.sqrt(E98**2-E99**2)), 1])

        N_bRdy = (E97 * C15 * 100 * f_y/y_M1)/1000

        # N_b,Rd,z
        K37 = wenn(H33 and J36, "b", wenn(I33 and J36, "a_0", False))
        K39 = wenn(H33 and J38, "c", wenn(I33 and J38, "a", False))
        K41 = wenn(H33 and J40, "c", wenn(I33 and J40, "a", False))
        K43 = wenn(H33 and J42, "d", wenn(I33 and J42, "c", False))

        D47 = wenn(K37, K37, wenn(K39, K39, wenn(K41, K41, wenn(K43, K43, False))))

        F47 = wenn(D47 == "a0", 0.13, wenn(D47 == "a", 0.21, wenn(D47 == "b", 0.34, wenn(D47 == "c", 0.49, 0.76))))

        I15 = ((1/1)*(t_f*c**3/3+(b-c)**2*t_f*c+1/12*(b-2*c)**3*t_f)+(2/2)*((h-2*t_f)*t_w**3/12+(r**4/12-r**4/2*(np.pi/8-8/(9*np.pi))+r**2*(1/2*t_w+r/2)**2-r**2*(1/4*np.pi)*(1/2*t_w+r-4/3*r/np.pi)**2)*4)+(3/3)*((b-2*c)**3*t_f/12))/10**4
   
        E108 = (E * I15 * 10**4 * np.pi**2 / (L * 1000)**2) / 1000

        E107 = np.sqrt(((C15 * 100 * f_y) / 1000) / E108)
     
        E106 = 0.5 * (1 + F47 * (E107 - 0.2) + E107**2)
 
        E105 = np.min([1 / (E106 + np.sqrt(E106**2 - E107**2)), 1])
   
        N_bRdz = (E105 * C15 * 100 * f_y/y_M1)/1000

        # N_b,Rd,T
        Iyy = ((1/1)*(c*t_f**3/3+c*t_f*(h-t_f)**2+(b-2*c)*t_f**3/48+t_f/4*(h-t_f/2)**2*(b-2*c))+(2/2)*((h-2*t_f)**3*t_w/12+(r**4/12-r**4/2*(np.pi/8-8/(9*np.pi))+r**2*(h/2-t_f-r/2)**2-r**2*(1/4*np.pi)*(h/2-t_f-r+4/3*r/np.pi)**2)*4)+(3/3)*((b-2*c)*t_f**3/48+(b-2*c)*t_f/4*(h-3*t_f/2)**2))/(10**4)
        Izz = ((1/1)*(t_f*c**3/3+(b-c)**2*t_f*c+1/12*(b-2*c)**3*t_f)+(2/2)*((h-2*t_f)*t_w**3/12+(r**4/12-r**4/2*(np.pi/8-8/(9*np.pi))+r**2*(1/2*t_w+r/2)**2-r**2*(1/4*np.pi)*(1/2*t_w+r-4/3*r/np.pi)**2)*4)+(3/3)*((b-2*c)**3*t_f/12))/10**4
        C27 = ((2/3)*(b-0.63*t_f)*t_f**3+(1/3)*(h-2*t_f)*t_w**3+2*(t_w/t_f)*(0.145+0.1*r/t_f)*(((r+t_w/2)**2+(r+t_f)**2-r**2)/(2*r+t_f))**4)/(10**4)
        C28 = ((t_f*b**3)/24*(h-t_f)**2)/(10**6)
        F16 = np.sqrt(Iyy/C15)
        I16 = np.sqrt(Izz/C15)
        E126 = E /(2*(1+v))
        E125 = F16**2 + I16**2
        E124 = (1/(E125 * 10**2)*(E126 * C27 * 10**4 + np.pi**2 * E * C28 * 10**6/(L_T * 1000)**2))/1000
        E132 = E125*10**2/(2*((F16*10)**2+(I16*10)**2))*(E100+E124-np.sqrt((E100+E124)**2-4*E100*E124*((F16*10)**2+(I16*10)**2)/(E125*10**2)))    
        E123 = np.min([E124,E132])
        E122 = np.sqrt(((C15 * 10**2 * f_y)/1000)/E123)
        E121 = F47
        E120 = 0.5*(1 + E121 * (E122 - 0.2) + E122**2)
        E119 = np.min([1 / (E120 + np.sqrt(E120**2-E122**2)), 1])
        N_bRdT = (E119 * C15 * 10**2 * f_y/y_M1)/1000

    except: 
        print("Error occured while computing static steel column!")
        raise
   
    return N_bRdy, N_bRdz, N_bRdT

def compute_static_steel_I_column(cls_type, f_y, y_M1, L, E, v, L_T, # Information about the element
                                  N_Ed): # Information about loads
    """
    Wrapper function to compute the static loads for a steel column of type I.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L (float): 
        L_T (float): 
        v (float): 
        N_Ed (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """

    N_bRdy, N_bRdz, N_bRdT = only_static_steel_I_column(cls_type,f_y, y_M1, L/1000, E, v, L_T/1000) # /1000 as input is in mm
    # Nachweise
    if N_Ed/N_bRdy > 1:
        return False
    if N_Ed/N_bRdz > 1:
        return False
    if N_Ed/N_bRdT > 1:
        return False
    return True

def only_static_steel_RH_column(cls_type, f_y, y_M1, E, L_Cry):
    """
    This function computes the static loads for a steel column of type RH.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L_Cry (float): 
    
    Returns:
        N_bRdy (float): Maximal load
        N_bRdz (float): Maximal load
    """

    try:
        qs_info = get_steel_RH_QS_info(cls_type)
        H = float(qs_info['H'])
        B = float(qs_info['B'])
        T = float(qs_info['T'])
        r_o = float(qs_info['ro'])
        r_i = float(qs_info['ri'])
        temp = str(qs_info['temp'])
        A = (2 * T * (B + H - 2 * T) - (4 - np.pi) * (r_o**2 - r_i**2)) / 100

        C31 = y_M1
        C30 = f_y
        C13 = A
        C34 = E
        C32 = L_Cry
        C7 = H
        C8 = B
        C9 = T
        C10 = r_o
        C11 = r_i
        D26 = wenn(temp == "w", wenn(C30 < 460, "a", wenn(C30 >= 460 and C30 <= 700,"a0", False)), "c")
        D27 = D26
        if not D26: raise ValueError(f"Maximal possible value for f_y is 700. However, the input was {C30}!")

        P17 = (1/3-np.pi/16-1/(3*(12-3*np.pi)))*C10**4
        P13 = (1-np.pi/4)*C10**2
        P15 = C7/2-(10-3*np.pi)/(12-3*np.pi)*C10
        P18 = (1/3-np.pi/16-1/(3*(12-3*np.pi)))*C11**4
        P14 = (1-np.pi/4)*C11**2
        P16 = (C7-2*C9)/2-(10-3*np.pi)/(12-3*np.pi)*C11
        F13 = 1/(10**4)*((C8*C7**3)/12-((C8-2*C9)*(C7-2*C9)**3)/12-4*(P17+P13*P15**2)+4*(P18+P14*P16**2))
        F26 = wenn(D26=="a0",0.13,wenn(D26=="a",0.21,wenn(D26=="b",0.34,wenn(D26=="c",0.49,0.76))))
        E60 = (C34*F13*10**4*np.pi**2/(C32*1000)**2)/1000
        E59 = np.sqrt(((C13*10**2*C30)/1000)/E60)
        E58 = 0.5*(1+F26*(E59-0.2)+E59**2)
        E57 = np.min([1/(E58+np.sqrt(E58**2-E59**2)),1])
        N_bRdy = (E57*C13*10**2*C30/C31)/1000

        S16 = (C8-2*C9)/2-(10-3*np.pi)/(12-3*np.pi)*C11
        S15 = C8/2-(10-3*np.pi)/(12-3*np.pi)*C10
        I13 = 1/(10**4)*((C7*C8**3)/12-((C7-2*C9)*(C8-2*C9)**3)/12-4*(P17+P13*S15**2)+4*(P18+P14*S16**2))
        F27 = wenn(D27=="a0",0.13,wenn(D27=="a",0.21,wenn(D27=="b",0.34,wenn(D27=="c",0.49,0.76))))
        E68 = (C34*I13*10**4*np.pi**2/(C32*1000)**2)/1000
        E67 = np.sqrt(((C13*10**2*C30)/1000)/E68)
        E66 = 0.5*(1+F27*(E67-0.2)+E67**2)
        E65 = np.min([1/(E66+np.sqrt(E66**2-E67**2)),1])
        N_bRdz = (E65*C13*10**2*C30/C31)/1000
    except: 
        print("Error occured while computing static steel Profile Hohl column!")
        raise
    return N_bRdy,N_bRdz
   
def compute_static_steel_RH_column(cls_type, f_y, y_M1, E, L_Cry, # Information about the element
                                   N_Ed): # Information about loads
    """
    Wrapper function to compute the static loads for a steel column of type RH.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L (float): 
        L_Cry (float): 
        N_Ed (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """

    N_bRdy, N_bRdz = only_static_steel_RH_column(cls_type, f_y, y_M1, E, L_Cry/1000) # /1000 because input is in mm
    # Nachweise
    if N_Ed/N_bRdy > 1:
        return False
    if N_Ed/N_bRdz > 1:
        return False
    return True

def only_static_steel_O_column(cls_type, f_y, y_M1, L_Cr, E, temp_w):
    """
    This function computes the static loads for a steel column of type O.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L_Cr (float): 
        temp_w (bool): Wheter warm or cold forged. True equals warm.
    
    Returns:
        N_bRd (float): Maximal load
    """

    try:
        ksl = wenn(temp_w, wenn(f_y < 460, "a", wenn(f_y >= 460 and f_y <= 700, "a0", "False")), "c")
        alpha = wenn(ksl == "a0", 0.13, wenn(ksl == "a", 0.21, wenn(ksl == "b", 0.34, wenn(ksl == "c", 0.49, 0.76))))

        qs_info = get_steel_O_QS_info(cls_type)
        D = float(qs_info['D'])
        t = float(qs_info['T'])

        C7 = D
        C8 = t
        C29 = E
        C28 = L_Cr
        C27 = y_M1
        C26 = f_y
        F10 = np.pi/64*((C7/10)**4-((C7-2*C8)/10)**4)
        F22 = alpha
        C10 = (np.pi*(C7/2)**2-np.pi*(C7/2-C8)**2)/(10**2)
        E48 = (C29*F10*10**4*np.pi**2/(C28*1000)**2)/1000
        E47 = np.sqrt(((C10*10**2*C26)/1000)/E48)
        E46 = 0.5*(1+F22*(E47-0.2)+E47**2)
        E45 = np.min([1/(E46+np.sqrt(E46**2-E47**2)),1])
        N_bRd = (E45*C10*10**2*C26/C27)/1000
    except: 
        print("Error occured while computing static steel Profile 0 column!")
        raise
    return N_bRd

def compute_static_steel_O_column(cls_type, f_y, y_M1, L_Cr, E, # Information about the element
                                  N_Ed): # Information about loads
    """
    Wrapper function to compute the static loads for a steel column of type O.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L_Cr (float): 
        N_Ed (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """

    split_string = cls_type.split()
    cls_type = split_string[0]
    temp_w = True if len(split_string) == 1 else split_string[1] == "(w)"

    N_bRd = only_static_steel_O_column(cls_type, f_y, y_M1, L_Cr/1000, E, temp_w) # /1000 because input is in mm
    # Nachweise
    if N_Ed/N_bRd > 1:
        return False
    return True

########################
#####  Steel Beam  #####
########################

def only_static_steel_I_beam(cls_type, f_y, y_M1, L, E, z_p, qkl):
    """
    This function computes the static loads for a steel beam of type I.
    
    Parameters:
        cls_type (str): Type of steel beam
        qkl (int): 
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L (float): 
        z_p (float): 
    
    Returns:
        M_bRd (float): Maximal load
        M_yRk (float): Maximal load
    """
    try:
        qs_info = get_steel_I_QS_info(cls_type)        

        t_w = float(qs_info['tw'])
        t_f = float(qs_info['tf'])
        h = float(qs_info['h'])
        b = float(qs_info['b'])
        r = float(qs_info['r'])
        c = (b - 2 * r - t_w) / 2
        k_C = 1
        beta = 0.75
        l_quer_LT_0 = 0.4

        Izz = ((1/1)*(t_f*c**3/3+(b-c)**2*t_f*c+1/12*(b-2*c)**3*t_f)+(2/2)*((h-2*t_f)*t_w**3/12+(r**4/12-r**4/2*(np.pi/8-8/(9*np.pi))+r**2*(1/2*t_w+r/2)**2-r**2*(1/4*np.pi)*(1/2*t_w+r-4/3*r/np.pi)**2)*4)+(3/3)*((b-2*c)**3*t_f/12))/10**4
        E160 = Izz
        I_T = ((2/3)*(b-0.63*t_f)*t_f**3+(1/3)*(h-2*t_f)*t_w**3+2*(t_w/t_f)*(0.145+0.1*r/t_f)*(((r+t_w/2)**2+(r+t_f)**2-r**2)/(2*r+t_f))**4)/(10**4)
        E159 = I_T
        I_W = ((t_f*b**3)/24*(h-t_f)**2)/(10**6)
        E158 = I_W
        E156 = z_p
        E161 = L
        E157 = (E158*10**6+0.039*(E161*1000)**2*E159*10**4)/(E160*10**4)
        I15 = Izz
        C62 = E
        E162 = (C62*I15*10**4*np.pi**2/(E161*1000)**2)/1000
        M_cr = (1.12 * E162 * 1000 * (np.sqrt(E157 + 0.25 * E156**2) + 0.5 * E156)) / (1000*1000) 


        E151 = k_C
        
        #E148 = WENN('Stabilität I kurz'!Q3=WAHR;'Stabilität I kurz'!C23;WENN('Stabilität I kurz'!Q4=WAHR;'Stabilität I kurz'!F24;'Stabilität I kurz'!I24))
        F51 = wenn(h/b<=2, True, False)
        D55 = wenn(F51,"b","c")
        F55 = wenn(D55=="a",0.21,wenn(D55=="b",0.34,wenn(D55=="c",0.49,0.76)))
        Iyy = ((1/1)*(c*t_f**3/3+c*t_f*(h-t_f)**2+(b-2*c)*t_f**3/48+t_f/4*(h-t_f/2)**2*(b-2*c))+(2/2)*((h-2*t_f)**3*t_w/12+(r**4/12-r**4/2*(np.pi/8-8/(9*np.pi))+r**2*(h/2-t_f-r/2)**2-r**2*(1/4*np.pi)*(h/2-t_f-r+4/3*r/np.pi)**2)*4)+(3/3)*((b-2*c)*t_f**3/48+(b-2*c)*t_f/4*(h-3*t_f/2)**2))/(10**4)
        E140 = F22 = ((1/1)*(t_f*c*2*(h-t_f)+(t_w+2*r)*t_f*(h/2-t_f/4))+(2/2)*((t_w*(h-2*t_f)**2)/4+(4-np.pi)/2*r**2*(h-2*t_f)+(3*np.pi-10)/3*r**3)+(3/3)*((t_w+2*r)*t_f*(h/2-3/4*t_f)))/(10**3)
        E141 = F17 = 2*Iyy/h*10
        E142 = wenn(qkl<=2,E140,E141)
        E147 = np.sqrt(E142*10**3*f_y/(M_cr*1000**2))
        E144 = 0.5*(1+F55*(E147-l_quer_LT_0)+beta*E147**2)

        E150 = np.min([1,1-0.5*(1-E151)*(1-2*(E147-0.8)**2)])
        E143 = np.min([1/(E144+np.sqrt(E144**2-beta*E147**2)),1,1/E147**2])
        
        E149 = np.min([1,E143/E150])
        M_bRd = (E149*E142*10**3*f_y/y_M1)/(1000**2)

        E180 = F22 
        G180 = F17 

        M_yRk = (wenn(qkl<3,E180,G180)*1000*f_y)/(1000**2)

    except:
        raise ValueError("Error occured while static beam!")
    return M_bRd, M_yRk

def compute_static_steel_I_beam(cls_type, f_y, y_M1, E, L, z_p,# Information about the element
                                N_Ed, M_Ed, M_yEd, M_zEd=0): # Information about loads
    """
    Wrapper function to compute the static loads for a steel beam of type I.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L (float): 
        z_p (float): 
        N_Ed (float): Load
        M_Ed (float): Load
        M_yEd (float): Load
        M_zEd (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """
    qkl = compute_QKL_steel_I(cls_type, N_Ed, M_yEd, f_y, M_zEd)
    M_bRd, M_yRk = only_static_steel_I_beam(cls_type, f_y, y_M1, L/1000, E, z_p, qkl)
    # Nachweise
    if M_Ed/M_bRd > 1:
        return False
    if M_Ed/M_yRk > 1: 
        return False
    return True

def only_static_steel_RH_beam(cls_type, qkl, y_M1, f_y, E, L_Cry, N_Ed, M_zEd, M_yEd):
    """
    This function computes the static loads for a steel beam of type RH.
    
    Parameters:
        cls_type (str): Type of steel beam
        qkl (int): 
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L_Cry (float): 
        N_Ed (float): Load
        M_zEd (float): Load
        M_yEd (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """

    try:
        qs_info = get_steel_RH_QS_info(cls_type)

        H = float(qs_info['H'])
        B = float(qs_info['B'])
        T = float(qs_info['T'])
        r_o = float(qs_info['ro'])
        r_i = float(qs_info['ri'])
        temp = qs_info['temp']
        C_mz = 0.9
        C_my = 0.9
        A = (2 * T * (B + H - 2 * T) - (4 - np.pi) * (r_o**2 - r_i**2)) / 100

        C35 = qkl
        C31 = y_M1
        G31 = N_Ed
        C30 = f_y        
        C34 = E
        C32 = L_Cry
        C7 = H
        C8 = B
        C9 = T
        C10 = r_o
        C11 = r_i
        C38 = C_my
        H112 = C_mz
        G33 = M_zEd
        G32 = M_yEd
        D26 = wenn(temp == "w", wenn(C30 < 460, "a", wenn(C30 >= 460 and C30 <= 700,"a0", False)), "c")
        if not D26: raise ValueError(f"Maximal possible value for f_y is 700. However, the input was {C30}!")
        D106 = 1

        C13 = (2*C9*(C8+C7-2*C9)-(4-np.pi)*(C10**2-C11**2))/10**2
        P17 = (1/3-np.pi/16-1/(3*(12-3*np.pi)))*C10**4
        P13 = (1-np.pi/4)*C10**2
        P15 = C7/2-(10-3*np.pi)/(12-3*np.pi)*C10
        P18 = (1/3-np.pi/16-1/(3*(12-3*np.pi)))*C11**4
        P14 = (1-np.pi/4)*C11**2
        P16 = (C7-2*C9)/2-(10-3*np.pi)/(12-3*np.pi)*C11
        F13 = 1/(10**4)*((C8*C7**3)/12-((C8-2*C9)*(C7-2*C9)**3)/12-4*(P17+P13*P15**2)+4*(P18+P14*P16**2))
        F26 = wenn(D26=="a0",0.13,wenn(D26=="a",0.21,wenn(D26=="b",0.34,wenn(D26=="c",0.49,0.76))))
        E60 = (C34*F13*10**4*np.pi**2/(C32*1000)**2)/1000
        E59 = np.sqrt(((C13*10**2*C30)/1000)/E60)
        E58 = 0.5*(1+F26*(E59-0.2)+E59**2)
        E57 = np.min([1/(E58+np.sqrt(E58**2-E59**2)),1])

        S16 = (C8-2*C9)/2-(10-3*np.pi)/(12-3*np.pi)*C11
        S15 = C8/2-(10-3*np.pi)/(12-3*np.pi)*C10
        I13 = 1/(10**4)*((C7*C8**3)/12-((C7-2*C9)*(C8-2*C9)**3)/12-4*(P17+P13*S15**2)+4*(P18+P14*S16**2))
        E68 = (C34*I13*10**4*np.pi**2/(C32*1000)**2)/1000

        E67 = np.sqrt(((C13*10**2*C30)/1000)/E68)
        D27 = D26
        F27 = wenn(D27=="a0",0.13,wenn(D27=="a",0.21,wenn(D27=="b",0.34,wenn(D27=="c",0.49,0.76))))
        E66 = 0.5*(1+F27*(E67-0.2)+E67**2)
        E65 = np.min([1/(E66+np.sqrt(E66**2-E67**2)),1])
        I17 = 1/(10**3)*((C7*C8**2)/4-((C7-2*C9)*(C8-2*C9)**2)/4-4*(P13*S15)+4*(P14*S16))
        I15 = 2*I13/C8*10
        F17 = 1/(10**3)*((C8*C7**2)/4-((C8-2*C9)*(C7-2*C9)**2)/4-4*(P13*P15)+4*(P14*P16))
        F15 = 2*F13/C7*10
        D100 = (C13*10**2*C30)/1000
        C87 = G31/(E57*D100/C31)

        D102 = (wenn(C35<3,I17,I15)*10**3*C30)/(1000**2)
        D114 = np.min([H112*(1+(E67-0.2)*G31/(E65*D100/C31)),H112*(1+0.8*G31/(E65*D100/C31))])
        D112 = 0.6 * D114
        D111 = np.min([C38*(1+(E59-0.2)*G31/(E57*D100/C31)),C38*(1+0.8*G31/(E57*D100/C31))])
        E111 = np.min([C38*(1+0.6*E59*G31/(E57*D100/C31)),C38*(1+0.6*G31/(E57*D100/C31))])
        E112 = np.min([H112*(1+0.6*E67*G31/(E65*D100/C31)),H112*(1+0.6*G31/(E65*D100/C31))])
        G87 = wenn(C35<3,D112,E112)*(G33/(D102/C31))
        D101 = (wenn(C35<3,F17,F15)*10**3*C30)/(1000**2)
        E87 = wenn(C35<3,D111,E111)*(G32/(D106*D101/C31))
        C87 = G31/(E57*D100/C31)

        E113 = wenn((G31 != 0 and G32 != 0 and G33==0),0,0.8*E111)
        E114 = np.min([H112*(1+0.6*E67*G31/(E65*D100/C31)),H112*(1+0.6*G31/(E65*D100/C31))])
        D113 = wenn((G31 != 0 and G32 != 0 and G33==0),0,0.6*D111)
        E90 = wenn(C35<3,D113,E113)*(G32/(D106*D101/C31))
        C90 = G31/(E65*D100/C31)
        G90 = wenn(C35<3,D114,E114)*(G33/(D102/C31))

        I87 = C87 + E87 + G87 

        I90 = C90+E90+G90

    except: 
        print("Error occured while computing static steel Profile Hohl beam!")

        raise
    return I87,I90
    
def compute_static_steel_RH_beam(cls_type, y_M1, f_y, E, L_Cry, # Information about the element
                                 N_Ed, M_zEd, M_yEd): # Information about loads
    """
    Wrapper function to compute the static loads for a steel beam of type RH.
    
    Parameters:
        cls_type (str): Type of steel beam
        y_M1 (float): 
        f_y (float): 
        E (float): 
        L_Cry (float): 
        N_Ed (float): Load
        M_zEd (float): Load
        M_yEd (float): Load
    
    Returns:
        bool: Whether the element can withstand the loads.
    """

    qs_info = get_steel_RH_QS_info(cls_type)
    qkl = compute_QKL_steel_RH(f_y, r_i=float(qs_info['ri']), t=float(qs_info['T']), h=float(qs_info['H']), b=float(qs_info['B']))
    M_N_1, M_N_2 = only_static_steel_RH_beam(cls_type, qkl, y_M1, f_y, E, L_Cry, N_Ed, M_zEd, M_yEd)

    if M_N_1 > 1:
        return False
    if M_N_2 > 1: 
        return False
    return True

'''
def only_static_steel_O_beam():
    pass # no calculations for steel O beam available

def compute_static_steel_O_beam():
    pass # no calculations for steel O beam available
'''
