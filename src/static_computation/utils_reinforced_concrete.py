import numpy as np
import pandas as pd

########################
##### Table Helper #####
########################

# Additional information for the usage of the table w_table_w
w_table_w_col_dict = { "d_2/d=0.05": {"w1_id": "1",
                                       "w2_id": "2",
                                        "sigma_sd": -435.7},
                        "d_2/d=0.1": {"w1_id": "3",
                                      "w2_id": "4",
                                      "sigma_sd": -435.3},
                        "d_2/d=0.15": {"w1_id": "5",
                                       "w2_id": "6",
                                       "sigma_sd": -434.9},
                        "d_2/d=0.2": {"w1_id": "7",
                                      "w2_id": "8",
                                      "sigma_sd": -388.9}}

def get_values_w_table_w(d2_d, w1=None, w2=None):
    """
    Function to retrieve data from the table w_table_w which are necessary for the computation of static loads. 
    For general information check documentations about static loads.

    Parameters:
        d2_d (float):
        w1 (float, optional): Defaults to None. Either w1 or w2 have to contain a value.
        w2 (float, optional): Defaults to None. Either w1 or w2 have to contain a value.

    Returns:
        value (float): Either w1 or w1 depending which input paramter was none.
        u_Eds (float):
        sigma_sd2 (float):
    """
    # Get table
    w_table_w = pd.read_csv("data/static_computation/reinforced_concrete/w_table_w.csv", sep=";")

    # Identify d2_d and get columns + sigma_sd2
    if d2_d <= 0.2:
        cat = 0.2
        if d2_d <= 0.15:
            cat = 0.15
            if d2_d <= 0.1:
                cat = 0.1
                if d2_d <= 0.05:
                    cat = 0.05
    else:
        raise ValueError(f"The paramter d_2/d must be smaller or equal 0.2, not {d2_d}!")
    w1_id, w2_id, sigma_sd2 = w_table_w_col_dict[f"d_2/d={cat}"].values()
    
    # If w1 was provided, find w2
    if w1:
        row = w_table_w[w_table_w[w1_id] >= w1]
        if not row.empty:
            row = row.iloc[0]
        else:
            row = w_table_w.iloc[-1]
            print(f"Info: get_values_w_table_w(d2_d={d2_d}, w1={w1}, w2={w2}), but maximal w1 is {w_table_w[w1_id]}")
        value = row[w2_id]
        u_Eds = row["u_Eds"]
    # If w2 was provided, find w1
    elif w2:
        row = w_table_w[w_table_w[w2_id] >= w2]
        if not row.empty:
            row = row.iloc[0]
        else:
            row = w_table_w.iloc[-1]
            print(f"Info: get_values_w_table_w(d2_d={d2_d}, w1={w1}, w2={w2}), but maximal w2 is {w_table_w.w2}")

        value = row[w1_id]
        u_Eds = row["u_Eds"]
    # If none were provided, raise error
    else:
        raise ValueError("Both w1 and w2 are none. One of them has to be set!")
    
    return value, u_Eds, sigma_sd2

def get_values_w_table_wo(w1=None, w1_sigma_sd=None):
    """
    Function to retrieve data from the table w_table_wo which are necessary for the computation of static loads. 
    For general information check documentations about static loads.

    Parameters:
        w1 (float, optional): Defaults to None. Either w1 or w1_sigma_sd have to contain a value.
        w1_sigma_sd (float, optional): Defaults to None. Either w1 or w1_sigma_sd have to contain a value.

    Returns:
        value (float): Either w1 or w1_sigma_sd depending which input paramter was none.
        u_Eds (float):
    """
    # Get table
    w_table_wo = pd.read_csv("data/static_computation/reinforced_concrete/w_table_wo.csv", sep=";")
    
    # If w1 was provided, find w1_sigma_sd
    if w1:
        row = w_table_wo[w_table_wo["w1"] >= w1]
        if not row.empty:
            row = row.iloc[0]
        else:
            row = w_table_wo.iloc[-1]
            print(f"Info: get_values_w_table_wo(w1={w1}, w1_sigma_sd={w1_sigma_sd}), but maximal w1 is {w_table_wo['w1']}")

        return row["w_1/sigma_sd"], row["u_Eds"]
    
    # If w1_sigma_sd was provided, find w1
    elif w1_sigma_sd:
        row = w_table_wo[w_table_wo["w_1/sigma_sd"] >= w1_sigma_sd]
        if not row.empty:
            row = row.iloc[0]
        else:
            row = w_table_wo.iloc[-1]
            print(f"Info: get_values_w_table_wo(w1={w1}, w1_sigma_sd={w1_sigma_sd}), but maximal w1_sigma_sd is {w_table_wo['w_1/sigma_sd']}")
        return row["w1"], row["u_Eds"]
    
    # If none was provided, raise an error
    else:
       raise ValueError("Both w1 and w_1/sigma_sd are none. One of them has to be set!") 

def get_values_concrete_cls(concrete_cls):
    """
    Get values from table with information regarding concrete classes

    Parameters:
        concrete_cls (str): concrete class, e.g., 'C20/25'
    Returns:
        "f_cd" (float):
        "v1*f_cd" (float):
    """
    fcd_table = pd.read_csv("data/static_computation/reinforced_concrete/concrete_cls_fcd.csv", sep=";")
    if concrete_cls in fcd_table[fcd_table["v1*f_cd"].notna()].concrete_cls.values:
        row = fcd_table[fcd_table["concrete_cls"] == concrete_cls].iloc[0]
        return row["f_cd"], row["v1*f_cd"]
    else:
        raise LookupError(f"Concrete information are not implmeneted for {concrete_cls}. Use one of the following: {fcd_table[fcd_table['v1*f_cd'].notna()].concrete_cls.to_list()}.")

# Function to get f_ck from concrete class name
def get_fck(concrete_cls):
    """
    Get f_ck from the concrete class. As the class alreday contains the variable, the string wil be simply split
    
    Parameters:
        concrete_cls (str): concrete class, e.g., 'C20/25'
    
    Returns:
        "f_ck" (float): 
    """
    return int(concrete_cls[1:3])

##############################
##### Static Computation #####
##############################
    
def compute_static_stb_beam(concrete_cls, strength1, h, b, As1_unten, As2_oben, d1_unten, d1_oben, asw_z, # Information about the element
                            M_Ed, V_Ed): # Information about loads
    """
    Compute static strength of a rc beam.
    
    Parameters:
        concrete_cls (str): Concrete class, e.g., 'C20/25'
        strength1 (float): Strength of steel material
        h (float): 
        b (float): 
        As1_unten (float): 
        As2_oben (float): 
        d1_unten (float): 
        d1_oben (float): 
        asw_z (float): 
        M_Ed (float): Load
        V_Ed (float): Load
    
    Returns:
        bool: Whether the beam can withstand the loads without failure
    """
    
    #-----Biegenachweis-----
    if d1_unten:
        d = h - d1_unten
    else:
        d = 0.9 * h
    f_cd, _ = get_values_concrete_cls(concrete_cls)

    w1_sigma_sd = As1_unten / (b * d * f_cd)
    w1, _ = get_values_w_table_wo(w1_sigma_sd=w1_sigma_sd)

    if w1 < 0.3:
        _, u_Eds = get_values_w_table_wo(w1_sigma_sd=w1_sigma_sd)
        M_Eds = u_Eds * b * d**2 * f_cd
        N_Ed = 0
    else:
        delta_sd = 436.8 # N/mm2
        w1 = As1_unten * delta_sd / (b * d * f_cd)
        d2 = d1_oben
        d2_d = d2 / d
        w2, u_Eds, sigma_sd2 = get_values_w_table_w(d2_d=d2_d, w1=w1)

        w2_option = As2_oben * sigma_sd2 / (b * d * f_cd)
        if w2_option >= w2:
            pass
        else:
            w2 = w2_option
            _, u_Eds, _ = get_values_w_table_w(d2_d=d2_d, w2=w2)

        M_Eds = u_Eds * b * d**2 * f_cd
        N_Ed = 0
    M_Rd = M_Eds / 1000000 # Newtonmillimeter into Kilonewtonmeter
    #---Nachweis---
    if M_Ed/M_Rd > 1:
        return False

    #-----Schubnachweis----- #
    cot = 1.2 # Maßeinheit
    tan = 1/1.2 # Maßeinheit
    z = 0.9 * d # mm
    #---NW der Betondruckstrebe---  
    _, v1_pl = get_values_concrete_cls(concrete_cls) # MN/m2 bzw. N/mm2
    V_Rdmax = v1_pl * b * z * (1 / (tan + cot))
    #---NW der Schubbewehrung---
    a_sw = asw_z / 100 # cm2/m? in mm2/m
    f_ywd = strength1 # N/mm2
    V_Rdsw = a_sw / (f_ywd * z * cot)
    V_Rd = np.min([V_Rdsw, V_Rdmax]) / 1000000 # Newtonmillimeter into Kilonewtonmeter

    #---Nachweis---
    if V_Ed/V_Rd > 1:
        return False
    return True

def compute_static_stb_slab(concrete_cls, h, b, As1_unten, As2_oben, d1_unten, d1_oben, # Information about the element
                            M_Ed, V_Ed): # Information about loads
    """
    Compute static strength of a rc slab. 
    
    Parameters:
        concrete_cls (str): Concrete class, e.g., 'C20/25'
        h (float): 
        b (float): 
        As1_unten (float): 
        As2_oben (float): 
        d1_unten (float): 
        d1_oben (float): 
        M_Ed (float): Load
        V_Ed (float): Load
    
    Returns:
        bool: Whether the slab can withstand the loads without failure
    """

    #-----Biegenachweis-----
    if d1_unten:
        d = h - d1_unten
    else:
        d = 0.9 * h
    f_cd, _ = get_values_concrete_cls(concrete_cls)

    w1_sigma_sd = As1_unten / (b * d * f_cd)
    w1, _ = get_values_w_table_wo(w1_sigma_sd=w1_sigma_sd)

    if w1 < 0.3:
        _, u_Eds = get_values_w_table_wo(w1_sigma_sd=w1_sigma_sd)
        M_Eds = u_Eds * b * d**2 * f_cd
        N_Ed = 0
    else:
        w1 = As1_unten * 436.8 / (b * d * f_cd)
        d2 = d1_oben
        d2_d = d2 / d
        w2, u_Eds, sigma_sd2 = get_values_w_table_w(d2_d=d2_d, w1=w1)

        w2_option = As2_oben * sigma_sd2 / (b * d * f_cd)
        if w2_option >= w2:
            pass
        else:
            w2 = w2_option
            _, u_Eds, _ = get_values_w_table_w(d2_d=d2_d, w2=w2)

        M_Eds = u_Eds * b * d**2 * f_cd

    N_Ed = 0
    M_Rd = M_Eds / 1000000 # Newtonmillimeter into Kilonewtonmeter

    #---Nachweis---
    if M_Ed/M_Rd > 1:
        return False

    #-----Schubnachweis-----
    C_Rdc = 0.15 / 1.5
    b_w = b
    f_ck = get_fck(concrete_cls)
    A_sl = As1_unten
    k = 1 + np.sqrt(200 / d) # <= 2.0 mit d in mm 
    rho_1 = A_sl / (b_w * d) # <= 0.02
    gamma_c = 1.5 # bis C50/60
    if d <= 600:
        V_min = (0.0525 / gamma_c) * k**(3/2) * f_ck**(1/2)
    elif d > 800:
        V_min = (0.0375 / gamma_c) * k**(3/2) * f_ck**(1/2)
    else:
        V_min1 = (0.0525 / gamma_c) * k**(3/2) * f_ck**(1/2)
        V_min2 = (0.0375 / gamma_c) * k**(3/2) * f_ck**(1/2)
        np.interp(d, [600, 800], [V_min1, V_min2])
    V_Rdc = np.min([(C_Rdc * k * (100 * rho_1 * f_ck)**(1/3)) * b_w * d,
                    (V_min * b_w * d)])
    
    V_Rdc = V_Rdc / 1000000 # Newtonmillimeter into Kilonewtonmeter
    #---Nachweis---
    if V_Ed/V_Rdc > 1:
        return False
    return True

def compute_static_stb_wall(concrete_cls, h, b, l, As_rechts, As_links, d1_links, d1_rechts, # Information about the element 
                            console, corbel_notch_av, # Information about a potential console
                            N_Ed, # Information about loads
                            phi_eff=1): # "Effektiver Kriechbeiwert"
    """
    Compute static strength of a rc wall.
    
    Parameters:
        concrete_cls (str): Concrete class, e.g., 'C20/25'
        strength1 (float): Strength of steel material
        l (float): 
        h (float): 
        b (float): 
        As_rechts (float): 
        As_links (float): 
        d1_links (float): 
        d1_rechts (float): 
        console (bool): Wheter a console is present
        corbel_notch_av (float): A_v. If console is False this value will not be used and can be set to any value.
        N_Ed (float): Load
        phi_eff (float): Effektiver Kriechbeiwert, default is 1. 
    
    Returns:
        bool: Whether the wall can withstand the loads without failure
    """
    h = h /1000 # turn mm into m
    b = b /1000 # turn mm into m
    l = l /1000 # turn mm into m
    As_rechts = As_rechts /1000**2 # turn mm^2 into m^2
    As_links = As_links /1000**2 # turn mm^2 into m^2
    d1_rechts = d1_rechts /1000 # turn mm into m
    d1_links = d1_links /1000 # turn mm into m
    corbel_notch_av = corbel_notch_av /1000 # turn mm into m
    
    #-----Knicknachweis-----
    beta = 1
    L_col = l 
    L_0 = beta * L_col
    i = 0.289 * h
    lambd = L_0 / i
    d = 0.9 * h

    if console:
        e_0 = (b + corbel_notch_av) * 0.5
    else:
        e_0 = np.max([2, b/30])
    f_ck = get_fck(concrete_cls)
    beta_phi = 0.35 + f_ck/200 - lambd/150
    K_phi = 1 + beta_phi * phi_eff

    K_1 = np.max([0, np.min([(0.1 * lambd - 2.5), 1])])
    K_2 = 1
    epsilon_yd = 0.0217 / 100 # percent per mille
    G40 = epsilon_yd
    G41 = 2 * K_2 * G40 / (0.9 * d)
    e_2 = K_phi * K_1 * K_2 * 0.1 * L_0**2 * G41

    theta_i = np.min([(1 / (100 * np.sqrt(L_col))), 1 / 200])
    e_a = theta_i * L_0 / 2

    e_tot = e_0 + e_a + e_2

    e = K_phi * e_tot
    M_Ed = e * N_Ed

    #-----M-N-Interaktion-----
    #---(1) maximal zulässige Druckkraft (wird negativ angegeben)---
    as_innen = As_rechts
    as_aussen = As_links
    as_ = np.min([as_innen, as_aussen])
    f_cd, _ = get_values_concrete_cls(concrete_cls)
    sigma_s = 435 * 1000 # N/mm2 in kN/m2
    a_c = b * 100 - (as_innen + as_innen)
    n_max = -(a_c*f_cd+as_ * sigma_s + as_ * sigma_s)

    #---(3) maximal zulässige Zugkraft (wird positiv angegeben)---
    n_min = as_ * sigma_s + as_ * sigma_s 

    #---(2) maximal zulässiges Moment & Druckkraft die gleichzeitig aufgenommen werden kann---
    d_1 = np.max([d1_links, d1_rechts])
    x = (h - d_1) * 3.5 / (3.5 + 2.17)
    sigma_s1 = 435 * 1000 # N/mm2 in kN/m2
    sigma_s2 = 435 * 1000 # N/mm2 in kN/m2
    F_c = 0.8 * x * f_cd
    F_s1 = as_ * sigma_s1
    F_s2 = as_ * sigma_s2
    n_Rd = -F_c - F_s1 + F_s2
    m_Rd = F_c * (x - d_1) + F_s1 * (h / 2 - d_1) + F_s2 * (h / 2 - d_1)

    #-----Nachweise-----
    if n_max > N_Ed:
        return False
    if n_min < N_Ed:
        return False
    # check demand.nEd>=g(M)
    g = n_max + ((n_Rd - n_max) / m_Rd) * M_Ed
    if N_Ed < g:
        return False
    # check demand.nEd<=f(M)
    f = n_min + ((n_Rd - n_min) / m_Rd) * M_Ed
    if N_Ed > f:
        return False
    return True

def compute_static_stb_column(concrete_cls, strength1, L_col, h, b, As_oben, As_unten, As_links, As_rechts, d1_oben, d1_unten, # Information about the element
                              console, corbel_notch_av, # Information about a potential console
                              N_Ed, # Information about loads
                              phi_eff=1): # "Effektiver Kriechbeiwert" 
    """
    Compute static strength of a rc beam.
    
    Parameters:
        concrete_cls (str): Concrete class, e.g., 'C20/25'
        strength1 (float): Strength of steel material
        L_col (float): 
        h (float): 
        b (float): 
        As_unten (float): 
        As_oben (float): 
        d1_unten (float): 
        d1_oben (float): 
        console (bool): Whether a console is present.
        corbel_notch_av (float): A_v. If console is False this value will not be used and can be set to any value.
        N_Ed (float): Load
        phi_eff (float): Effektiver Kriechbeiwert, default is 1. 
    
    Returns:
        bool: Whether the column can withstand the loads without failure
    """
    L_col = L_col/1000 # turn mm into m
    h = h/1000 # turn mm into m
    b = b/1000 # turn mm into m
    As_oben = As_oben/1000000 # turn mm2 into m2
    As_unten = As_unten/1000000 # turn mm2 into m2
    As_links = As_links/1000000 # turn mm2 into m2
    As_rechts = As_rechts/1000000 # turn mm2 into m2
    d1_oben = d1_oben/1000 # turn mm into m
    d1_unten = d1_unten/1000 # turn mm into m
    corbel_notch_av = corbel_notch_av/1000 # turn mm into m
    N_Ed = -N_Ed

    #-----Knicknachweis-----
    #---(1) maximal zulässige Druckkraft (wird negativ angegeben)---
    A_c = b * h - (As_oben + As_unten + As_links + As_rechts)
    A_s = As_oben + As_unten + As_links + As_rechts
    f_cd, _ = get_values_concrete_cls(concrete_cls)
    sigma_s = strength1

    N_max = -(A_c * f_cd + A_s * sigma_s + A_s * sigma_s)

    #---(3) Maximal zulässige Zugkraft (Zug wird positiv angegeben)---
    N_min = A_s * sigma_s 

    #---(2) Bestimmung des Maximal zulässiges Moment und Druckkraft die gleichzeitig aufgenommen werden kann in Abhängigkeit der Knickrichtung---
    beta = 1
    L_0 = beta * L_col
    i = 0.289 * h
    lambd = L_0 / i
    d = 0.9 * h

    if console:
        e_0 = (b + corbel_notch_av) * 0.5
    else:
        e_0 = np.max([0.02, b/30])

    f_ck = get_fck(concrete_cls)
    beta_phi = 0.35 + f_ck/200 - lambd/150
    K_phi = 1 + beta_phi * phi_eff

    K_1 = np.max([0, np.min([(0.1 * lambd - 2.5), 1])])
    K_2 = 1
    epsilon_yd = 0.00217 # percent per mille
    G40 = epsilon_yd
    G41 = 2 * K_2 * G40 / (0.9 * d)
    e_2 = K_phi * K_1 * K_2 * 0.1 * L_0**2 * G41

    theta_i = np.min([(1 / (100 * np.sqrt(L_col))), 1 / 200])
    e_a = theta_i * L_0 / 2

    e_tot = e_0 + e_a + e_2

    e = K_phi * e_tot
    M_Ed = e * N_Ed

    d_1 = np.max([d1_oben, d1_unten])
    x = (h - d_1) * 3.5 / (3.5 + 2.17)
    sigma_s1 = strength1
    sigma_s2 = strength1
    F_c = 0.8 * x * f_cd * b
    As = np.min([As_oben, As_unten])
    F_s1 = As * sigma_s1
    F_s2 = As * sigma_s2
    N_Rd = -F_c - F_s1 + F_s2
    M_Rd = F_c * (x - d_1) + F_s1 * (h / 2 - d_1) + F_s2 * (h / 2 - d_1)

    #-----Nachweise-----
    g = N_max + ((N_Rd - N_max) / M_Rd) * M_Ed
    f = N_min + ((N_Rd - N_min) / M_Rd) * M_Ed
    # check demand.nEd>=g(M)
    if N_Ed < g:
        return False
    # check demand.nEd<=f(M)
    if N_Ed > f:
        return False
    
    return True
