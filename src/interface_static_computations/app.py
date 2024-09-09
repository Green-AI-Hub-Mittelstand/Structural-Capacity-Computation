from flask import Flask, render_template, request
import numpy as np
import pandas as pd

import sys
sys.path.append("src/")

from static_computation.utils_steel import *
from static_computation.utils_reinforced_concrete import *


app = Flask(__name__)
global parameters
parameters = {}
parameters['wall'] = {
        "concrete_cls": ["Concrete class, e.g., 'C20/25'","","str"],    
        "h": ["Height","mm","float"], 
        "b": ["Width","mm","float"],
        "l": ["Length", "mm","float"],
        "As_rechts": ["As (outside)","mm2/m","float"],
        "As_links": ["As (inside)","mm2/m","float"],
        "d1_links": ["d1 (inside) ","mm","float"], 
        "d1_rechts":["d1 (outside)","mm","float"], 
        "console": ["Whether a console is present","","bool"],
        "corbel_notch_av": ["A_v. If console is False this value will not be used and can be set to any value.","mm","float"],
        "N_Ed": ["NEd", "kN","float"],
        #"phi_eff": ["Effektiver Kriechbeiwert, default is 1.","","float"]
    }

parameters['slabs'] = {
        "concrete_cls": ["Concrete class, e.g., 'C20/25'","","str"],     
        "h": ["Height","mm","float"], 
        "b": ["Width","mm","float"],
        "As1_unten": ["As (bottom) ","mm2/m","float"],
        "As2_oben": ["As (top)","mm2/m","float"],
        "d1_unten": ["d1 (inside) ","mm","float"], 
        "d1_oben":["d1 (outside)","mm","float"], 
        "M_Ed": ["MEd", "kNm","float"],
        "V_Ed": ["VEd", "kNm","float"]
}
parameters['column_RC'] = {
        "concrete_cls": ["Concrete class, e.g., 'C20/25'","","str"],    
        "strength1": ["Steel strength", "N/mm2","float"], 
        "L_col": ["Height", "mm","float"],
        "h": ["rc_height ","mm","float"], 
        "b": ["rc_width","mm","float"],
        "As_oben": ["As (top)","mm2/m","float"],
        "As_unten": ["As (bottom)","mm2/m","float"],
        "As_links": ["As (left)","mm2","float"],
        "As_rechts": ["As (right)","mm2","float"],       
        "d1_oben": ["d1 (top) ","mm","float"], 
        "d1_unten":["d1 (bottom)","mm","float"], 
        "console": ["Whether the column has a console" ,"","bool"],
        "corbel_notch_av": ["av of console. If no console exists, input any number","mm","float"],
        "N_Ed": ["NEd", "kN","float"],
        #"phi_eff": ["Effektiver Kriechbeiwert, default is 1.","","float"]
    }

parameters['column_I'] = {
        "cls_type": ["Steel class type, e.g. “HEA 160”","","str"],
        "strength1": ["Steel strength", "N/mm2","float"],  
        #"y_M1": ["","","float"],
        "L": ["Height","mm","float"],
        "E": ["Elasticity modulus","N/mm2","float"],
        #"v": ["","","float"], 
        #"L_T":["height","mm","float"], 
        "N_Ed": ["NEd", "kN","float"],

}


parameters['column_O'] = {
        "cls_type": ["Steel class type, e.g. “219,1x5 (w)","","str"],   
        "strength1": ["Steel strength", "N/mm2","float"],   
        #"y_M1": ["","","float"],
        "L_Cr": ["Height","mm","float"],
        "E": ["Elasticity modulus","N/mm2","float"],
        "N_Ed": ["NEd", "kN","float"],

}

parameters['column_RH'] = {
        "cls_type": ["Steel class type, e.g. “200x100x5 (w)” ","","str"],   
        "strength1": ["Steel strength", "N/mm2","float"],  
        #"y_M1": ["","","float"],
        "E": ["Elasticity modulus","N/mm2","float"],
        "L_Cry": ["Height","mm","float"],
        "N_Ed": ["NEd", "kN","float"],

}

parameters['beam_RC'] = {
        "concrete_cls": ["Concrete strength class, e.g., e.g., 'C20/25'","","str"],   
        "strength1": ["Steel strength", "N/mm2","float"],   
        "h": ["h ","mm","float"], 
        "b": ["b","mm","float"],        
        "As1_unten": ["As (bottom)","mm2/m","float"],
        "As2_oben": ["As (top)","mm2/m","float"],     
        "d1_unten":["d1 (bottom)","mm","float"], 
        "d1_oben": ["d1 (top) ","mm","float"], 
        "asw_z": ["asw_z ","mm","float"], 
        "M_Ed": ["MEd", "kNm","float"],
        #"V_Ed": ["VEd","kNm","float"]
    }


parameters['beam_I'] = {
        "cls_type": ["Steel class type, e.g. “HEA 160”","","str"],   
        "strength1": ["Steel strength", "N/mm2","float"],  
        #"y_M1": ["","","float"],
        "E": ["Elasticity modulus","N/mm2","float"],
        "L": ["Length","mm","float"],
        #"z_p": ["","","float"], 
        #"N_Ed": ["compressive force NEd", "kN","float"],
        "M_Ed": ["MEd", "kNm","float"],
        "M_yEd": ["MyEd", "kNm","float"],

}

parameters['beam_RH'] = {
        "cls_type": ["Steel class type, e.g.“200x100x5 (w)”","","str"],   
        "strength1": ["Steel strength", "N/mm2","float"],  
        #"f_y": ["","","float"], 
        "E": ["Elasticity modulus","N/mm2","float"],
        "L_Cry": ["Length","mm","float"],
        #"N_Ed": ["compressive force NEd", "kN","float"],
        "M_zEd": ["MEd", "kNm","float"],
        "M_yEd": ["MyEd", "kNm","float"],

}



@app.route('/', methods=['GET', 'POST'])
def index():    
    return render_template('index_interface_static.html',category = 'start')


@app.route('/get_parameters', methods=['GET', 'POST'])
def get_parameters():
    category = request.args.get('category')
    material = ""        
    if category == 'wall' or category == 'slabs' or category == 'column_RC' or category == 'beam_RC':
        material = "reinforced concrete"
    else:
        material = "steel"
    parameters_category = parameters[category]
    return render_template('index_interface_static.html',parameters = parameters_category,category=category,material=material)

@app.route('/calculate', methods=['GET', 'POST'])
def calculate():
    try: 
        input_parameters = {}
        category = request.args.get('category')
        material = request.args.get('material')
        parameters_category = parameters[category]
        for p in parameters_category.keys():
            if parameters_category[p][2] == "float":
                input_parameters[p] = float(request.form.get(p))
            elif parameters_category[p][2] == "bool":
                input_parameters[p] = request.form.get(p) == 'on'
            else:
                input_parameters[p] = request.form.get(p)
    
        if category == 'wall':
            result = compute_static_stb_wall(*input_parameters.values(),phi_eff=1)
        elif category == 'slabs':
            result = compute_static_stb_slab(*input_parameters.values())
        elif category == 'column_RC':
            result = compute_static_stb_column(*input_parameters.values(),phi_eff=1)
        elif category == 'column_I':                                                
            result = compute_static_steel_I_column(cls_type=input_parameters["cls_type"],
                                                   f_y=input_parameters["strength1"], 
                                                   y_M1=1.1, 
                                                   L=input_parameters["L"],
                                                   E=input_parameters["E"], 
                                                   v=0.3, 
                                                   L_T=input_parameters["L"],
                                                   N_Ed=input_parameters["N_Ed"])
    
        elif category == 'column_O': 
            result = compute_static_steel_O_column(cls_type=input_parameters["cls_type"],
                                                    f_y=input_parameters["strength1"], 
                                                    y_M1=1.1,
                                                    L_Cr=input_parameters["L_Cr"],
                                                    E=input_parameters["E"],
                                                    N_Ed=input_parameters["N_Ed"])
        elif category == 'column_RH':
            

            result = compute_static_steel_RH_column(cls_type=input_parameters["cls_type"],
                                                    f_y=input_parameters["strength1"],
                                                    y_M1=1.1,
                                                    E=input_parameters["E"],
                                                    L_Cry=input_parameters["L_Cry"],
                                                    N_Ed=input_parameters["N_Ed"])
        elif category == 'beam_RC':
            result = compute_static_stb_beam(*input_parameters.values(), V_Ed=0)
        elif category == 'beam_I':
            
            result = compute_static_steel_I_beam(cls_type=input_parameters["cls_type"], 
                                                f_y=input_parameters["strength1"], 
                                                y_M1=1.1, 
                                                L=input_parameters["L"], 
                                                E=input_parameters["E"], 
                                                #v=0.3, 
                                                z_p=-80, # in mm
                                                N_Ed=0,
                                                M_Ed=input_parameters["M_Ed"],
                                                M_yEd=input_parameters["M_yEd"]
                                                )
        elif category == 'beam_RH':
        
            print(input_parameters)
            result = compute_static_steel_RH_beam(cls_type=input_parameters["cls_type"],
                                                f_y=input_parameters["strength1"],
                                                y_M1=1.1,
                                                E=input_parameters["E"],
                                                L_Cry=input_parameters["L_Cry"],
                                                N_Ed=0,
                                                M_zEd=input_parameters["M_zEd"], 
                                                M_yEd=input_parameters["M_yEd"])
        else:
            raise ValueError(f"Category has no valid value!")   
      
        return render_template('index_interface_static.html',parameters = parameters_category,category=category,material=material, result= result, values=input_parameters)
    except Exception as e:
        error_message = str(e)
        raise e
        return render_template('index_interface_static.html',parameters = parameters_category,category=category, material=material, error_message = error_message, values=input_parameters)



if __name__ == '__main__':
    app.run(debug=False)