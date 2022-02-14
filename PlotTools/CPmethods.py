methods_translator = {
    "polarimetric": "PV + PV",
    "polarimetricGEF": "PV + PV with GEF",
    "decayplane": "NP + NP",
    "impactparameter": "IP + IP",
    "DPIP": "DP + IP",
    "PVIP": "PV + IP"
}

all_methods = ['polarimetric','polarimetricGEF','decayplane','impactparameter','DPIP','PVIP']

def build_methods(channel,genlevel=False):
    methods = all_methods[:]
    if channel == 'a1a1':
        methods = all_methods[:3]
        if(genlevel==True):
            methods = all_methods[:3:2]
    if channel == 'rhorho':
        methods = all_methods[:3:2]
    if channel == 'pipi':
        methods = all_methods[:4:3]
    if channel == 'a1pi':
        methods = all_methods[:2]+all_methods[4:]
        if(genlevel==True):
            methods = all_methods[:1]+all_methods[4:]
    if channel == 'a1rho':
        methods = all_methods[:3:2]
    if channel == 'pirho':
        methods = all_methods[:1]+all_methods[4:]
    if channel == 'a1mu':
        methods = all_methods[4:]
    if channel == 'rhomu':
        methods = all_methods[4:]
    if channel == 'pimu':
        methods = all_methods[0]
    return methods
