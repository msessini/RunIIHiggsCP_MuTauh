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
    if channel == 'A1A1':
        methods = all_methods[:3]
        if(genlevel==True):
            methods = all_methods[:3:2]
    if channel == 'RHORHO':
        methods = all_methods[:3:2]
    if channel == 'PIONPION':
        methods = all_methods[:4:3]
    if channel == 'A1PION':
        methods = all_methods[:2]+all_methods[4:]
        if(genlevel==True):
            methods = all_methods[:1]+all_methods[4:]
    if channel == 'A1RHO':
        methods = all_methods[:3:2]
    if channel == 'RHOPION':
        methods = all_methods[:1]+all_methods[4:]
    if channel == 'A1MU':
        methods = all_methods[4:]
    if channel == 'RHOMU':
        methods = all_methods[4:]
    if channel == 'PIONMU':
        methods = all_methods[0]
    return methods
