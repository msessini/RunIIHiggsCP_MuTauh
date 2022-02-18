algos_translator = {
    "MTT": "FastMTT",
    "SVFit": "SVFit",
    "Mixed": "Combined",
    "GEF": "GEF"
}

all_algos = ['MTT','SVFit','Mixed','GEF']

def build_algos(channel):
    algos = all_algos[:]
    if channel == 'a1a1':
        algos = all_algos[2:]
    if channel == 'rhorho':
        algos = all_algos[:3]
    if channel == 'pipi':
        algos = all_algos[:3]
    if channel == 'a1pi':
        algos = all_algos[2:]
    if channel == 'a1rho':
        algos = all_algos[:3]
    if channel == 'pirho':
        algos = all_algos[:3]
    if channel == 'a1mu':
        algos = all_algos[2:]
    if channel == 'rhomu':
        algos = all_algos[:3]
    if channel == 'pimu':
        algos = all_algos[:3]
    return algos
