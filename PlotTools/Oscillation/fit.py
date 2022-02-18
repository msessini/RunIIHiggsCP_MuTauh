from ROOT import TF1, TMath

fit_style = [
    {'style': 1, 'color': 1, 'width': 1},
    {'style': 9, 'color': 8, 'width': 1},
    {'style': 5, 'color': 9, 'width': 1},
    {'style': 6, 'color': 13, 'width': 1}
]

fit_function = TF1('fit','[0]*cos(x)+[1]',0,2*TMath.Pi())
