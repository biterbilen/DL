import inspect

myself = lambda: inspect.stack()[1][3]

chrstrip = lambda x: int(x.strip('chr').replace('X','23'))

flatten = lambda l: [item for sublist in l for item in sublist]

nuc2ord = lambda s: {'A':[1,0,0,0], 'C':[0,1,0,0], 'G':[0,0,1,0], 'T':[0,0,0,1]}[s]

label2int = lambda l: {'B':1, 'U':0}[l]

listComplement = lambda l,lx: [ i for i in l if not i in lx ]