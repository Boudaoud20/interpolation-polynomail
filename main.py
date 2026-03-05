import argparse
from ex1 import Lagrange
from ex2 import Newton
import numpy as np

def main():

    parser = argparse.ArgumentParser(description="Lagrange and Newton method to extract polynomial form")
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-l","--Lagrange",help="Lagrange method",action="store_true")
    group.add_argument("-n","--Newton",help="Newton method",action="store_true")

    parser.add_argument("-f",help="using function",action="store_true")
    parser.add_argument("-x",help="x coordinate",type=float,nargs="+",required=True)
    parser.add_argument("-y",help="y coordinate",type=float,nargs="+")
    parser.add_argument("-i",help="interval",type=float,nargs="+",required=True)

    arg=parser.parse_args()
    if arg.f:
        if arg.y is not None:
            parser.error("When using -f you must NOT provide -y")
        elif len(arg.i) !=2 or arg.i[0] > arg.i[1]:
            parser.error("you must write the correct  interval")
        else:
            points = tuple(arg.x)
            if arg.Lagrange:
                lagrange = Lagrange(points=points,use_function=True,lower=arg.i[0], upper=arg.i[1])
            else:
                newton = Newton(points=np.array(arg.x),use_function=True,lower=arg.i[0], upper=arg.i[1])
        
    else:
        if arg.y is None:
            parser.error("you need to provide y coordinate in this case")
        elif len(arg.i) !=2 or arg.i[0] > arg.i[1]:
            parser.error("you must write the correct interval")
        elif len(arg.y) != len(arg.x):
            parser.error("the length of x and y coordinate must be equal")
        else:
            points = dict()
            for i in range(len(arg.y)):
                points[arg.x[i]]=arg.y[i]
            if arg.Lagrange:
                lagrange = Lagrange(points=points,use_function=False,lower=arg.i[0], upper=arg.i[1] )
            else:
                newton = Newton(points=points,use_function=False,lower=arg.i[0], upper=arg.i[1])
    if arg.Lagrange:
        lagrange.run()
    else:
        newton.run()

            

if __name__ == "__main__":
    main()