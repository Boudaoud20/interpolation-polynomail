import numpy as np
from numpy import polyval
import matplotlib.pyplot as plt
import random
class Lagrange:
    colors = [
    "red", "green", "blue", "cyan", "magenta", "yellow", "black", "white",
    "orange", "purple", "brown", "pink", "gray", "olive", "lime", "teal",
    "navy", "gold", "coral", "salmon", "turquoise", "indigo", "violet",
    "maroon", "chocolate", "crimson", "darkblue", "darkcyan", "darkgreen",
    "darkorange", "darkred", "deeppink", "dodgerblue", "firebrick",
    "forestgreen", "fuchsia", "gainsboro", "hotpink", "khaki",
    "lavender", "lightblue", "lightcoral", "lightgreen", "lightsalmon",
    "lightseagreen", "lightskyblue", "limegreen", "mediumblue",
    "mediumorchid", "mediumpurple", "mediumseagreen", "midnightblue",
    "orangered", "orchid", "palegreen", "peachpuff", "peru",
    "plum", "powderblue", "rosybrown", "royalblue", "sandybrown",
    "seagreen", "sienna", "skyblue", "slateblue", "springgreen",
    "steelblue", "tan", "thistle", "tomato", "wheat"]
    
    def __init__(self,points:tuple | dict ,use_function:bool,lower: float| int, upper:float| int):
        self.points = points
        self.use_function = use_function
        self.division_number:int = 0
        self.multiplication_number:int = 0
        self.addition:int = 0
        self.subtraction:int = 0
        self.lower= lower
        self.upper = upper
        self.lks:list=list()

    def run(self):

        coffs = self.pk()
        self.display_function(coffs)
        self.plot()
        print(f"f(2.9)= {self.pk_values(interval=[2.9])}")
        print(f"f(5.25)= {self.pk_values(interval=[5.25])}")

    def f(self,x: float| np.ndarray | int) -> float| np.ndarray:
        """
        this function substitute x value to mathematical function

        :param x: the value that will be substitute in the mathematical function
        :type x: float | ndarray | int
        :return: f(x)
        :rtype: float | ndarray
        :raises TypeError: throwed when we try to access to dict and we have tuple
        :raises ValueError: raised when x (or any element of x) does not exist in the points dict
        """
        if self.use_function:
            self.subtraction += 1
            return np.log(x) - ((2 * (x - 1)) / x)
        else:
            if isinstance(self.points,dict):
                if isinstance(x,(float,int)):
                    y = self.points.get(x)
                    if y is None:
                        raise ValueError("element not exist")
                    return y
                else:
                    subtituted_values = []
                    for item in x:
                        y=self.points.get(item)
                        if y is None:
                            raise ValueError("element not exist")
                        subtituted_values.append(y)
                    return np.array(subtituted_values)
            raise TypeError()
        
        

    def lk(self, k: int) -> tuple:
        """
        this funtion will calculate Lk(x) form of lagrange interpolation

        :param points: container of polynomial nodes
        :type points: float

        :param k: refer to iteration number k
        :type k: int

        :return: this function will return denominators with roots of the numerator 
        :rtype: tuple

        :raises TypeError: raised when we can not extract x values

        """
        self.division_number += 1
        if self.use_function and isinstance(self.points,tuple):
            points = tuple(self.points)
        elif isinstance(self.points,dict):
            points = tuple(key for key in self.points)
        else:
            raise TypeError()
        roots = np.array([])
        denominator = 1.0
        for j, point in enumerate(points):
            if j != k:
                self.multiplication_number += 1
                self.subtraction+=1
                denominator *= (points[k] - point)
                roots = np.append(roots, point)
        return denominator, roots

    def pk(self) -> np.ndarray:
        """
        this function will calculate the polynomail formula

        :param points: polynomail node
        :type points: float

        :return: will return array of coefficients
        :rtype: np.ndarray

        :raises TypeError: raised when we can not extract x values
        """
        self.lks = [] 
        if self.use_function and isinstance(self.points,tuple):
            points = tuple(self.points)
        elif isinstance(self.points,dict):
            points = [key for key in self.points]
        else:
            raise TypeError()
        coffs_list =[]
        for k, point in enumerate(points):
            denominator, roots = self.lk(k=k)
            self.lks.append(roots)
            coff = np.poly(roots)
            coffs_list.append(coff * (self.f(point) / denominator))
            self.division_number += 1
            self.multiplication_number+=1
        self.addition += len(coffs_list) 
        return np.sum(coffs_list, axis=0)

    def pk_values(self, interval:list) -> np.ndarray:
        """
        this function will substitute list of values in the calculated polynomial

        :param points: this refere to polynomial coefficients
        :type points: float

        :param interval: list of x values
        :type interval: list

        :return: Results of the substituted values
        :rtype: np.ndarray
        """
        return polyval(self.pk(), interval)

    def display_function(self,coeffs: np.ndarray) -> None:
        """
        this function will display the polynomial formula

        :param coeffs: represent coefficients values of the polynomial
        :type coeffs:  ndarray

        :return: none
        :rtype: None
        
        """
        if self.use_function:
            points = self.points
        else:
            points = tuple([item for item in self.points])
        for i,term in enumerate(self.lks):
            print(f"L{i}=[",end=" ")
            for sub_term in term:
                print(f"(x-{sub_term:.2e}) *  ",end=" ")
            print("/",end="")
            for j in range(len(self.lks[0])):
                if i!=j:
                    print(f"({points[i]-points[j]})",end=" ")
                    ...
            print("]\n",end=" ")

        print("\n\nsimple P(x) = ", end="")

        for i, point in enumerate(reversed(coeffs)):
            if point != 0:
                print(f"({point:.2e})*x^{i} + ",end=" ")
        print(f"\nDivision count for Lagrange interpolation: {self.division_number}")
        print(f"Subtraction count for Lagrange interpolation: {self.subtraction}")
        print(f"Multiplication count for Lagrange interpolation: {self.multiplication_number}")
        print(f"Addition count for Lagrange interpolation: {self.addition}")
        print(f"Total count of operations for Lagrange interpolation: {self.division_number + self.subtraction + self.multiplication_number + self.addition}")

    def lk_substituted(self,x:float | np.ndarray,k:int)->float | np.ndarray:
        """
        this function take scalar or vector and sunstitute in the in Li(k) forms

        :param x: refers the x coordinate value
        :type x: float | np.ndarray

        :param k: refers to the order of Li(x)
        :type k: int

        :return: this function will return y coordinate as scalar or vector
        :rtype: float | np.ndarray:

        """
        if self.use_function:
            points = self.points
        else:
            points = tuple([item for item in self.points])
        denominator:np.ndarray | float= 1.0
        numerator:np.ndarray | float =1.0
            
        for j in range(len(self.lks[k])):
            if k!=j:
                self.multiplication_number += 2
                self.subtraction+=1
                denominator*=points[k]-points[j]
                numerator*=(x-points[j])
        return numerator/denominator
                


    def plot(self) -> None:
        """
        this function will plot the real f(x), p(x) and polynomail node and the error + y=0

        :return: none
        :rtype: None
        
        """
        x = np.arange(self.lower,self. upper, 0.1)

        coffs = self.pk()
        #coffs[np.abs(coffs) < 1e-30] = 0
        y_poly = polyval(coffs, x)
        if self.use_function:
            y = self.f(x)
            error = np.abs(y - y_poly)
        self.display_error_range(error)
              
        if self.use_function and isinstance(self.points,tuple):
            nodes_x = np.array(list(self.points))
        elif isinstance(self.points,dict):
            nodes_x = np.array([key for key in self.points])
        nodes_y = self.f(nodes_x) 

        plt.figure(figsize=(10, 6))

        if self.use_function:
            plt.plot(x, y, c="black", linestyle="--", linewidth=2, label="f(x) = ln(x) - 2(x-1)/x")
            plt.xticks(np.arange(self.lower, self.upper + 1, 0.5))


        

        plt.plot(x, y_poly, c="blue", linewidth=1.5, linestyle="-", label="P(x) interpolation")

        plt.scatter(nodes_x, nodes_y, c="red", marker="x", s=80, zorder=5, label="Interpolation nodes")

        if self.use_function:
            plt.plot(x, error, c="orange", linestyle="--", linewidth=1, label="|f(x) - P(x)| error")
        plt.axhline(0, c="green", linestyle="-", linewidth=1, label="y=0")

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Lagrange Interpolation of f(x) = ln(x) - 2(x-1)/x")
        plt.xticks(np.arange(self.lower, self.upper + 1, 0.5))
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.figure(figsize=(10,6))
        color = random.sample(Lagrange.colors,k=len(self.lks))
        for i in range(len(self.lks)):
            plt.plot(x,self.lk_substituted(x,i),linestyle="-",linewidth=1,label=f"L{i}(x)",c=color[i])
        plt.title("Li(x)")
        plt.xlabel("x")
        plt.ylabel("y")

            
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    def display_error_range(self,interval:np.ndarray)->None:
        """
        this function will display the error range of the interpolation

        :param interval: list of y values
        :type interval: np.ndarray

        :return: none
        :rtype: None
        
        """
        min = np.min(interval)
        max = np.max(interval)
        print(f"Error range: [{min:.2e}, {max:.2e}]")
        