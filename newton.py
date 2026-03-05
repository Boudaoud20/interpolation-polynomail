import numpy as np
import matplotlib.pyplot as plt
import random

class Newton:
    

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    def __init__(self,points:np.ndarray | dict ,use_function:bool,lower: float| int, upper:float| int):
        self.points=points
        self.use_function= use_function
        self.lower=lower
        self.upper=upper
        self.division_number:int = 0
        self.multiplication_number:int = 0
        self.addition:int = 0
        self.subtraction:int = 0

    def run(self):
        self.display_polynomial()  
        self.print_error_statistics()
        self.plot()
        print("addition count for newton interpolation & expanding newton polynomial: ",self.addition)
        print("subtraction count for newton interpolation & expanding newton polynomial: ",self.subtraction)
        print("multiplication count for newton interpolation & expanding newton polynomial: ",self.multiplication_number)
        print("division count for newton interpolation & expanding newton polynomial: ",self.division_number)
        print("total count of operations for newton interpolation & expanding newton polynomial: ",self.addition + self.subtraction + self.multiplication_number + self.division_number)


    def print_error_statistics(self):
        """
        Print error statistics for the current interpolation points.
        """
        if not self.use_function:
            print("Error statistics only available for function-based interpolation")
            return
            
        if isinstance(self.points, dict):
            pts = np.array(list(self.points.keys()))
        else:
            pts = self.points
            
        coeffs = self.divided_differences(pts)
        
        x_eval = np.arange(self.lower, self.upper, 0.05)
        y_true = self.f(x_eval)
        y_interp = np.array([self.eval_newton(xi, pts, coeffs) for xi in x_eval])
        
        errors = np.abs(y_interp - y_true)
        print(f"\nERROR STATISTICS for {len(pts)} interpolation points:")
        print(f"Error range: [{errors.min():.2e}, {errors.max():.2e}]")
        print(f"Root mean square error: {np.sqrt(np.mean(errors**2)):.2e}")
        print(f"Mean absolute error: {np.mean(errors):.2e}")


    def f(self,x: float | np.ndarray | int) -> float | np.ndarray:
        """
        Runge function.

        :param x: evaluation point(s)
        :type x: float | np.ndarray
        :return: f(x) = 1 / (1 + x^2)
        :rtype: float | np.ndarray
        """
        if self.use_function:
            #return 1 / (1 + (x**2))
            self.division_number += 1
            self.multiplication_number += 2
            self.addition += 1
            return 1 / (1 + (25*(x**2)))
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
                    
                    
                

    def divided_differences(self,points) -> np.ndarray:
        """
        Simpler divided differences computation using a 1D array.

        :return: Newton coefficients
        :rtype: np.ndarray
        """

        if self.use_function:
            values = self.f(points)
        else:
            if isinstance(self.points, dict):
                values = np.array([self.points[p] for p in points])
            else:
                values = self.points

        n = len(points)
        coeffs = values.astype(float).copy()

        for j in range(1, n):
            for i in range(n - 1, j - 1, -1):
                self.division_number += 1
                self.subtraction += 2
                coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (points[i] - points[i - j])

        return coeffs


    def eval_newton(self,x: float, points: np.ndarray, coeffs: np.ndarray) -> float:
        """
        Evaluate Newton interpolation polynomial at x, (substitutes the value x into the Newton interpolation polynomial and computes)

        :param x: evaluation point
        :type x: float

        :param coeffs: divided differences coefficients
        :type coeffs: np.ndarray

        :return: P(x)
        :rtype: float
        """

        px = coeffs[0]
        product = 1.0

        for i in range(1, len(points)):
            self.multiplication_number += 2
            self.subtraction += 1
            self.addition += 1
            product *= (x - points[i - 1])
            px += coeffs[i] * product
        return px


    def eval_newton_basis(self, x: float | np.ndarray, points: np.ndarray, k: int) -> float | np.ndarray:
        """
        Evaluate the k-th Newton basis function at x.
        N_k(x) = (x - x_0)(x - x_1)...(x - x_{k-1})

        :param x: evaluation point(s)
        :type x: float | np.ndarray
        :param points: interpolation nodes
        :type points: np.ndarray
        :param k: degree of basis function
        :type k: int
        :return: basis function value(s)
        :rtype: float | np.ndarray
        """
        if k == 0:
            if isinstance(x, np.ndarray):
                return np.ones_like(x)
            else:
                return 1.0
        
        result = 1.0
        for i in range(k):
            self.multiplication_number += 1
            self.subtraction += 1
            result *= (x - points[i])
        return result


    def expand_newton_to_polynomial(self, points: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
        """
        Expand Newton form to simple polynomial form.
        Returns coefficients in descending order of powers.

        :param points: interpolation nodes
        :type points: np.ndarray
        :param coeffs: divided differences coefficients
        :type coeffs: np.ndarray
        :return: polynomial coefficients [a_n, a_(n-1), ..., a_1, a_0]
        :rtype: np.ndarray
        """
        n = len(points)
        poly = np.array([coeffs[0]], dtype=float)
        
        for i in range(1, n):
            self.multiplication_number += 1
            self.addition += 3
            self.subtraction += 1
            new_poly = np.zeros(len(poly) + 1)
            new_poly[:-1] += poly
            new_poly[1:] += -points[i-1] * poly
            new_poly[-1] += coeffs[i]
            poly = new_poly

        return poly


    def eval_horner(self, coeffs: np.ndarray, x: float) -> float:
        result = coeffs[0]
        for i in range(1, len(coeffs)):
            result = result * x + coeffs[i]
        return result

    def get_horner_string(self, coeffs):
        if len(coeffs) == 1:
            return f"{coeffs[0]:.2e}"
        return f"{coeffs[0]:.2e} + x * ({self.get_horner_string(coeffs[1:])})"


    def display_polynomial(self) -> None:
        """
        Display the Newton polynomial in both complex and simple forms.

        :return: None
        :rtype: None
        """

        if isinstance(self.points, dict):
            points = np.array(list(self.points.keys()))
        else:
            points = self.points if isinstance(self.points, np.ndarray) else np.array(self.points)
        
        coeffs = self.divided_differences(points)
        
        print("NEWTON INTERPOLATION POLYNOMIAL")
        
        print("\nCOMPLEX FORM - Newton's Divided Difference Form:")
        print("P(x) = ", end="")
        
        terms = []
        for i in range(len(coeffs)):
            if i == 0:
                terms.append(f"a₀")
            else:
                product_str = " * ".join([f"(x - {points[j]:.2e})" for j in range(i)])
                terms.append(f"a{i} * {product_str}")
        
        print(" + ".join(terms))
        
        print("\nWhere divided differences coefficients are:")
        for i, c in enumerate(coeffs):
            if i == 0:
                print(f"  a₀ = {c:.10f}")
            else:
                product_str = " * ".join([f"(x - {points[j]:.2e})" for j in range(i)])
                print(f"  a{i} = {c:.10f}")
        
        print("\nExpanded Complex Form:")
        print("P(x) = ", end="")
        terms = []
        for i in range(len(coeffs)):
            if coeffs[i] != 0:
                if i == 0:
                    terms.append(f"{coeffs[i]:.2e}")
                else:
                    product_str = " * ".join([f"(x - {points[j]:.2e})" for j in range(i)])
                    terms.append(f"{coeffs[i]:.2e} * {product_str}")
        
        print(" + ".join(terms))

        print("\nSIMPLE FORM - Expanded Polynomial:")
        simple_coeffs = self.expand_newton_to_polynomial(points, coeffs)
        
        print("P(x) = ", end="")
        degree = len(simple_coeffs) - 1
        terms = []
        for i, c in enumerate(simple_coeffs):
            power = degree - i
            if abs(c) > 1e-10: 
                if power == 0:
                    terms.append(f"{c:.2e}")
                elif power == 1:
                    terms.append(f"{c:.2e}*x")
                else:
                    terms.append(f"{c:.2e}*x^{power}")
        
        print(" + ".join(terms).replace(" + -", " - "))

        print("HORNER'S FORM - Nested Multiplication:")
        print("P(x) = ", end="")
        horner_str = self.get_horner_string(simple_coeffs)
        print(horner_str)
        


    def chebyshev_nodes(self, n: int) -> np.ndarray:
        """
        Generate Chebyshev interpolation nodes on [lower, upper].

        :param n: polynomial degree
        :type n: int

        :return: Chebyshev nodes
        :rtype: np.ndarray
        """
        k: np.ndarray = np.arange(0, n + 1)
        return 0.5 * (self.lower + self.upper) + 0.5 * (self.upper - self.lower) * np.cos(((2 * k + 1) / (2 * (n + 1))) * np.pi)


    def plot_newton(self,points: np.ndarray, label: str, color: str, x: np.ndarray, ax) -> None:
        """
        Plot Newton interpolation polynomial.

        :param label: curve label
        :type label: str

        :param color: plot color
        :type color: str

        :param x: evaluation grid
        :type x: np.ndarray

        :param ax: matplotlib axis
        """
        if self.use_function:
            values = self.f(points)
        else:
            if isinstance(self.points, dict):
                values = np.array([self.points[p] for p in points])
            else:
                values = self.points
        coeffs = self.divided_differences(points=points)
        y_newton = np.array([self.eval_newton(xi, points, coeffs) for xi in x])

        ax.plot(x, y_newton, color=color, linewidth=1.5, label=label)
        ax.scatter(points, values, color=color, s=60, zorder=5)


    def plot(self) -> None:
        """
        Plot comparison between equidistant and Chebyshev interpolation
        using one column of subplots.  When a function is provided we also
        draw the true graph of f(x) and error plots; when only data points are supplied
        the plotting is limited to the interpolating polynomials.

        For function mode the rows correspond to a sequence of degrees, e.g.
        n=5,7,9,11; each row shows two subplots: left shows the original function
        together with the Newton polynomials built from equidistant and Chebyshev
        nodes; right shows the absolute interpolation errors for both node types.
        When operating on explicit points the same structure is used but the function
        graph and error plot are omitted.
        """
        degrees = [5, 7, 9, 11]

        if self.use_function:
            x = np.arange(self.lower, self.upper, 0.05)
            y = self.f(x)

            for deg in degrees:
                fig, (ax_interp, ax_error) = plt.subplots(1, 2, figsize=(16, 4.5))

                pts_equ = np.linspace(self.lower, self.upper, deg + 1)
                pts_cheb = self.chebyshev_nodes(deg)

                coeffs_equ = self.divided_differences(pts_equ)
                coeffs_cheb = self.divided_differences(pts_cheb)
                
                y_newton_equ = np.array([self.eval_newton(xi, pts_equ, coeffs_equ) for xi in x])
                y_newton_cheb = np.array([self.eval_newton(xi, pts_cheb, coeffs_cheb) for xi in x])
                
            
                error_equ = np.abs(y_newton_equ - y)
                error_cheb = np.abs(y_newton_cheb - y)

                ax_interp.plot(x, y, "k--", linewidth=2, label="f(x)=1/(1+25x^2)")
                ax_interp.plot(x, y_newton_equ, "r-", linewidth=1.5, label=f"Equidistant (n={deg})")
                ax_interp.plot(x, y_newton_cheb, "b-", linewidth=1.5, label=f"Chebyshev (n={deg})")
                ax_interp.scatter(pts_equ, self.f(pts_equ), color="red", s=40, zorder=5, marker="o")
                ax_interp.scatter(pts_cheb, self.f(pts_cheb), color="blue", s=40, zorder=5, marker="s")
                ax_interp.set_title(f"Degree {deg} - Interpolants")
                ax_interp.set_xlabel("x")
                ax_interp.set_ylabel("y")
                ax_interp.legend()
                ax_interp.grid(True, alpha=0.3)
                ax_interp.set_xticks(np.arange(self.lower, self.upper + 1, 1))

                ax_error.semilogy(x, error_equ, "r-", linewidth=1.5, label=f"Error Equidistant (n={deg})")
                ax_error.semilogy(x, error_cheb, "b-", linewidth=1.5, label=f"Error Chebyshev (n={deg})")
                ax_error.scatter(pts_equ, np.abs(self.f(pts_equ) - self.f(pts_equ)), color="red", s=40, zorder=5, marker="o")
                ax_error.scatter(pts_cheb, np.abs(self.f(pts_cheb) - self.f(pts_cheb)), color="blue", s=40, zorder=5, marker="s")
                ax_error.set_title(f"Degree {deg} - Absolute Error")
                ax_error.set_xlabel("x")
                ax_error.set_ylabel("Absolute Error")
                ax_error.legend()
                ax_error.grid(True, alpha=0.3, which="both")
                ax_error.set_xticks(np.arange(self.lower, self.upper + 1, 1))

        else:
            if isinstance(self.points, dict):
                all_pts = np.array(list(self.points.keys()))
            else:
                all_pts = np.array(self.points)

            x = np.arange(self.lower, self.upper, 0.05)

            for deg in degrees:
                if len(all_pts) >= deg + 1:
                    nodes = all_pts[: deg + 1]
                else:
                    nodes = all_pts.copy()

                coeffs = self.divided_differences(nodes)
                y_newton = np.array([self.eval_newton(xi, nodes, coeffs)
                                     for xi in x])

                if isinstance(self.points, dict):
                    values = np.array([self.points[p] for p in nodes])
                else:
                    values = nodes

                fig, ax = plt.subplots(figsize=(10, 4.5))
                ax.plot(x, y_newton, color="blue", linewidth=2,
                        label=f"Newton interpolant (n={len(nodes)-1})")
                ax.scatter(nodes, values, color="red", s=80, zorder=5,
                           marker="x",
                           label="data points")

                ax.set_title(f"Interpolation using first {len(nodes)} points")
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.legend()
                ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()