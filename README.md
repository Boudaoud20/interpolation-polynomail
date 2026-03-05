# Numerical Analysis Interpolation Tool

This program implements Lagrange and Newton polynomial interpolation methods for numerical analysis. It can interpolate data points or functions and generate plots of the results.

## Requirements

- Python 3.6 or higher
- NumPy
- Matplotlib

## Installation

1. Clone or download the project files.
2. Create a virtual environment (optional but recommended):
   ```
   python -m venv venv
   ```
3. Activate the virtual environment:
   - On Windows: `venv\Scripts\activate`
   - On macOS/Linux: `source venv/bin/activate`
4. Install the required packages:
   ```
   pip install numpy matplotlib
   ```

## Usage

Run the program from the command line using the following syntax:

```
python main.py [METHOD] [OPTIONS]
```

### Methods

- `-l, --Lagrange`: Use Lagrange interpolation method
- `-n, --Newton`: Use Newton interpolation method

### Options

- `-f`: Use a mathematical function instead of providing y-coordinates
- `-x`: X-coordinates (required, space-separated values)
- `-y`: Y-coordinates (required if not using -f, space-separated values)
- `-i`: Interval for plotting (required, two values: lower upper)

### Examples

#### Using Lagrange with data points
```
python main.py -l -x 1 2 3 4 -y 1 4 9 16 -i 0 5
```
This interpolates the points (1,1), (2,4), (3,9), (4,16) using Lagrange method and plots over the interval [0,5].

#### Using Newton with a function
```
python main.py -n -f -x 1 2 3 -i 0 5
```
This uses Newton method to interpolate a function at points x=1,2,3 and plots over [0,5].

## Output

The program will:
- Display the polynomial equation
- Generate a plot saved in the `images/` directory
- Print function values at x=2.9 and x=5.25

## Notes

- When using `-f`, do not provide `-y` coordinates
- The number of x and y coordinates must be equal when not using `-f`
- The interval must be specified as two numbers: lower bound followed by upper bound
- Plots are saved as image files in subdirectories under `images/`
