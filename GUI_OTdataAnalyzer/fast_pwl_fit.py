# piecewise_linear_fit.py

# -- coding: utf-8 --
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cProfile
import pstats
import io
from numba import njit, prange
import argparse
import sys
import csv

class FastPWLFit:
    def __init__(self, x, y):
        """
        Initialize the FastPWLFit class with x and y data.

        Parameters:
        - x: array-like, independent variable.
        - y: array-like, dependent variable.
        """
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.n = len(x)
        self.breakpoints = None
        self.coefs = None

    @staticmethod
    @njit
    def compute_cumulative_sums(x, y):
        """
        Compute cumulative sums required for efficient error computation.

        Returns:
        - cumulative_x, cumulative_y, cumulative_xy, cumulative_xx, cumulative_yy
        """
        n = len(x)
        cumulative_x = np.zeros(n)
        cumulative_y = np.zeros(n)
        cumulative_xy = np.zeros(n)
        cumulative_xx = np.zeros(n)
        cumulative_yy = np.zeros(n)
        cumulative_x[0] = x[0]
        cumulative_y[0] = y[0]
        cumulative_xy[0] = x[0] * y[0]
        cumulative_xx[0] = x[0] * x[0]
        cumulative_yy[0] = y[0] * y[0]
        for i in range(1, n):
            cumulative_x[i] = cumulative_x[i-1] + x[i]
            cumulative_y[i] = cumulative_y[i-1] + y[i]
            cumulative_xy[i] = cumulative_xy[i-1] + x[i] * y[i]
            cumulative_xx[i] = cumulative_xx[i-1] + x[i] * x[i]
            cumulative_yy[i] = cumulative_yy[i-1] + y[i] * y[i]
        return cumulative_x, cumulative_y, cumulative_xy, cumulative_xx, cumulative_yy

    @staticmethod
    @njit(parallel=True)
    def fit_numba(n, n_segments, x, y, cumulative_x, cumulative_y, cumulative_xy, cumulative_xx, cumulative_yy):
        """
        Perform the dynamic programming algorithm to compute the optimal segmentation.

        Returns:
        - total_error: 2D array of cumulative errors.
        - segments: 2D array of segment breakpoints.
        """
        # Initialize total_error and segments arrays
        total_error = np.full((n_segments, n), 1e20)
        segments = np.full((n_segments, n), -1, dtype=np.int32)

        # Base case: first segment
        for j in prange(n):
            # Compute error for segment [0, j]
            sum_x = cumulative_x[j]
            sum_y = cumulative_y[j]
            sum_xy = cumulative_xy[j]
            sum_xx = cumulative_xx[j]
            sum_yy = cumulative_yy[j]
            n_points = j + 1

            denom = n_points * sum_xx - sum_x * sum_x
            if denom == 0.0:
                coef = 0.0
                intercept = sum_y / n_points
            else:
                coef = (n_points * sum_xy - sum_x * sum_y) / denom
                intercept = (sum_y - coef * sum_x) / n_points

            # Compute sum of squared errors
            error = sum_yy + coef * coef * sum_xx + n_points * intercept * intercept \
                    - 2 * coef * sum_xy - 2 * intercept * sum_y + 2 * coef * intercept * sum_x

            total_error[0, j] = error
            segments[0, j] = -1  # No previous segment

        # Dynamic programming for segments 2 to n_segments
        for s in range(1, n_segments):
            for j in prange(n):
                if j < s:  # Not enough points to have s+1 segments
                    continue
                min_error = 1e20
                best_k = -1
                for k in range(s-1, j-1):  # Ensure at least two points in the segment
                    current_error = total_error[s-1, k]  # Cost up to k with s segments

                    # Compute error for segment [k+1, j]
                    sum_x = cumulative_x[j] - (cumulative_x[k] if k >=0 else 0.0)
                    sum_y = cumulative_y[j] - (cumulative_y[k] if k >=0 else 0.0)
                    sum_xy = cumulative_xy[j] - (cumulative_xy[k] if k >=0 else 0.0)
                    sum_xx = cumulative_xx[j] - (cumulative_xx[k] if k >=0 else 0.0)
                    sum_yy = cumulative_yy[j] - (cumulative_yy[k] if k >=0 else 0.0)
                    n_points = j - k

                    denom = n_points * sum_xx - sum_x * sum_x
                    if denom == 0.0:
                        coef = 0.0
                        intercept = sum_y / n_points
                    else:
                        coef = (n_points * sum_xy - sum_x * sum_y) / denom
                        intercept = (sum_y - coef * sum_x) / n_points

                    # Compute sum of squared errors
                    error = sum_yy + coef * coef * sum_xx + n_points * intercept * intercept \
                            - 2 * coef * sum_xy - 2 * intercept * sum_y + 2 * coef * intercept * sum_x

                    total_cost = current_error + error

                    if total_cost < min_error:
                        min_error = total_cost
                        best_k = k

                # Assign the best cost and breakpoint
                total_error[s, j] = min_error
                segments[s, j] = best_k

        return total_error, segments

    def fit_model(self, n_segments, min_segment_length=2, max_segment_length=1000):
        """
        Fit the piecewise linear model with the specified number of segments.

        Parameters:
        - n_segments: int, number of segments to fit.
        - min_segment_length: int, minimum number of data points per segment.
        - max_segment_length: int, maximum number of data points per segment (optional).
        """
        if n_segments < 1:
            raise ValueError("Number of segments must be at least 1.")
        if n_segments > self.n:
            raise ValueError("Number of segments cannot exceed number of data points.")

        n = self.n

        # Compute cumulative sums using Numba
        cumulative_x, cumulative_y, cumulative_xy, cumulative_xx, cumulative_yy = self.compute_cumulative_sums(self.x, self.y)

        # Compute total_error and segments using Numba-compiled function
        total_error, segments = self.fit_numba(n, n_segments, self.x, self.y, cumulative_x, cumulative_y, cumulative_xy, cumulative_xx, cumulative_yy)

        # Backtrack to find breakpoints
        self.breakpoints = []
        self.coefs = []

        idx = n -1
        for s in range(n_segments-1, -1, -1):
            k = segments[s, idx]
            self.breakpoints.insert(0, idx)
            # Compute coef and intercept for this segment
            if k == -1:
                start = 0
            else:
                start = k +1
            end = idx +1  # Because Python slices are exclusive at the end
            # Ensure start < end and enforce minimum segment length
            if start >= end:
                # Empty segment, assign intercept as y[start-1] if possible
                intercept = self.y[start-1] if start-1 >=0 else 0.0
                coef = 0.0
            else:
                xi = self.x[start:end]
                yi = self.y[start:end]
                n_points = len(xi)
                if n_points < min_segment_length:
                    # Enforce minimum segment length
                    intercept = yi[0]
                    coef = 0.0
                else:
                    A = np.vstack([xi, np.ones(n_points)]).T
                    coef, intercept = np.linalg.lstsq(A, yi, rcond=None)[0]
            self.coefs.insert(0, (intercept, coef))
            idx = k

        self.breakpoints.insert(0, -1)  # Add -1 as the starting point

    def predict(self, x_vals):
        """
        Predict y values for the given x values based on the fitted model.

        Parameters:
        - x_vals: array-like, x values to predict y for.

        Returns:
        - y_pred: array-like, predicted y values.
        """
        if self.breakpoints is None or self.coefs is None:
            raise ValueError("Model has not been fitted yet. Call fit_model() first.")

        x_vals = np.asarray(x_vals)
        y_pred = np.zeros_like(x_vals)
        for i in range(len(self.breakpoints)-1):
            start = self.breakpoints[i] +1
            end = self.breakpoints[i+1] +1
            if start >= len(self.x):
                start_x = self.x[-1]
            else:
                start_x = self.x[start]
            if end > len(self.x):
                end_x = self.x[-1]
            else:
                end_x = self.x[end-1]
            intercept, coef = self.coefs[i]
            idx = (x_vals >= start_x) & (x_vals <= end_x)
            y_pred[idx] = intercept + coef * x_vals[idx]
        return y_pred

    def plot_fit(self):
        """
        Plot the original data and the fitted piecewise linear model.
        """
        if self.breakpoints is None or self.coefs is None:
            raise ValueError("Model has not been fitted yet. Call fit_model() first.")

        plt.figure(figsize=(8, 6))
        plt.scatter(self.x, self.y, s=5, color='black', label='Data')
        x_fit = self.x
        y_fit = self.predict(x_fit)
        plt.plot(x_fit, y_fit, color='yellow', label='Piecewise Linear Fit')
        plt.xlabel('Time(s)', fontsize=12)
        plt.ylabel('Basepairs(bp)', fontsize=12)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def save_results_csv(self, filename):
        """
        Save the fit results to a CSV file.

        The CSV will contain:
        - Segment Number
        - Start Index
        - End Index
        - Start X
        - End X
        - Slope
        - Intercept
        - R2
        - Duration
        - Fragment Length
        - SSR (Sum of Squared Residuals)
        - SST (Total Sum of Squares)

        Parameters:
        - filename: str, name of the output CSV file.
        """
        if self.breakpoints is None or self.coefs is None:
            raise ValueError("Model has not been fitted yet. Call fit_model() first.")

        results = []
        for i in range(len(self.coefs)):
            start = self.breakpoints[i] +1
            end = self.breakpoints[i+1] +1
            if end > len(self.x):
                end = len(self.x)
            xi = self.x[start-1:end]
            yi = self.y[start-1:end]
            if len(xi) ==0:
                slope = np.nan
                intercept = np.nan
                R2 = np.nan
                duration = np.nan
                fragment_length = np.nan
                ssr = np.nan
                sst = np.nan
            else:
                y_pred = self.predict(xi)
                ssr = np.sum((yi - y_pred) ** 2)
                sst = np.sum((yi - np.mean(yi)) ** 2)
                R2 = 1 - ssr / sst if sst !=0 else np.nan
                duration = xi[-1] - xi[0]
                fragment_length = yi[-1] - yi[0]
                slope = self.coefs[i][1]
                intercept = self.coefs[i][0]
            result = {
                'Segment': i+1,
                'Start Index': start-1,
                'End Index': end-1,
                'Start X': self.x[start-1] if start-1 < len(self.x) else np.nan,
                'End X': self.x[end-1] if end-1 < len(self.x) else np.nan,
                'Slope': slope,
                'Intercept': intercept,
                'R2': R2,
                'Duration': duration,
                'Fragment Length': fragment_length,
                'SSR': ssr,
                'SST': sst
            }
            results.append(result)

        # Convert to DataFrame and save
        df = pd.DataFrame(results)
        df.to_csv(filename, index=False)
        print(f"Fit results saved to {filename}")

def main():
    """
    Main function to execute the piecewise linear fit.

    Usage:
    python piecewise_linear_fit.py --input your_data.xlsx --format excel --segments 3 --output fit_results.csv
    """
    parser = argparse.ArgumentParser(description='Fast Piecewise Linear Fit')
    parser.add_argument('--input', type=str, default='your_data.xlsx', help='Input data file (Excel or CSV)')
    parser.add_argument('--output', type=str, default='fit_results.csv', help='Output CSV file for fit results')
    parser.add_argument('--segments', type=int, required=True, help='Number of segments to fit')
    parser.add_argument('--format', type=str, choices=['excel', 'csv'], default='excel', help='Format of input file')
    args = parser.parse_args()

    # Load your data
    input_filename = args.input
    try:
        if args.format == 'excel':
            data = pd.read_excel(input_filename)
        elif args.format == 'csv':
            data = pd.read_csv(input_filename)
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading '{input_filename}': {e}")
        sys.exit(1)

    # Check if required columns exist
    required_columns = ['time_pol', 'basepairs_pol']
    for col in required_columns:
        if col not in data.columns:
            print(f"Error: '{col}' column not found in the input file.")
            sys.exit(1)

    x = data['time_pol'].values
    y = data['basepairs_pol'].values

    # Number of segments to fit
    segment_number = args.segments

    if segment_number < 1:
        print("Error: Number of segments must be at least 1.")
        sys.exit(1)
    if segment_number > len(x):
        print("Error: Number of segments cannot exceed the number of data points.")
        sys.exit(1)

    # Initialize and fit the model
    pwlf = FastPWLFit(x, y)
    pwlf.fit_model(segment_number)

    # Predict and plot
    pwlf.plot_fit()

    # Save results to CSV
    pwlf.save_results_csv(args.output)

    # Optionally, get the slopes and intercepts
    slopes = [coef[1] for coef in pwlf.coefs]
    intercepts = [coef[0] for coef in pwlf.coefs]
    breakpoints = [pwlf.x[idx] for idx in pwlf.breakpoints[1:-1]]  # Exclude first and last indices

    # Print segment statistics
    print("\nSegment Statistics:")
    for i in range(len(slopes)):
        start = pwlf.breakpoints[i] + 1
        end = pwlf.breakpoints[i+1] +1
        if end > len(x):
            end = len(x)
        xi = x[start-1:end]
        yi = y[start-1:end]
        if len(xi) ==0:
            print(f"Segment {i+1}: Empty segment.")
            continue
        y_pred = pwlf.predict(xi)
        ssr = np.sum((yi - y_pred) ** 2)
        sst = np.sum((yi - np.mean(yi)) ** 2)
        R2 = 1 - ssr / sst if sst !=0 else np.nan
        duration = xi[-1] - xi[0]
        length = yi[-1] - yi[0]
        print(f"Segment {i+1}:")
        print(f"  Slope: {slopes[i]}")
        print(f"  Intercept: {intercepts[i]}")
        print(f"  R^2: {R2}")
        print(f"  Duration: {duration}")
        print(f"  Fragment Length: {length}")
        print(f"  SSR: {ssr}")
        print(f"  SST: {sst}")
        print()

def profile_main():
    """
    Profile the main function to identify bottlenecks.
    """
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(profiler, stream=s).sort_stats(sortby)
    ps.print_stats(10)  # Print top 10 functions
    print("\nProfiling Results:\n")
    print(s.getvalue())

if __name__ == "__main__":
    main()
    # To run with profiling, uncomment the following line:
    # profile_main()
