
import numpy
from numpy.linalg import lstsq


def _ramp(u):
    return numpy.maximum(u, 0)


def _step(u):
    return (u > 0).astype(float)


def fit_segmented_lineal_model(x_values, y_values, breakpoints, num_max_iterations=10):
    # Taken from
    # https://datascience.stackexchange.com/questions/8457/python-library-for-segmented-regression-a-k-a-piecewise-regression/32833#32833

    num_max_iterations = 10

    breakpoints = numpy.sort(numpy.array(breakpoints) )

    dt = numpy.min(numpy.diff(x_values))
    ones = numpy.ones_like(x_values)

    for i in range(num_max_iterations):
        # Linear regression:  solve A*p = Y
        Rk = [_ramp(x_values - xk) for xk in breakpoints]
        Sk = [_step(x_values - xk) for xk in breakpoints]
        A = numpy.array([ones, x_values] + Rk + Sk)
        p = lstsq(A.transpose(), y_values, rcond=None)[0] 

        # Parameters identification:
        a, b = p[0: 2]
        ck = p[2: 2 + len(breakpoints)]
        dk = p[2 + len(breakpoints):]

        # Estimation of the next break-points:
        new_breakpoints = breakpoints - dk / ck 

        # Stop condition
        if numpy.max(numpy.abs(new_breakpoints - breakpoints)) < dt / 5:
            break

        breakpoints = new_breakpoints
        print(breakpoints)
    else:
        print( 'maximum iteration reached' )

    # Compute the final segmented fit:
    x_solution = numpy.insert(numpy.append(breakpoints, max(x_values)), 0, min(x_values))
    ones = numpy.ones_like(x_solution) 
    Rk = [c * _ramp(x_solution - x0) for x0, c in zip(breakpoints, ck)]

    y_solution = a * ones + b * x_solution + numpy.sum(Rk, axis=0)

    return {'x_segments_edges': x_solution,
            'x_segments_edges': y_solution}


if __name__ == '__main__':

    x_values1 = numpy.linspace(0, 1)
    x_values2 = numpy.linspace(1, 2)
    y_values1 = x_values1 * 2
    y_values2 = x_values2 * 1 + 1
    x_values = numpy.array(list(x_values1) + list(x_values2))
    y_values = numpy.array(list(y_values1) + list(y_values2))
    breakpoints = numpy.array([0.1])

    fit_segmented_lineal_model(x_values, y_values, breakpoints)
