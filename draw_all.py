from numpy import *
from matplotlib.pyplot import *
from general_fit import general_fit

def data_to_function_raw(data_x, data_y, z_input):
    if z_input <= data_x[0] or z_input >= data_x[-1]:
        return 0.0
    left_index = where(data_x <= z_input)[0][-1]
    right_index = left_index + 1

    x_1 = data_x[left_index]
    x_2 = data_x[right_index]
    y_1 = data_y[left_index]
    y_2 = data_y[right_index]

    return y_1 + (y_2 - y_1 + 0.0) / (x_2 - x_1 + 0.0) * (z_input - x_1)

def width_at_half_height(x, y, peak_x, peak_y):
    left_indices = where(x < peak_x)[0]
    right_indices = where(x > peak_x)[0]

    x_left = x[left_indices]
    y_left = y[left_indices]
    left = x_left[where(y_left <= peak_y/2.0)[0][-1]]

    x_right = x[right_indices]
    y_right = y[right_indices]
    right = x_right[where(y_right <= peak_y/2.0)[0][0]]
    return left, right

def draw_fitting(model, measured, offset, graph_title, location):
    x = model[0]
    y = (model[1] - min(model[1])) / (max(model[1]) - min(model[1]) + 0.0)

    data_x = measured[0]
    data_y = (measured[1] - offset) / (max(measured[1]) - offset + 0.0)

    function = vectorize(lambda z : data_to_function_raw(x, y, z))

    function_total = lambda parameters, data: function(parameters[0] + parameters[1] * data)
    # parameters[0] = multiplying factor
    # parameters[1] = offset
    parameter_ideal, parameter_error, p_value, SSR = general_fit(data_x, data_y, function_total, [0.0, 1.0])

    peak_x = data_x[argmax(data_y)]
    peak_y = data_y[argmax(data_y)]
    left, right = width_at_half_height(data_x, data_y, peak_x, peak_y)
    whm = (right - left)/parameter_ideal[1]
    s_whm = (right - left) * parameter_error[1] / parameter_ideal[1]**2

    
    plot(x, y, 'b,-', label = 'Reference')
    plot((data_x - parameter_ideal[0] + 0.0)/parameter_ideal[1], data_y, 'k.', label = 'Measured signal')
    legend(loc = 2, fontsize = 17)

    peak = (data_x[argmax(data_y)] - parameter_ideal[0] + 0.0) / parameter_ideal[1]
    s_peak = sqrt(parameter_error[0]**2 / (parameter_ideal[1]**2 + 0.0) + (data_x[argmax(data_y)] -parameter_ideal[0])**2 * parameter_error[1]**2 / (parameter_ideal[1]**4 +0.0))
    result = '    Peak at\n' + ('%4.1f' % peak) + r'$\pm$' + ('%3.1f' % s_peak) + r' nm'
    if p_value > 0.99:
        result = result + '\n p-value > 0.99'
    else:
        result = result + ('\np-value: %3.2f' % p_value)
    text(peak+20, 0.8, result, fontsize = 18)
    text(peak-140, 0.5, '\nWHM = ' + ('%3.1f' % whm) + r'$\pm$' + ('%3.1f' % s_whm) + r' nm', fontsize = 18)
    xlim((650, 900))
    ylim((0.0, 1.0))
    title(graph_title, fontsize = 24)
    ylabel('Relative Power', fontsize = 24)
    xlabel('Wavelength in [nm]', fontsize = 24)
    savefig(location)
    show()


# So, there are two plots we have to draw in total.
# One for the mode-locked one and one for the continuous one.

# Loading files
model_CW_file = open('161202_cw')
model_ML_file = open('161202_modelocked')

measured_CW_file = open('CW_Messung_90')
measured_ML_file = open('ML_Messung_90')

# Extracting data
model_CW = loadtxt(model_CW_file)
model_ML = loadtxt(model_ML_file)

measured_CW = loadtxt(measured_CW_file)
measured_ML = loadtxt(measured_ML_file)

offset = median(measured_CW[1])

draw_fitting(model_CW, measured_CW, offset, 'Fitted spectral distribution, CW-Mode', 'CW.pdf')
draw_fitting(model_ML, measured_ML, offset, 'Fitted spectral distribution, ML-Mode', 'ML.pdf')
