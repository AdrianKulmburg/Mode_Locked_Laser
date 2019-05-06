from numpy import *
from matplotlib.pyplot import *
from general_fit import *

def parabola(parameters, data):
    a = parameters[0]
    b = parameters[1]
    c = parameters[2]

    return a*data**2 + b*data + c

def line(parameters, data):
    a = parameters[0]
    b = parameters[1]

    return a*data + b

CW_Gal = flip(array([1750.,1600.,1300.,1040.,740.,615.,275.,183.,26.5,4.2,0.57]))
s_CW_Gal = flip(array([50.,20.,30.,30.,20.,30.,3,3,0.5,0.1,0.03]))
CW_Sil = flip(array([6800.,6010.,5000.,3450.,2640.,2000.,1000.,620.,90.,11.,1.04]))
s_CW_Sil = flip(array([300.,200.,200.,50.,30.,100.,100.,30.,3.,1.,0.01]))
ML_Gal = array([0.45,1.38,6.39,27.20,36.05,47.25,51.90,54.95,61,63,61.4])
s_ML_Gal = array([0.01,0.01,0.02,1,1,1,1,0.5,2,2,1]) + 0.0

attenuations = array([0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.0,2.0,3.0,4.0])
T = pow(10, -attenuations)

R = 10. # in kOhm

# ML-case
total_power = 513. # in mW or mV
s_total_power = 2.

errorbar(flip(T) * total_power, ML_Gal/R, xerr = flip(T)*s_total_power, yerr = s_ML_Gal/R, fmt = 'k.', label = 'Data')
title('Measured photodiode current for GaP diode, ML-regime', fontsize = 20)
xlabel(r'Power $P$ of the incoming laser beam in [mW]', fontsize = 18)
ylabel(r'Photodiode current $i_P$ in [$\mu$A]', fontsize = 18)

odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR = general_fit(flip(T) * total_power, ML_Gal/R, parabola, polyfit(flip(T) * total_power, ML_Gal/R,2))
x = linspace(min(flip(T) * total_power), max(flip(T) * total_power), 1000)
plot(x, parabola(odr_parameter_ideal, x), 'b--', label = 'Fitted parabola')
legend(loc = 'best', fontsize = 14)
print odr_p_value
print odr_parameter_ideal[0] , '+-', odr_parameter_error[0]
savefig('ML_Gal_attenuations.pdf')
show()

# CW-case
total_power = 554.2
s_total_power = 10.0

x = linspace(min(flip(T) * total_power), max(flip(T) * total_power), 1000)

errorbar(flip(T) * total_power, CW_Gal/R, xerr = flip(T)*s_total_power, yerr = s_CW_Gal/R, fmt = 'k.', label = 'Data')
title('Measured photodiode current for GaP diode, CW-regime', fontsize = 20)
xlabel(r'Power $P$ of the incoming laser beam in [mW]', fontsize = 18)
ylabel(r'Photodiode current $i_P$ in [$\mu$A]', fontsize = 18)
odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR = general_fit(flip(T) * total_power, CW_Gal/R, line, polyfit(flip(T) * total_power, CW_Gal/R,1))
print odr_p_value
print odr_parameter_ideal[0] , '+-', odr_parameter_error[0]
plot(x, line(odr_parameter_ideal, x), 'b--', label = 'Fitted line')
legend(loc = 'best', fontsize = 14)
savefig('CW_Gal_attenuations.pdf')
show()

T = T * pow(10, -3.0)
x = linspace(min(flip(T) * total_power), max(flip(T) * total_power), 1000)
errorbar(flip(T) * total_power, CW_Sil/R, xerr = flip(T)*s_total_power, yerr = s_CW_Sil/R, fmt = 'k.', label = 'Data')
title('Measured photodiode current for Si diode, CW-regime', fontsize = 20)
xlabel(r'Power $P$ of the incoming laser beam in [mW]', fontsize = 18)
ylabel(r'Photodiode current $i_P$ in [$\mu$A]', fontsize = 18)
odr_parameter_ideal, odr_parameter_error, odr_p_value, odr_SSR = general_fit(flip(T) * total_power, CW_Sil/R, line, polyfit(flip(T) * total_power, CW_Sil/R,1))
print odr_p_value
print odr_parameter_ideal[0] , '+-', odr_parameter_error[0]
plot(x, line(odr_parameter_ideal, x), 'b--', label = 'Fitted line')
legend(loc = 'best', fontsize = 14)
savefig('CW_Sil_attenuations.pdf')
show()
