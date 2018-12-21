import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from math import floor, ceil


def create_line(point_1, point_2):
    m = (point_2[1] - point_1[1]) / (point_2[0] - point_1[0])
    b = point_2[1] - point_2[0] * m

    return sp.S('m*x + b').subs([(sp.S('m'), m), (sp.S('b'), b)])


def interpolate(x1, y1, x2, y2, x3):
    return (y2 - y1) * (x3 - x2) / (x2 - x1) + y2


def find_equilibrium_y(x):
    try:
        y = equilibrium_data['y'][equilibrium_data.loc[equilibrium_data['x'] == x].index.values[0]]

    except IndexError:
        index_1 = floor(x / 0.001001)
        index_2 = ceil(x / 0.001001)
        y = interpolate(equilibrium_data['x'][index_1], equilibrium_data['y'][index_1],
                        equilibrium_data['x'][index_2], equilibrium_data['y'][index_2], x)

    return y


def find_equilibrium_x(y):
    try:
        x = equilibrium_data['x'][equilibrium_data.loc[equilibrium_data['y'] == y].index.values[0]]

    except IndexError:
        equilibrium_points = equilibrium_data.ix[(equilibrium_data['y'] - y).abs().argsort()[:2]]
        equilibrium_points.index = [1, 0]
        x = interpolate(equilibrium_points['y'][0], equilibrium_points['x'][0],
                        equilibrium_points['y'][1], equilibrium_points['x'][1], y)

    return x


# Could determine using automated mass balance system
compositions = {'Distillate': 0.854, 'Bottoms': 0.0094, 'Feed': 0.044}
flows = {'Distillate': 519.6,
         'Bottoms': 12214.5,
         'Feed': 12734.1}

# Could automate calc
q = 1.11

equilibrium_data = pd.read_excel('CHE411_Project.xlsx', 'xy')

q_slope = q / (q - 1)
q_intercept = -compositions['Feed'] / (q - 1)
q_line = create_line((compositions['Feed'], compositions['Feed']),
                     (0.8, q_slope * 0.8 + q_intercept))
count = 0

for count in range(len(equilibrium_data)):
    if equilibrium_data['y'][count] - q_line.subs(sp.S('x'), equilibrium_data['x'][count]) < 0:
        break

linear_interpolation = create_line((equilibrium_data['x'][count - 1], equilibrium_data['y'][count - 1]),
                                   (equilibrium_data['x'][count], equilibrium_data['y'][count]))
q_line_eqlbm = sp.solve(linear_interpolation - q_line)[0]

R_min_line = create_line((q_line_eqlbm, q_line.subs(sp.S('x'), q_line_eqlbm)),
                         (compositions['Distillate'], compositions['Distillate']))
R_min_intercept = R_min_line.subs(sp.S('x'), 0)

x_check = compositions['Distillate']


while x_check > 0.4:
    equilibrium_value = find_equilibrium_y(x_check)

    R_min_value = R_min_line.subs(sp.S('x'), x_check)

    if R_min_value <= equilibrium_value:
        x_check -= 0.05

    else:
        R_min_intercept -= 0.01
        R_min_line = create_line((0, R_min_intercept), (compositions['Distillate'], compositions['Distillate']))
        x_check = compositions['Distillate']

    if R_min_intercept < 0:
        print('Could not ascertain R min')
        break

R_min = compositions['Distillate'] / R_min_intercept - 1
R = 2.8 * R_min

flows['Enriching Liquid'] = R * flows['Distillate']
flows['Enriching Gas'] = flows['Enriching Liquid'] + flows['Distillate']
flows['Stripping Liquid'] = q * flows['Feed'] + flows['Enriching Liquid']
flows['Stripping Gas'] = (q - 1) * flows['Feed'] + flows['Enriching Gas']

strip_intercept = -flows['Bottoms'] * compositions['Bottoms'] / flows['Stripping Gas']
strip_slope = flows['Stripping Liquid'] / flows['Stripping Gas']
stripping_line = create_line((0, strip_intercept),
                             (0.8, strip_slope * 0.8 + strip_intercept))

x_check = compositions['Distillate']
previous_x_check = x_check
x_checks = [x_check]
y_check = interpolate(0, compositions['Distillate'] / (R + 1),
                      compositions['Distillate'], compositions['Distillate'], x_check)
y_checks = [y_check, y_check]
stages = 0
intersect = sp.solve(stripping_line - q_line)[0]

while x_check > compositions['Bottoms']:

    print(x_check)

    if x_check == previous_x_check:
        x_check = find_equilibrium_x(y_check)
        x_checks.append(x_check)
        stages += 1

    else:
        previous_x_check = x_check

        if x_check > intersect:
            y_check = interpolate(0, compositions['Distillate'] / (R + 1),
                                  compositions['Distillate'], compositions['Distillate'], x_check)

        else:
            y_check = interpolate(0, strip_intercept,
                                  0.8, strip_slope * 0.8 + strip_intercept, x_check)

        y_checks.append(y_check)

equilibrium_data['OpLine'] = pd.Series()
equilibrium_data['q'] = pd.Series()
equilibrium_data['r_min'] = pd.Series()

for data_point in equilibrium_data.index:
    if 0 <= equilibrium_data['x'][data_point] <= intersect:
        equilibrium_data['OpLine'][data_point] = interpolate(0, strip_intercept,
                                                             0.8, strip_slope * 0.8 + strip_intercept,
                                                             equilibrium_data['x'][data_point])

    else:
        equilibrium_data['OpLine'][data_point] = interpolate(0, compositions['Distillate'] / (R + 1),
                                                             compositions['Distillate'], compositions['Distillate'],
                                                             equilibrium_data['x'][data_point])

    if compositions['Feed'] <= equilibrium_data['x'][data_point] <= intersect:
        equilibrium_data['q'][data_point] = q_line.subs(sp.S('x'), equilibrium_data['x'][data_point])

    else:
        equilibrium_data['q'][data_point] = None

    equilibrium_data['r_min'][data_point] = R_min_line.subs(sp.S('x'), equilibrium_data['x'][data_point])

plt.plot(equilibrium_data['x'], equilibrium_data['y'], 'k')
plt.plot(equilibrium_data['x'], equilibrium_data['OpLine'], 'b')
plt.plot(equilibrium_data['x'], equilibrium_data['q'], 'g')
plt.plot(equilibrium_data['x'], equilibrium_data['x'], 'k')
plt.plot(equilibrium_data['x'], equilibrium_data['r_min'], 'b--', linewidth=0.5)

plt.step(x_checks, y_checks, 'r')

plt.xlim([0, 1])
plt.ylim([0, 1])
plt.xlabel('Water Mol Fraction')
plt.ylabel('Ethanol Mol Fraction')
plt.title('McCabe-Thiele Diagram for Ethanol-Water Distillation Column')

plt.show()

print('')
print(R / R_min, stages)
