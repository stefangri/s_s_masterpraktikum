from pandas import Series, DataFrame
import pandas as pd
import collections
import numpy
import uncertainties
import pint
from uncertainties import ufloat
from pint import UnitRegistry
import os.path
ureg = UnitRegistry()
Q_ = ureg.Quantity



def return_int(num):
    num_str = str(num)
    
    num_str = num_str.split('.')[1]
    num_str = num_str[0:1]
    return int(num_str)

def abs_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)


def sign_digits(num):
    chars = f'{num}'
    chars = chars.split('+/-')[-1]
    chars = list(chars)
    while '.' in chars: chars.remove('.')
    if 'e' in chars:
        return chars[0] + chars[1]
    else:
        return chars[-2] + chars[-1]

def best_value(num):
    value = f'{num}'.split('+/-')[0]
    return value




class Latexdocument(object):
    def __init__(self, filename):
        self.name = filename
        self.data = DataFrame(columns=(['tex', 'var']))
    def tabular(self, data, header, places, caption, label):
        with open(self.name, 'w') as f:
            f.write('\\begin{table} \n\\centering \n\\caption{' + caption + '} \n\\label{tab: ' + label + '} \n\\begin{tabular}{')
            for i in range(0, len(data)):
                if type(data[i][0]) == uncertainties.core.Variable:
                    f.write('S[table-format=' + str(places[i][0]) + ']@{${}\pm{}$} S[table-format=' + str(places[i][1]) + '] ')
                else:
                    f.write('S ')

            f.write('} \n\\toprule  \n')

            for i in range(0, len(data)):
                if i == len(data) - 1:
                    if type(data[i][0]) == uncertainties.core.Variable:
                        f.write('\multicolumn{2}{c}{$' + header[i][0:header[i].find('/')] +  '\:/\: \si{' + header[i][header[i].find('/')+1:] + '}$} \\\ \n')
                    else:
                        f.write('{$' + header[i][0:header[i].find('/')] +  '/ \si{' + header[i][header[i].find('/')+1:] + '}$} \\\ \n')
                else:
                    if type(data[i][0]) == uncertainties.core.Variable:
                        f.write('\multicolumn{2}{c}{$' + header[i][0:header[i].find('/')] +  '\:/\: \si{' + header[i][header[i].find('/')+1:] + '}$} & ')
                    else:
                        f.write('{$' + header[i][0:header[i].find('/')] +  '/ \si{' + header[i][header[i].find('/')+1:] + '}$} & ')


            f.write('\\midrule  \n')
            for i in range(0, len(data[0])):
                for j in range(0, len(data)):
                    if type(data[j][0]) == uncertainties.core.Variable:
                        if j == len(data) - 1:
                            f.write(('{:.' + str(return_int(places[j][0])) + 'f} ' + '& {:.' + str(return_int(places[j][1])) + 'f}' + '\\\ \n').format(data[j][i].n, data[j][i].s))
                        else:
                            f.write(('{:.' + str(return_int(places[j])) + 'f} ' + '& {:.' + str(return_int(places[j][1])) + 'f}'+ ' & ').format(data[j][i].n, data[j][i].s))
                    else:
                        if j == len(data) - 1:
                            f.write(('{:.' + str(places[j]) + 'f}' + '\\\ \n').format(data[j][i]))
                        else:
                            f.write(('{:.' + str(places[j]) + 'f}' + ' & ').format(data[j][i]))
            f.write('\\bottomrule \n\\end{tabular} \n\\end{table}')

    def add_result(self, name, value):
            if (type(value.magnitude) == uncertainties.core.Variable or type(value.magnitude) == uncertainties.core.AffineScalarFunc):
                latex_value = f'{value.magnitude:+.2uS}'
                uncert = latex_value.split('(')[1]
                uncert = uncert.split(')')[0]
                uncert = uncert.replace('.', '')
                latex_value = latex_value.split('(')[0] + f'({uncert})' + latex_value.split(')')[-1]
                latex_unit = (f'{value:Lx}'.split('}{'))[1].split('}')[0]
                value = value.magnitude
                df = DataFrame({'var': pd.Series(value, index = [name]),
                'tex': '\SI{' + latex_value + '}{' + latex_unit + '}'})

            else:
                latex_unit = (f'{value:Lx}'.split('}{'))[1].split('}')[0]
                latex_value = str(value.magnitude)
                df = DataFrame({'var': pd.Series(value, index = [name] ),
                'tex': '\SI{' + latex_value + '}{' + latex_unit + '}'})

            self.data = self.data.append(df, sort = True)
            with open(abs_path('results/result_' + name.replace('\\', '') + '.tex'), 'w') as f:
                #f.write('\\begin{equation} \n')
                f.write(self.data['tex'][name] + '\n')
                #f.write('\label{eq: result_' +  name +  '}\n')
                #f.write(r'\end{equation}')


    
