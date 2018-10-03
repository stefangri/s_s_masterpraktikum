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

    def app(self, name, value):
            if (type(value.magnitude) == uncertainties.core.Variable or type(value.magnitude) == uncertainties.core.AffineScalarFunc):
                s = '{:Lx}'.format(Q_(2, value.units)) + '~'
                value = value.magnitude
                df = DataFrame({'var': pd.Series(value, index = [name] ),
                'tex': f'{name} = \SI{{{best_value(value)}({sign_digits(value)})}}{{' + s[s.index('}{') + 2:s.index('~')]})

            else:
                s = '{:Lx}'.format(Q_(2, value.units)) + '~'
                value = value.magnitude
                df = DataFrame({'var': pd.Series(value, index = [name] ),
                'tex': f'{name} = \SI{{{value}}}{{' + s[s.index('}{') + 2:s.index('~')]})
                
            self.data = self.data.append(df, sort = True)
            with open(abs_path(f'results/result_{name}.tex'), 'w') as f:
                f.write('\\begin{equation} \n')
                f.write('\t' + self.data['tex'][name] + '\n')
                f.write('\label{eq: result_' +  name +  '}\n')
                f.write(r'\end{equation}')


    def makeresults(self):
        print(self.data['var'])
        with open(abs_path(f'results/{self.name}'), 'w') as f:
            for i in self.data['tex']:
                f.write(i + '\n')
