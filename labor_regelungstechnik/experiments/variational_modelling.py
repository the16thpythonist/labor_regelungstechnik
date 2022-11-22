import os

import sympy
import sympy as sp
import sympy.physics.mechanics as spd

import numpy as np
import matplotlib.pyplot as plt
from pycomex.experiment import Experiment
from pycomex.util import Skippable

from labor_regelungstechnik.utils import TEMPLATE_ENV
from labor_regelungstechnik.utils import render_latex, latex_math


BASE_PATH = os.getcwd()
NAMESPACE = 'variational_modeling'
DEBUG = True
with Skippable(), (e := Experiment(base_path=BASE_PATH, namespace=NAMESPACE, glob=globals())):

    # https://stackoverflow.com/questions/25346132
    spd.init_vprinting()

    def save_expression(expression, name: str, template_name='math.tex.j2'):
        latex = spd.mlatex(expression)
        latex = latex_math(latex, template_name=template_name)
        pdf_path = os.path.join(e.path, name)
        render_latex({'content': latex}, pdf_path)

    def save_equations(equations: list, name: str, template_name='equations.tex.j2'):
        latex_list = [spd.mlatex(eq) for eq in equations]
        latex = latex_math(latex_list, template_name=template_name)
        pdf_path = os.path.join(e.path, name)
        render_latex({'content': latex}, pdf_path)

    # - DEFINING THE ENERGY FUNCTION
    e.info('starting out with the variational modeling...')

    m_x = sp.Symbol('m_x')
    m_y = sp.Symbol('m_y')
    y_max = sp.Symbol('y_max')
    g = sp.Symbol('g')
    c_varphi = sp.Symbol('c_varphi')
    c_x = sp.Symbol('c_x')
    k_x = sp.Symbol('k_x')
    k_l = sp.Symbol('k_l')
    v_x = sp.Symbol('v_x')
    v_l = sp.Symbol('v_l')
    l_0 = sp.Symbol('l_0')

    k_vx = sp.Symbol('k_vx')
    k_vl = sp.Symbol('k_vl')
    k_phi = sp.Symbol('k_phi')

    params_default_map = {
        'm_x': 30,
        'm_y': 3,
        'y_max': 1.2,
        'g': 9.81,
        'c_varphi': 0.5,
        'c_x': 0.5,
        'k_x': 200,
        'k_l': 200,
        'l_0': 0.1,
        'k_vx': 4,
        'k_vl': 2,
        'k_phi': 0.4
    }

    t = sp.Symbol('t')
    x, l, phi = spd.dynamicsymbols(r'x l \varphi')
    X, L, Phi = spd.dynamicsymbols(r'X L \phi')
    #l = sp.Symbol('l')

    l_sum = (l_0 + l)

    general_coordinates_map = {
        str(x): x,
        str(l): l,
        str(phi): phi
    }
    coordinate_derivatives_map = {
        str(x): X,
        str(l): L,
        str(phi): Phi
    }
    coordinate_rhs_map = {
        str(x): k_x * (v_x - x.diff(t)),
        str(l): k_l * (v_l - l_sum.diff(t)),
        str(phi): - m_y * g * sp.sin(phi) * l_sum - c_varphi * phi.diff(t)
    }

    energy = (sp.Rational(1, 2) * (m_x + m_y) * x.diff(t)**2) \
    + (sp.Rational(1, 2) * m_y * (l_sum.diff(t) * sp.sin(phi) - x.diff(t) + l_sum * phi.diff(t) * sp.cos(phi))**2) \
    + (sp.Rational(1, 2) * m_y * (-l_sum.diff(t) * sp.cos(phi) + l_sum * phi.diff(t) * sp.sin(phi))**2)

    # energy = sp.simplify(energy)

    save_expression(energy, '_energy.pdf')

    # - IMPLEMENT LAGRANGE EQUATION
    equation_map = {}
    for variable in general_coordinates_map.values():
        file_name = str(variable).replace('\\', '').replace('(t)', '')
        energy_derivative = sp.Derivative(energy, variable)
        save_expression(energy_derivative, f'_energy_derivative_{file_name}.pdf')

        variable_dot = variable.diff(t)
        energy_dot_derivative = sp.Derivative(energy, variable_dot)
        save_expression(energy_dot_derivative, f'_energy_dot_derivative_{file_name}.pdf')

        lagrange_term = sp.simplify(sp.Derivative(energy_dot_derivative, t) - energy_derivative)
        #lagrange_term = sp.Derivative(energy_dot_derivative, t) - energy_derivative
        rhs_term = coordinate_rhs_map[str(variable)]
        lagrange_equation = sp.Eq(lagrange_term, rhs_term)
        equation_map[str(variable)] = lagrange_equation

        save_expression(lagrange_equation, f'_lagrange_equation_{file_name}.pdf')

    e.info(f'produced {len(equation_map)} equations for each of the variables: {equation_map.keys()}')

    # - CONVERTING TO ODE SYSTEM

    # Doing all the substitutions so that there are no derivatives of the second order anymore
    substituted_equation_map = {}
    final_equations = []
    for name, equation in equation_map.items():
        substituted_equation = equation.copy()

        for variable_name in general_coordinates_map.keys():
            variable = general_coordinates_map[variable_name]
            variable_d = variable.diff(t)  # first derivative w.r.t time
            variable_dd = variable.diff(t).diff(t)  # second derivative w.r.t. time

            # ~ Substituting all higher derivatives of the base coordinates with explicit variables
            substitute = coordinate_derivatives_map[variable_name]
            substitute_d = substitute.diff(t)

            substituted_equation = substituted_equation.subs(variable_d, substitute)
            # substituted_equation = substituted_equation.subs(variable_dd, substitute_d)

        substituted_equation = sp.simplify(substituted_equation)
        substituted_equation_map[name] = substituted_equation

        final_equations.append(substituted_equation)

        variable = general_coordinates_map[name]
        variable_d = variable.diff(t)
        substitute = coordinate_derivatives_map[name]
        substitution_equation = sp.Eq(variable_d, substitute)
        final_equations.append(substitution_equation)

        file_name = name.replace('\\', '').replace('(t)', '')
        save_expression(substituted_equation, f'substituted_equation_{file_name}.pdf')

    save_equations(final_equations, 'final_equations.pdf')

    # '# Now we have to solve these equations for the derivative term
    # final_equations = []
    # for name, equation in substituted_equation_map.items():
    #     variable = general_coordinates_map[name]
    #     variable_d = variable.diff(t)
    #
    #     substitute = coordinate_derivatives_map[name]
    #     substitute_d = substitute.diff(t)
    #
    #     # Obviously we can only solve for the term if the term is actually contained in the expression,
    #     # which depending on the system structure is not always the case for all of the variables...
    #     contains_derivative = substitute_d in equation.lhs.atoms(sp.Expr)
    #     if contains_derivative:
    #         solved = sp.solve(equation, substitute_d)
    #         solved = sp.simplify(solved[0])
    #         solved_equation = sp.Eq(substitute_d, solved)
    #         substitution_equation = sp.Eq(variable_d, substitute)
    #         final_equations.append(substitution_equation)
    #         final_equations.append(solved_equation)
    #
    #         file_name = name.replace('\\', '').replace('(t)', '')
    #         save_expression(solved, f'solved_{file_name}.pdf')
    #
    # save_equations(final_equations, f'system.pdf')'

    # - CLEANING UP FOR CODE GENERATION
    # So now that we essentially have created our system of ordinary differential equations we would like
    # to simulate this numerically. For that to happen we would need this system of equations not in
    # symbolic form but instead as python code that we can pipe into numpy.
    # Sympy supports automatic code generation which we can use to achieve just that, but to make this
    # work out, we first have to get rid of all the symbolic "derivative" terms by replacing them with
    # just plain old symbols.

    # First of all we create symbols which replace the time functions
    variable_symbol_map = {}
    plain_symbols = []
    derivative_symbols = []
    index = 0
    for variable in [*coordinate_derivatives_map.values(), *general_coordinates_map.values()]:
        variable_d = variable.diff(t)
        variable_string = str(variable).replace('\\', '').replace('(t)', '')

        symbol_d = sp.Symbol(f'd_{variable_string}')
        variable_symbol_map[variable_d] = symbol_d
        derivative_symbols.append(symbol_d)

        symbol = sp.Symbol(variable_string)
        variable_symbol_map[variable] = symbol
        plain_symbols.append(symbol)

    numeric_equations = []
    for equation in final_equations:
        if equation != True:
            numeric_equation = equation.copy()

            for variable, symbol in variable_symbol_map.items():
                numeric_equation = sp.Eq(
                    numeric_equation.lhs.subs(variable, symbol),
                    numeric_equation.rhs.subs(variable, symbol)
                )

            numeric_equations.append(numeric_equation)

    save_equations(numeric_equations, 'numeric_system.pdf')

    solved = sp.solve(numeric_equations, derivative_symbols, dict=True)
    solved_dict = {symbol: expression
                   for symbol, expression in sorted(solved[0].items(), key=lambda v: str(v[0]))}
    solved_equations = [sp.Eq(symbol, sp.simplify(expression))
                        for symbol, expression in solved_dict.items()]

    try:
        save_equations(solved_equations, 'solved_numeric_system.pdf')
    except ChildProcessError:
        e.info('could not render the final solved system!')

    equations_map = {sp.pycode(symbol): sp.pycode(expression, fully_qualified_modules=False)
                     for symbol, expression in solved_dict.items()}

    code_template = TEMPLATE_ENV.get_template('system.py.j2')
    code_string = code_template.render({
        'equations_map': equations_map,
        'params_default_map': params_default_map
    })
    code_path = os.path.join(e.path, 'system.py')
    with open(code_path, mode='w') as file:
        file.write(code_string)

    # == MATLAB CODE GENERATION
    matlab_expression_map = {sp.octave_code(symbol).replace('d_', ''): sp.octave_code(expression)
                             for symbol, expression in solved_dict.items()}

    matlab_system_template = TEMPLATE_ENV.get_template('system.m.j2')
    matlab_code = matlab_system_template.render({
        'class_name': 'System',
        'property_value_map': params_default_map,
        'state_expression_map': matlab_expression_map,
        'input_names': ['v_x', 'v_l'],
        'output_names': ['x_out', 'l_out', 'phi_out']
    })
    matlab_code_path = os.path.join(e.path, 'system.m')
    with open(matlab_code_path, mode='w') as file:
        file.write(matlab_code)

