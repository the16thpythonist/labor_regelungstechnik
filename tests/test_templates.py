import os

from labor_regelungstechnik.utils import TEMPLATE_ENV
from .util import ASSETS_PATH, LOG


def test_matlab_system_template():
    template = TEMPLATE_ENV.get_template('system.m.j2')
    content = template.render({
        'class_name': 'StateTest',
        'input_names': ['u1', 'u2'],
        'output_names': ['y1', 'y2'],
        'property_value_map': {
            'm': 1.0,
            'C': -1.0,
        },
        'state_expression_map': {
            'x': 'm * x + C * z',
            'z': 'C * z + u1'
        }
    })

    matlab_path = os.path.join(ASSETS_PATH, 'system.m')
    with open(matlab_path, mode='w') as file:
        file.write(content)


