import os
import pathlib
import subprocess
import tempfile
import typing as t

import shutil
import jinja2 as j2

PATH = pathlib.Path(__file__).parent.absolute()
MEASUREMENTS_PATH = os.path.join(PATH, 'measurements')
TEMPLATES_PATH = os.path.join(PATH, 'templates')




TEMPLATE_ENV = j2.Environment(
    loader=j2.FileSystemLoader(TEMPLATES_PATH),
    autoescape=j2.select_autoescape(),
)
TEMPLATE_ENV.globals.update(**{
    'zip': zip,
    'int': int,
    'enumerate': enumerate,
})

def add_prefix(inp, prefix: str):

    if isinstance(inp, str):
        return prefix + inp
    else:
        return [prefix + string for string in inp]

    return func


TEMPLATE_ENV.filters['add_prefix'] = add_prefix


def latex_math(content: str,
               template_name: str = 'math.tex.j2') -> str:
    template = TEMPLATE_ENV.get_template(template_name)
    latex_string = template.render({'content': content})
    return latex_string


def render_latex(kwargs: dict,
                 output_path: str,
                 template_name: str = 'article.tex.j2'
                 ) -> None:
    with tempfile.TemporaryDirectory() as temp_path:
        # First of all we need to create the latex file on which we can then later invoke "pdflatex"
        template = TEMPLATE_ENV.get_template(template_name)
        latex_string = template.render(**kwargs)
        latex_file_path = os.path.join(temp_path, 'main.tex')
        with open(latex_file_path, mode='w') as file:
            file.write(latex_string)

        # Now we invoke the system "pdflatex" command
        command = (f'pdflatex  '
                   f'-interaction=nonstopmode '
                   f'-output-format=pdf '
                   f'-output-directory={temp_path} '
                   f'{latex_file_path} ')
        proc = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            raise ChildProcessError(f'pdflatex command failed! Maybe pdflatex is not properly installed on '
                                    f'the system? Error: {proc.stdout.decode()}')

        # Now finally we copy the pdf file - currently in the temp folder - to the final destination
        pdf_file_path = os.path.join(temp_path, 'main.pdf')
        shutil.copy(pdf_file_path, output_path)