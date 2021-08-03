from skbuild import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="VolcGases",
    packages=['VolcGases'],
    version='2.3',
    license='MIT',
    install_requires=['numpy','numba','scipy'],
    author = 'Nicholas Wogan',
    author_email = 'nicholaswogan@gmail.com',
    description = 'Python program that calculates the '+\
                  'gases produced by a volcano.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>3.6',
    url = "https://github.com/Nicholaswogan/VolcGases",
    cmake_args=['-DSKBUILD=ON']
    )