#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

#requirements = ['py3dmol', 'ipywidgets', 'biopython' ]
requirements = ['py3dmol', 'biopython', 'pandas', 'scikit-learn'] # deps can cause trouble in google colab install

test_requirements = [ ]

setup(
    author="Ajasja Ljubetic",
    author_email='ajasja.ljubetic@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="INSRTR -- Designed allosteric protein logic",
    entry_points={
        'console_scripts': [
            'insrtr=insrtr.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='insrtr',
    name='insrtr',
    packages=find_packages(include=['insrtr', 'insrtr.*']),
    package_data={'insrtr': ['models/*']},
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ajasja/insrtr',
    version='0.1.0',
    zip_safe=False,
)
