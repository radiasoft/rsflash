# -*- coding: utf-8 -*-
u"""rsflash setup script

:copyright: Copyright (c) 2021 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from pykern import pksetup

pksetup.setup(
    name='rsflash',
    author='RadiaSoft LLC',
    author_email='pip@radiasoft.net',
    description='Tools for performing and analyzing FLASH simulations of capillary discharge plasmas',
    install_requires=[
        'pykern',
    ],
    license='http://www.apache.org/licenses/LICENSE-2.0.html',
    url='https://github.com/radiasoft/rsflash',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python',
        'Topic :: Utilities',
    ],
)
