import os
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from glob import glob
# from distutils.command.install import install
import platform

system = platform.system()

libdir = os.path.dirname(os.path.abspath(__file__))
# cfitslib = os.path.join(libdir, 'csst_mci_psf', 'lib')
# cfitsinc = os.path.join(libdir, 'csst_mci_psf', 'include')
psflib = os.path.join(libdir, 'HybPSF')
# Define the extension module


if system == 'Linux':
    shared_ext = 'so'
    psffit_library_name = 'libpsffit.so'
    star_extrac_library_name = 'star_extrac.so'
elif system == 'Darwin':
    shared_ext = 'dylib'
    psffit_library_name = 'libpsffit.dylib'
    star_extrac_library_name = 'star_extrac.dylib'
else:
    raise ValueError(f'Unsupported platform: {system}')

try:
    os.mkdir(os.path.join(os.path.expanduser('~'), '.cfitsio'))
except OSError:
    pass

prefix = '--prefix='+os.path.join(os.path.expanduser('~'), '.cfitsio')
prefixdir = os.path.join(os.path.expanduser('~'), '.cfitsio')


with open("requirements.txt", "r") as f:
    requirements = [
        req.strip()
        for req in f.readlines()
        if not req.startswith("#") and req.__contains__("==")
    ]


def match(fname, str1):
    if str1 in open(fname).read():
        return (True)
    else:
        return (False)


if system == 'Linux':
    bashrc = os.path.join(os.path.expanduser('~'), '.bashrc')
    str1 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+os.path.join(os.path.expanduser('~'), '.cfitsio', 'lib\n')
    if match(bashrc, str1) is not True:
        with open(bashrc, "a") as file:
            print(file.write(str1))
            file.write(str1)
        file.close()
    else:
        pass
    subprocess.run(['source', '~/.bashrc'], shell=True)
elif system == 'Darwin':
    bashrc = os.path.join(os.path.expanduser('~'), '.bash_profile')
    str1 = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+os.path.join(os.path.expanduser('~'), '.cfitsio', 'lib\n')
    if match(bashrc, str1) is not True:
        with open(bashrc, "a") as file:
            print(str1)
            file.write(str1)
        file.close()
    else:
        pass
    subprocess.run(['source', '~/.bash_profile'], shell=True)


class BuildCfitsio(build_ext):
    def run(self):
        # Run the installation process for cfitsio library
        self.install_cfitsio()

        # Call the parent run() method to build the extensions
        super().run()

    def install_cfitsio(self):
        # Run the installation commands for cfitsio library
        subprocess.check_call(['./configure',
                               prefix], cwd=os.path.join(psflib, 'cfitsio'))
        # print(['configure', prefix])
        subprocess.check_call(['make'], cwd=os.path.join(psflib, 'cfitsio'))
        subprocess.check_call(['make', 'install'], cwd=os.path.join(psflib, 'cfitsio'))


extensions = [
    Extension('HybPSF.source.libpsffit',
              sources=[os.path.join('HybPSF', 'source',
                                    os.path.basename(name)) for name in glob('HybPSF/source/*.c')],
              include_dirs=[os.path.join(prefixdir, 'include'), '/usr/local/include'],
              library_dirs=[os.path.join(prefixdir, 'lib'), '/usr/local/lib'],
              libraries=['cfitsio'],
              # extra_link_args=[f'-install_name @rpath/{psffit_library_name}'],
              extra_compile_args=['-fPIC'],
              language="c",
              # extension_name=f'libpsffit.{shared_ext}',
              ),
    Extension('HybPSF.staridf.star_extrac',
              sources=[os.path.join('HybPSF', 'staridf',
                                    os.path.basename(name)) for name in glob('HybPSF/staridf/*.c')],
              include_dirs=[os.path.join(prefixdir, 'include'), '/usr/local/include'],
              library_dirs=[os.path.join(prefixdir, 'lib'), '/usr/local/lib'],
              libraries=['cfitsio'],
              extra_compile_args=['-fPIC'],
              language="c",
              # extension_name=f'star_extrac.{shared_ext}',
              ),
]


VERSION = '0.0.1' 
DESCRIPTION = 'PSF module of JWST NIRCam image '
LONG_DESCRIPTION = 'create Point Spread Function (PSF) module for the JWST NIRCam image image'


setup(
    name="HybPSF",
    version="1.0.0",
    author="JWST_revised team @shao",
    author_email="linn@mail.ustc.edu.cn",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license='MIT',
    classifiers= [
        "Development Status :: beta",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    cmdclass={
        'build_ext': BuildCfitsio
    },
    ext_modules=extensions,
    packages=find_packages(),
    package_dir={'HybPSF':'./HybPSF',
                     'cfitsio':'./cfitsio',
                     'config':'./config',
                     'source':'./source',
                     'staridf':'./staridf'},
    package_data={
        '': ['./*.so'],
        'HybPSF': ['config/*', 'include/*', 'lib/*', 'source/*', 'staridf/*'],
    }

)