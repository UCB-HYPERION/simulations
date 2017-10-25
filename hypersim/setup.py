from setuptools import setup
import glob

setup(name='hypersim',
        version='0.0.1',
        description='HYPERION simulation tools',
        author='Kara Kundert',
        package_dir={'hypersim':'src'},
        packages=['hypersim'],
        scripts=glob.glob('scripts/*')
)
