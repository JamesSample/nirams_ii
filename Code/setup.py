from distutils.core import setup
	  
setup(name='niramsii',
      version=0.1,
	packages=['niramsii',],
      author='James Sample',
      author_email='james.e.sample@gmail.com',
      url='https://github.com/JamesSample/nirams_ii/',
      description='The Nitrogen Risk Assessment Model for Scotland (version 2; NIRAMS II).',
	long_description=open('README.txt').read(),
      package_data={'': ['LICENSE.txt']},
      include_package_data=True,
      classifiers=['Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: MIT License',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python'])