from setuptools import setup
setup(
  name = 'strain3D',         # How you named your package folder (MyLib)
  packages = ['strain3D'],   # Chose the same as "name"
  version = '1.0.0',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'LV myocardium 3D strain from BSF model',   # Give a short description about your library
  author = 'Yu Zheng',                   # Type in your name
  author_email = 'jorry.zhengyu@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/JorryZ/strain3D',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/JorryZ/strain3D/archive/v1.0.0.tar.gz',    # I explain this later on
  keywords = ['BSF model', 'motion', '3D strain', 'myocardium', 'ECHO'],   # Keywords that define your package best
  install_requires=['numpy','medImgProc','scipy','trimesh','motionSegmentation'],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package    
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',    
    'License :: OSI Approved :: MIT License',   # Again, pick a license    
    'Programming Language :: Python :: 3.6',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.7',
  ],
)
