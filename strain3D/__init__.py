# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:56:00 2019

@author: JorryZ

File: __init__.py
Description: load strainMyocardium

History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: jorry.zhengyu@gmail.com         20NOV2019             -V1.0.0 test version
                                                        -strainMyocardium version 1.0.0  
  Author: jorry.zhengyu@gmail.com         07Dec2019             -V1.0.1 test version
                                                        -strainMyocardium version 1.0.1                                                        
  Author: jorry.zhengyu@gmail.com         03JAN2020             -V1.0.2 test version
                                                        -strainMyocardium version 1.0.2
  Author: renmeifeng05@gmail.com          19JUN2020             -V1.0.3 test version
                                                        -strainMyocardium version 1.0.3
  Author: jorry.zhengyu@gmail.com         25JUN2020             -V1.0.4 test version
                                                        -strainMyocardium version 1.0.4   
  Author: jorry.zhengyu@gmail.com         09July2020            -V1.0.5 test version
                                                        -strainMyocardium version 1.0.5  
  Author: jorry.zhengyu@gmail.com         09July2020            -V1.0.6 test version
                                                        -strainMyocardium version 1.0.6    
  Author: jorry.zhengyu@gmail.com         09July2020            -V1.0.7 test version
                                                        -strainMyocardium version 1.0.7 
  Author: jorry.zhengyu@gmail.com         24July2021            -V1.0.8 test version
                                                        -strainMyocardium version 1.0.8                                                         
  Author: jorry.zhengyu@gmail.com         06Dec2021             -V1.0.9 test version
                                                        -strainMyocardium version 1.0.9
 
Requirements:
    numpy
    scipy
    motionSegmentation
    medImgProc
    trimesh
    json (optional)
    pickle (optional)
All rights reserved.
"""
_version='1.0.9'
print('strain3D version',_version)

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import strain3D.strainMyocardium as strainMyocardium
