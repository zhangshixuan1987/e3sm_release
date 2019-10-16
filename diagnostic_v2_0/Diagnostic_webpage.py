#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make a website of CLUBB outputs

@author: Zhun Guo guozhun@lasg.iap.ac.cn guozhun@uwm.edu
"""

import os
import Common_functions
from shutil import copyfile, move

# Generate the main web page
def main_web(casename,casedir):

    replacements={'CASENAME_here':casename}
    print(replacements.items)
    Common_functions.replace_string("CLUBB_diagnostics_template.html","CLUBB_diagnostics.html",\
                                             replacements)

    move("CLUBB_diagnostics.html",casedir+"/CLUBB_diagnostics.html")

# Generate the webpage for the output files
def sets_web(casename,casedir,varstoplot,setnum,setname,width,height):

    # Replace relevant headers in the webpages
    replacements={'CASENAME_here':casename, 'SETNAME_here':setname}
    Common_functions.replace_string("set_template.htm","set_working.htm",\
                                             replacements)

    with open("set_working.htm","r") as in_file:
        buf = in_file.readlines()

    with open("set_working.htm","w") as out_file:
        for line in buf:
            if line == "<TABLE>\n":
                for v in range(0,len(varstoplot)):
                    thevar=varstoplot[v]
                    line = line + '\n<img src="'+setnum+'/'+thevar+\
                    '.png" style="width:'+width+'px;height:'+height+'px;">'

            out_file.write(line)

    move("set_working.htm",casedir+"/"+setnum+".htm")
