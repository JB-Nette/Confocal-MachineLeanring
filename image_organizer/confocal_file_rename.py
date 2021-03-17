        # -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:19:15 2021

@author: bjeeranan
"""


import xml.etree.ElementTree as et 
import os
import tkinter as tk
from tkinter import filedialog
import shutil
import re


def spaceremoval(path):
    """
    Remove all spaces in file name
    """
    for f in os.listdir(path):       
        os.rename(os.path.join(path,f),os.path.join(path,f.replace(' ','')))
        
   
def imagingrounddetection(path):
    """
    Detect primary patterns in the file in the folder
    Return a list of primary patterns
    """
    # confocal image file name follows a pattern: 1st cycle (name)_F(ii).ims,
    # 2nd cycle (name)_1_F(ii).ims  and so on.  
    #2nd cycles onward image file pattern     
    ims_pattern = re.compile(r"([a-zA-Z0-9_-]+)_(\d+)_F(\d+).(ims)")
    #1st cycle image file pattern
    ims_pattern1 = re.compile(r"([a-zA-Z0-9_-]+)_F(\d+).(ims)")
    
    pattern = []
    pattern1 = []
    # final_pattern shows the primary image file pattern.
    final_pattern = []
    spaceremoval(path)
    for f in os.listdir(path):
        match = re.match(ims_pattern,f)
        match1 = re.match(ims_pattern1,f)
        if match:
            if match.group(1) not in pattern:
                pattern.append(match.group(1))
        if match1:
            if match1.group(1) not in pattern1:
                pattern1.append(match1.group(1))
    if len(pattern) >= 1:
        final_pattern = pattern
    if len(pattern) == 0:
        final_pattern = pattern1

    return final_pattern

def spaceremovalxml(path):
    """
    Remove the space in the file names under IMG_REGEX inside the xml file
    """
    spaceremoval(path)
    pattern = imagingrounddetection(path)
    for pat in pattern:
        xml_pattern = re.compile(r'(' + str(pat) + ').(xml)')
        xml1_pattern = re.compile(r'(' + str(pat) + ')_(\d+).(xml)')
    
        for f in os.listdir(path):
            match = re.match(xml_pattern, f)

            if match:
                xtree = et.parse(os.path.join(path,f))
                xroot = xtree.getroot()
                for node in xroot:
                    if node.tag == 'STACKS':
                        for i in node:
                            d = i.attrib
                            for k,v in d.items():
                                if k == 'IMG_REGEX':
                                    v = v.replace(' ','')
                                    i.attrib[k] = v
                xtree.write(os.path.join(path,f))
        for f in os.listdir(path):
            match1 = re.match(xml1_pattern,f)                
            if match1:
                xtree = et.parse(os.path.join(path,f))
                xroot = xtree.getroot()
                for node in xroot:
                    if node.tag == 'STACKS':
                        for i in node:
                            d = i.attrib
                            for k,v in d.items():
                                if k == 'IMG_REGEX':
                                    v = v.replace(' ','')
                                    i.attrib[k] = v
                xtree.write(os.path.join(path,f))   

def propernaming(path):
    """
    Rename the files of the 1st imaging cycle from 
    1. (name)_F(ii).ims to (name)_0_F(ii).ims 
    2. (name).xml to (name)_0.xml
    3. (name)_metadata.txt to (name)_0_metadata.text
    so that it has the same naming convention as the rest of imaging cycle
    (files with FusionStitcher in name and the .bat files are ignored)
    
    The corresponding ims file name under IMG_REGEX in xml file is renamed as well.
    
    """
    
    pattern = imagingrounddetection(path)
    all_rename = {}
    for pat in pattern:                
        xml_pattern = re.compile(r'(' + str(pat) + ').(xml)')
        ims_pattern = re.compile(r'(' + str(pat) + ')_F(\d+).(ims)')
        txt_pattern = re.compile(r'(' + str(pat) + ')_metadata.(txt)')
        for f in os.listdir(path):
            match = re.match(ims_pattern, f)
            if match:
                new_f = match.group(1) + '_0_F' + match.group(2) + '.' + match.group(3)
                all_rename[f] = new_f
    
        for f in os.listdir(path):   
            match1 = re.match(xml_pattern,f)

            if match1:
                new_f = match1.group(1) + '_0.' + match1.group(2)
                all_rename[f] = new_f
                try:
                    xtree = et.parse(os.path.join(path,f))
                    xroot = xtree.getroot()
                    for node in xroot:
                        if node.tag == 'STACKS':
                            for i in node:
                                d = i.attrib
                                for k,v in d.items():
                                    if k == 'IMG_REGEX':
                                        v = all_rename[v]
                                        i.attrib[k] = v
                    xtree.write(os.path.join(path,f))   
                except KeyError:
                    pass
        for f in os.listdir(path):           
            match2 = re.match(txt_pattern,f)                    
            if match2:
                new_f = match2.group(1) + '_0_metadata.' + match2.group(2)
                all_rename[f] = new_f
        for k,v in all_rename.items():
            os.rename(os.path.join(path,k),os.path.join(path,v))

def findhyb(xml_file):
    """
    

    Parameters
    ----------
    xml_file : The .xml Davil script used for confocal imaging
    Davil script has to set the washing buffer as 'Flow 55percent Buffer'
    
    Returns
    -------
    hyb : Return a list of index containing the bleach cycle. 

    """
    xtree = et.parse(xml_file)
    xroot = xtree.getroot()
    nodes = []
    
    for node in xroot:
        for i,e in enumerate(node):
            n_list = []
            n_list.append(e.tag)
            n_list.append(e.text)
            nodes.append(n_list)
    img_index = []
    hyb_index = []
    for i,e in enumerate(nodes):
        if 'Flow 55percent Buffer' in e:
            hyb_index.append(i)
        if 'DAFusion' in e:
            img_index.append(i)
            
    hyb = []
    base = img_index[0]
    
    for i in range(len(img_index)-1):
        for j in hyb_index:
            if j > base and j < img_index[i+1]:
                hyb.append(i+1)
        base = img_index[i+1]
    return hyb


def namereplacement(xml_file,path,new_path):
    """

    Parameters
    ----------
    xml_file : The .xml Davil script used for confocal imaging
    Davil script has to set the washing buffer as 'Flow 55percent Buffer'
    path : The folder of the images 
    new_path : Path of the newly created folder "Renamed"

    Renamed all bleached cycle according to the hyb_id from findhyb(xml_file) 
    Copy all the renamed file and the unchanged file in the newly created folder "Renamed"

    """
    spaceremovalxml(path)
    propernaming(path)
    hyb_id = findhyb(xml_file)
    pattern = imagingrounddetection(path)
    all_rename = {}
    for pat in pattern:
        xml_pattern = re.compile(r'(' + str(pat) + ')_(\d+).(xml)')
        ims_pattern = re.compile(r'(' + str(pat) + ')_(\d+)_F(\d+).(ims)')
        txt_pattern = re.compile(r'(' + str(pat) + ')_(\d+)_metadata.(txt)')
        for f in os.listdir(path):
           match = re.match(ims_pattern, f)
           match1 = re.match(xml_pattern,f)
           match2 = re.match(txt_pattern,f)
           if match:
               for i,hyb in enumerate(hyb_id):
                    if int(match.group(2)) == int(hyb):
                        new_f = match.group(1) + '_' + match.group(2) + '_bleached_F'+ match.group(3) + '.' + match.group(4)                        
                        all_rename[f] = new_f
                        #print(all_rename)
                    if int(match.group(2)) == int(hyb)-1:
                        new_f = f
                        all_rename[f] = new_f
           # print(all_rename)             
           if match1:
               for i,hyb in enumerate(hyb_id):
                    if int(match1.group(2)) == int(hyb):
                        new_f = match1.group(1) + '_' + match1.group(2) + '_bleached.' + match1.group(3)
                        all_rename[f] = new_f
                    if int(match1.group(2)) == int(hyb)-1:
                        new_f = f
                        all_rename[f] = new_f
           if match2:
               for i,hyb in enumerate(hyb_id):
                    if int(match2.group(2)) == int(hyb):
                        new_f = match2.group(1) + '_' + match2.group(2) + '_bleached_metadata.' + match2.group(3)
                        all_rename[f] = new_f
                    if int(match2.group(2)) == int(hyb)-1:
                        new_f = f
                        all_rename[f] = new_f
        for k,v in all_rename.items():
            shutil.copy(os.path.join(path,k),os.path.join(new_path,v))         
        input('Renaming is done. Press Enter To Exit')
        exit(0)


def main():
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askdirectory(title="Please select image data directory")
    xml_file = filedialog.askopenfilename(title = 'Please choose your Davil script file')
    root.destroy()
    new_path = path +'/Renamed'
    if not os.path.isdir(new_path):
        os.makedirs(new_path)    
    namereplacement(xml_file,path,new_path)



if __name__ == '__main__':
    main()
