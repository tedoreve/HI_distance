#!/usr/bin/python
#Filename: ys_main.py

import  sys
from    astropy.io  import  ascii
import  numpy       as      np
#import spec_fit
#import gaussian
from    astropy.table       import Table, Column
#import spec_outflow
#import spcint_img
import extra_spec
#import spcint_img_1st

obj_list    = ascii.read(sys.argv[1]) #open hii_in_chimps_namelist04.txt
n_obj       = len(obj_list)


CRED = '\033[91m'  ##print word in red
CEND = '\033[0m'
print (CRED+'open the input data'+ sys.argv[1]+CEND)



for i in range(15,n_obj):  #the first one i=0 is out of the data range
    #############################################
    #   Initializing the keywords
    obj_name    = obj_list.field('name')[i]
    lon         = obj_list.field('l_deg')[i]
    lat         = obj_list.field('b_deg')[i]
    source_size = obj_list.field('Reff')[i] #clump size
    ####the following are for readcube function###
    fits_in13    = obj_list.field('chimps_13co')[i] #for 13co data
    fits_in18     = obj_list.field('chimps_c18o')[i] #for c18o data
    fits_atlasgal= obj_list.field('atlasgal_870')[i]

    print (CRED+'start to extra spectra of '+'%d'%i+ obj_name+CEND)
    ####end data preparation for readcube.py#######
    #####the following preparation is for spec_fit.py##
    path13co   = './data/13CO_txt/'
    pathc18o   = './data/C18O_txt/'

    ####data preparation for spcint_img.py##########################
    path13 = '/Volumes/Mac_Passport/Galactic_survey_data/CHIMPS/13CO/'
    path18 = '/Volumes/Mac_Passport/Galactic_survey_data/CHIMPS/C18O/'
    path_atlasgal = '/Volumes/Mac_Passport/Galactic_survey_data/ATLASGAL/'
    ####application for the function#########

    extra_spec.extra_data(  name  =  obj_name,
                      infits=  path13+fits_in13,
                      l_deg =  lon,
                      b_deg =  lat,
                      source_size = source_size,
                      outfile = path13co + obj_name +'_13CO_vall.txt')

    extra_spec.extra_data(  name  =  obj_name,
                            infits=  path18+fits_in18,
                            l_deg =  lon,
                            b_deg =  lat,
                            source_size = source_size,
                            outfile = pathc18o + obj_name +'_C18O_vall.txt')


print(CRED + 'read cube and extract CO spectra successfully' + CEND)



