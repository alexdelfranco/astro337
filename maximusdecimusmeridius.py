# TESTING

def import_all():
  '''
  Test
  '''
  print('hi')
  import numpy as np
  from tqdm import tqdm
  import matplotlib.pyplot as plt
  import pandas as pd
  from astropy.visualization import ZScaleInterval
  from astropy.io import fits 
  import os 
  import glob
  import scipy.ndimage.interpolation as interp


def fits_headerinfo(pathlist,headerTag='IMAGETYP'):
  '''
  Input: A list of fits file paths and a header tag
  Output: A list of information from the specified tag of each file in the path list
  Description: Returns a list of information from the specified tag of each file in the path list
  '''
  info_list = []
  for filename in pathlist:
    header = fits.getheader(filename)
    filetype = header[headerTag]
    info_list.append(filetype)
  return info_list

def fits_filesorter(filename, foldername, fitskeyword_to_check, keyword):
    '''
    Indputs: This function takes four inputs, including a file path to the specified file,
    the name of the folder in which the files are stored, a keyword to specify the expected
    type of fits file, and a second keyword to check that the information in the fits
    file matches expectation
    Outputs: Returns nothing.
    Description: Checks to see if the specified fits file matches specific properties specified in the
    inputs. Moves the file into a new folder, creates the folder if it doesn't already exist
    '''
    # Check to see if the file exists, if it doesn't then print an error
    # and return nothing. If it does exist, proceed
    if os.path.exists(filename):
        pass
    else:
        # print(filename + " does not exist or has already been moved.")
        return
    
    header = fits.getheader(filename)
    fits_type = header[keyword]
    
    # Check to see if the folder exists, if it doesn't then make the directory
    # print what is being done. If it does exist, proceed
    if os.path.exists(foldername):
        pass
    else:
        if '/' not in foldername:
          print("Making new directory: " + foldername)
          os.mkdir(foldername)
        elif '/' in foldername:
          print("Making new directories: " + foldername)
          os.makedirs(foldername)
    
    # Check to see if the fits file matches our predicted input file type
    # Create a destination string by concatenating the foldername, filename,
    # and the file path
    if fits_type == fitskeyword_to_check:
        destination = foldername + '/'
        # print("Moving " + filename + " to: ./" + destination + filename)
        os.rename(filename, destination + filename)  
    return

def organize_calibration():
  '''
  Input: None
  Output: None
  Description: Organizes calibration data into separate folders within a calibration folder
  '''
  fits_pathlist = glob.glob('*.fit')
  # Sort calibration data into separate folders:
  for fitsfile in tqdm(fits_pathlist):
      fits_filesorter(fitsfile, 'calibration/flats', 'Flat Field', 'IMAGETYP')
      fits_filesorter(fitsfile, 'calibration/biasframes', 'Bias Frame', 'IMAGETYP')
      fits_filesorter(fitsfile, 'calibration/darks', 'Dark Frame', 'IMAGETYP')

def mediancombine(filelist):
    '''
    Input: A list of file paths for fits files
    Output: A single fits frame that contains the median pixel valeues of all the
    input files
    Description: Combines an array of fits images into a single median image
    '''
    # Find the length of the filelist input (the number of fits files)
    n = len(filelist)
    
    # Take the data from the first file in the list
    first_frame_data = fits.getdata(filelist[0])
    
    # Take the dimensions of the first fits data matrix
    imsize_y, imsize_x = first_frame_data.shape
    
    # Create a 3d matrix of zeros the size of all the fits file data
    fits_stack = np.zeros((imsize_y, imsize_x , n)) 
    
    # Loop through each fits file and add its data to the matrix cube
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        fits_stack[:,:,ii] = im
        
    # Collapse the 3d matrix into a 2d array of median values      
    med_frame = np.median(fits_stack, axis = 2)
    
    return med_frame

def subtract(filename, path_to_bias, prefix, write=False):
    '''
    Input: A fits filename and the path to our master bias
    Output: No output, but writes a fits file with the prefix 'b_'
    Description: Subtracts the data from the master bias fits file from an input file
    '''
    # Get file data and header
    file_data = fits.getdata(filename)
    file_header = fits.getheader(filename)
    # Get bias data
    bias_data = fits.getdata(path_to_bias)
    # Subtract bias from the file data
    bsub_data = file_data - bias_data
    # If we want to write to the local directory
    if write:
      filepath = filename.split('/')
      filepath[-1] = prefix+filepath[-1]
      filename = os.path.join(*filepath)
      fits.writeto(filename, bsub_data, file_header, overwrite=True)
    # Return the bias subtracted data
    return bsub_data

def norm_combine_flats(filelist):
    '''
    Takes in a filelist, normalizes each image and then takes the median value from each image for each pixel 
    returns a new image that is the median value for each pixel
    '''
    # Number of files in filelist
    n = len(filelist)
    
    # Gets the data from the first file in file list
    first_frame_data = fits.getdata(filelist[0])
    
    # Assignes the shape of the two dimensional array from fits to new variables (y and x) for new array
    imsize_y, imsize_x = first_frame_data.shape
    
    # Creates a two-dimensional array of the shape above that is all zeros n times
    fits_stack = np.zeros((imsize_y, imsize_x , n))
    
    # runs through the range of images in filelist, pulls each ones data, normalizes using by dividing by median and then stores the data in fits_stack
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        norm_im =  im / np.median(im) # finish new line here to normalize flats
        fits_stack[:,:,ii] = norm_im
        
    # Takes the median for each pixel from the values of the given pixel in each image. Creates med_frame which is an image of each median value for each pixel
    med_frame = np.median(fits_stack, axis=2)
    return med_frame

def darksort():
  '''
  Input: None
  Output: A list of the unique exposure times as strings
  Description: Sorts the dark files by exposure time into separate folders
  '''
  os.chdir('calibration/darks/')
  darks = glob.glob('*.fit')
  # Find the exposure times of the darks
  dark_exp = list(set(fits_headerinfo(darks,'EXPTIME')))
  dark_exp_str = []
  for expTime in dark_exp:
    dark_exp_str.append(str(int(expTime)) + 'sec')
  for fitsfile in darks:
    for i in range(len(dark_exp)):
        fits_filesorter(fitsfile, dark_exp_str[i], dark_exp[i], 'EXPTIME')
  os.chdir('../..')
  print(dark_exp_str)
  return dark_exp_str

def flatsort():
  '''
  Input: None
  Output: None
  Description: Sorts the flats by band into separate folders
  '''
  os.chdir('calibration/flats')
  flats = glob.glob('*.fit')
  # sort the flats into their respective subfolders (feel free to add cells as needed)
  for fitsfile in flats:
    fits_filesorter(fitsfile, 'Rband', 'Red', 'FILTER')
    fits_filesorter(fitsfile, 'Vband', 'Visual', 'FILTER')
    fits_filesorter(fitsfile, 'Bband', 'Blue', 'FILTER')
  os.chdir('../..')

def master_bias():
  '''
  Input: None
  Output: None
  Description: Create and write a master bias
  '''
  os.chdir('calibration/biasframes/')
  biasfiles = glob.glob('*.fit')
  # Take the median of the 3-D stack: 
  med_bias = mediancombine(biasfiles)
  # Save your new FITS file!
  fits.writeto('Master_Bias.fit', med_bias, fits.getheader(biasfiles[0]), overwrite=True)
  master_bias_path = os.getcwd() + '/Master_Bias.fit'
  os.chdir('../..')

def master_darks(dark_exp_str):
  '''
  Input: A list of exposure times as strings
  Output: None
  Description: Create and write a master dark for each exposure time
  '''
  master_bias_path = 'calibration/biasframes/Master_Bias.fit'
  for expTime in tqdm(dark_exp_str):
    darkfiles = glob.glob('calibration/darks/'+expTime+'/*.fit')
    bs_dark_data = []
    for filepath in darkfiles:
      bs = subtract(filepath, master_bias_path, 'b_', True)
      bs_dark_data.append(bs)

    x,y = bs_dark_data[0].shape
    fits_stack = np.zeros((x,y,len(bs_dark_data)))
    for ii in range(len(bs_dark_data)):
      fits_stack[:,:,ii] = bs_dark_data[ii]

    dark_header = fits.getheader(darkfiles[0])
    master_dark_name = 'calibration/darks/'+expTime+'/Master_Dark_'+expTime+'.fit'
    master_dark = np.median(fits_stack,axis=2)

    fits.writeto(master_dark_name, master_dark, dark_header, overwrite=True)

def master_darkpath(dark_exp_str):
  '''
  Input: None
  Output: A dictionary of paths to master darks
  Description: Create a dictionary of paths to master darks of different exposure times
  '''
  master_dark_path = {}
  for i in range(len(dark_exp_str)):
    master_dark_path[dark_exp_str[i]] = 'calibration/darks/'+dark_exp_str[i]+'/Master_Dark_'+dark_exp_str[i]+'.fit'
  return master_dark_path

def master_flats(dark_exp_str):
  '''
  Input: None
  Output: None
  Description: Create and write master flats for each band
  '''
  master_bias_path = 'calibration/biasframes/Master_Bias.fit'
  master_dark_path = master_darkpath(dark_exp_str)
  for band in tqdm(['Bband','Rband','Vband']):
    flatfiles = glob.glob('calibration/flats/'+band+'/*.fit')
    # Create a list for the bias subtracted data
    bs_flat_data = []
    for filepath in flatfiles:
      bs = subtract(filepath, master_bias_path, 'b_', True)
      bs_flat_data.append(bs)

    bs_flatfiles = glob.glob('calibration/flats/'+band+'/b_*.fit')
    dbs_flat_data = []
    for filepath in bs_flatfiles:
      # Figure out which master dark to use
      expTime = fits.getheader(filepath)['EXPTIME']
      dbs = subtract(filepath, master_dark_path[str(int(float(expTime)))+'sec'], 'd', True)
      dbs_flat_data.append(dbs)

    # Prepare to median combine the data into one file
    db_flats = glob.glob('calibration/flats/'+band+'/db_*.fit')
    master_norm = norm_combine_flats(db_flats)
    flat_header = fits.getheader(db_flats[0])
    master_flat_name = 'calibration/flats/'+band+'/Master_'+band+'_Flat.fit'

    fits.writeto(master_flat_name, master_norm, flat_header, overwrite=True)

def keyword_sort(filename, foldername, keyword):
    '''
    Indputs: This function takes four inputs, including a file path to the specified file,
    the name of the folder in which the files are stored, a keyword to specify the expected
    type of fits file, and a second keyword to check that the information in the fits
    file matches expectation
    Outputs: Returns nothing.
    Description: Checks to see if the specified fits file matches specific properties specified in the
    inputs. Moves the file into a new folder, creates the folder if it doesn't already exist
    '''
    # Check to see if the file exists, if it doesn't then print an error
    # and return nothing. If it does exist, proceed
    if os.path.exists(filename):
        pass
    else:
        # print(filename + " does not exist or has already been moved.")
        return
    
    # Check to see if the folder exists, if it doesn't then make the directory
    # print what is being done. If it does exist, proceed
    if os.path.exists(foldername):
        pass
    else:
        if '/' not in foldername:
          print("Making new directory: " + foldername)
          os.mkdir(foldername)
        elif '/' in foldername:
          print("Making new directories: " + foldername)
          os.makedirs(foldername)
    
    # Check to see if the fits file matches our predicted input file type
    # Create a destination string by concatenating the foldername, filename,
    # and the file path
    if keyword in filename:
        destination = foldername + '/'
        # print("Moving " + filename + " to: ./" + destination + filename)
        os.rename(filename, destination + filename)  
    return

def sort_extras(img_numlist):
  '''
  Input: An array of fits image numbers of science images
  Output: A list of paths to the science images
  Description: Sorts the science images and returns their paths while
    sorting the nonscience images into an 'extra' folder
  '''
  # Retrieve all file paths from the chosen directory
  object_full_list = glob.glob('*.fit')
  object_paths = []
  # Select the preferred image numbers and save them to object_paths
  for img_num in img_numlist:
    for path in object_full_list:
      if str(img_num) in path:
        object_paths.append(path)
  # Move all the rest of the images into an "extra" folder
  for path in object_full_list:
    if path not in object_paths:
      foldername = 'extra/'
      # If the extra folder doesn't exist, make it
      if os.path.exists(foldername):
          pass
      else:
        print("Making new directory: " + foldername)
        os.mkdir(foldername)
      # Move the file to the extra folder
      os.rename(path, foldername + path)
  return object_paths

def sortraw(folder_names,file_strings):
  sc_images = glob.glob('*.fit')
  for img_path in sc_images:
    for index,folder in enumerate(folder_names):
      keyword_sort(img_path,folder,file_strings[index])

def find_exptime(image_range):
  '''
  Input: None
  Output: The single exposure time of the set of science images
  Description: Finds the exposure time of a set of images in the local folder
    and returns an exposure time if they all have the same one, or prints a statement
    if they include multiple exposure times
  '''
  object_paths = sort_extras(image_range)
  expTime = list(set(fits_headerinfo(object_paths,'EXPTIME')))
  if len(expTime) == 1:
    expTime = float(expTime[0])
    if expTime%1 == 0: expTime = int(expTime)
    print('The exposure time is '+str(expTime)+' seconds')
    return expTime
  else:
    print('This folder has images at multiple exposure times. Please choose one manually:')
    print(expTime)
    return 'null'

def db_subtract(expTime):
  '''
  Input: None
  Output: None
  Description: Bias and dark subtracts a list of images in the current directory and
    writes the new fits files
  '''
  master_bias_path = '../calibration/biasframes/Master_Bias.fit'
  master_dark_path = '../calibration/darks/'+str(expTime)+'sec/Master_Dark_'+str(expTime)+'sec.fit'
  
  object_rawpaths = glob.glob('*.fit')
  print('Bias Subtracting:')
  for filepath in tqdm(object_rawpaths):
    subtract(filepath, master_bias_path, 'b_', True)
  object_bpaths = glob.glob('b_*.fit')
  print('Dark Subtracting:')
  for filepath in tqdm(object_bpaths):
    subtract(filepath, master_dark_path, 'd', True)
  object_dbpaths = glob.glob('db_*.fit')

  object_paths = glob.glob('*.fit')
  print('Sorting')
  for img_path in object_paths:
    keyword_sort(img_path, 'db_sub', 'db_')
    keyword_sort(img_path, 'b_sub', 'b_')
    keyword_sort(img_path, 'raw', '')

def flat_fielding():
  '''
  Input: None
  Output: None
  Description: Flat fields science images and writes them as new fits files
  '''
  # Create a dictionary of master flats by band
  mflats = {}
  for band in ['Bband','Vband','Rband']:
    master_flat = fits.getdata('../calibration/flats/'+band+'/Master_'+band+'_Flat.fit')
    mflats[band] = master_flat
  db_sub = glob.glob('db_sub/*.fit')
  # For each path, divide it by its master flat
  print('Flat Fielding:')
  for path in tqdm(db_sub):
    data = fits.getdata(path)
    header = fits.getheader(path)
    band = header['Filter']
    div_data = data / mflats[band[0]+'band']
    # Create a name for the new image
    plist = path.split('/')
    plist[-2] = 'reduced'
    plist[-1] = 'f'+plist[-1]
    # foldername = os.path.join(*plist[])+'/reduced'
    # Create a new folder if necessary
    if os.path.exists('reduced'):
        pass
    else:
      print("Making new directory: " + 'reduced')
      os.mkdir('reduced')
    # Write the new flat fielded image
    fits.writeto(os.path.join(*plist), div_data, header, overwrite=True)

def centriod_math(image):
  '''
  Input: An image matrix
  Output: An array of the centroid coordinates of the main object in the image
  Description: Computes double sums to find the centroid of a matrix
  '''
  xsum,ysum,dsum = 0,0,0

  for xindex in range(image.shape[1]):
    for yindex in range(image.shape[0]):
      xsum += xindex*image[yindex,xindex]
      ysum += yindex*image[yindex,xindex]
  dsum = np.nansum(np.nansum(image,axis=0))

  return(np.array([xsum/dsum,ysum/dsum]))

def calc_centroid(im1,im2):
  '''
  Input: Two image arrays
  Output: The x and y offsets between the two images
  Description: Returns the centroid offset between two image arrays
  '''
  im1_grey = im1.astype('float')
  im2_grey = im2.astype('float')

  im1_cent = centriod_math(im1_grey)
  im2_cent = centriod_math(im2_grey)

  offset = im1_cent - im2_cent

  return offset[0],offset[1]

def register(imlist,region):
  '''
  Input: A list of images and a region to examine for centroiding
  Output: A list of shifts relative to the first image
  Description: Returns the shifts between the first image and the rest in a series of images
  '''
  xcrop = region[0]
  ycrop = region[1]
  
  im1,hdr1 = fits.getdata(imlist[0],header=True)

  im1 = processing(im1,region)
  shifts = [(0,0)]
  for path in imlist[1:]:
    im,hdr = fits.getdata(path,header=True)
    im = processing(im,region)
    shifts.append(calc_centroid(im1,im))
  
  return shifts

def sclip(im,times):
  '''
  Input: An image matrix and the number of times to run this function
  Output: An image matrix with a certain number of values over three standard deviations
  from the median clipped out
  Description: Sigma clips an array a specific number of times
  '''
  if times == 0: return im
  bkg = np.where(im > (np.nanmedian(im)+3*np.nanstd(im)), np.nan, im)
  return(sclip(bkg,times-1))

def processing(im,region):
  '''
  Input: An image matrix and a region within that image
  Output: A matrix of of values that deviate from the image more than 3 standard deviations
  Description: Crops, cleans, and removes the background from a star image
  '''
  xcrop,ycrop = region[0],region[1]
  im = im.astype('float')
  im[~np.isfinite(im)] = np.nan
  im = im[ycrop[0]:ycrop[1],xcrop[0]:xcrop[1]]
  bkg = sclip(im,5)
  bstd = np.nanstd(bkg)
  dev = np.where(im > (np.nanmean(bkg) + 3*bstd), im, 0)
  return dev

def implot(im1):
  '''
  Input: An image matrix
  Output: None
  Description: Plots an image the way I like
  '''
  fig,ax = plt.subplots(figsize=(20,10))
  interval = ZScaleInterval() # This decides on an intelligent scaling for the image automatically
  vmin, vmax = interval.get_limits(im1) # Use the data to set min and max for the scale
  im = ax.imshow(im1, vmin=vmin, vmax=vmax)

def centroid(reduced_paths,xyrange):
  '''
  Input: A list of reduced fits file paths of science images
  Output: The necessary shifts in those images relative to the first
  Description: This function uses centroiding to calculate the necessary shifts to align
    a list of images with the first image in the list
  '''
  imlist = []
  imrange = np.append(np.arange(67,72),(63))
  for path in reduced_paths:
    for num in imrange:
      if '000000'+str(num) in path:
        imlist.append(path)
  imlist = sorted(imlist)
  return register(imlist,xyrange)

def shift_image(image,xshift,yshift):
    '''
    Input: An image matrix and the necessary x and y pixel shifts
    Output: A shifted matrix array
    Description: Shift a matrix array by rolling pixels over the x and y axes
    '''
    # Note that this will not do any trimming, 
    # so we'll want to  trim later the edges of the image using the maximum shift.
    return np.roll(np.roll(image,int(yshift),axis=1), int(xshift), axis=0)

def centcrop(imsize,radius):
  '''
  Input: The image diameter and its radius of the center crop
  Output: A list of lists containting the x and y bounds of the crop
  Description: Find the bounds for a crop with a specific radius and image size
  '''
  imrad = imsize/2
  return [[int(imrad-radius),int(imrad+radius)],[int(imrad-radius),int(imrad+radius)]]

def regstack(datadir,imlist,offsets,padsize,objname='null'):
  '''
  Input: A list of fits images, their offsets, a padding size, and an optional object name for writing filepaths
  Output: A registered stack of images from the image list
  Description: A function that registers and combines a stack of fits images and writes them to a file
  '''
  im1,hdr1 = fits.getdata(imlist[0],header=True)
  image_stack = np.zeros([im1.shape[0]+2*padsize,im1.shape[1]+2*padsize,len(imlist)])
  if os.path.isdir(datadir + '/registered') == False: os.mkdir(datadir + '/registered')
  if os.path.isdir(datadir + '/stacked') == False: os.mkdir(datadir + '/stacked')
  for index,filename in tqdm(enumerate(imlist)):
    im,hdr = fits.getdata(filename,header=True)
    padded = np.pad(im,padsize,'constant', constant_values=-0.001)
    shifted = interp.shift(padded,(offsets[index][1],offsets[index][0]), cval=-0.001)
    shifted[shifted <= -0.0001] = np.nan
    image_stack[:,:,index] = shifted
    fits.writeto(datadir + '/registered/'+'reg_'+filename.split('/')[-1],shifted,hdr1,overwrite=True)
  median_image = np.nanmedian(image_stack,axis=2)
  if objname == 'null': objname = hdr1['OBJECT']
  fits.writeto(datadir + '/stacked/' + objname + '_' + hdr1['FILTER'] + 'Stack.fits', median_image, overwrite=True)
  return median_image

def getfiltlist(filtername,path):
  '''
  Input: A filter band and a path to a folder containting fits files
  Output: A list of fits files imaged in that filter (band)
  Description: Find which files in a folder were imaged in a chosen band
  '''
  imlist = glob.glob(path)
  filtlist = []
  for fitsfile in imlist:
    header = fits.getheader(fitsfile)
    filter = header['FILTER']
    if filtername == filter:
      filtlist.append(fitsfile)
  return filtlist

def masters():
  '''
  Input: None
  Output: None
  Description: Creates and sorts all necessary calibration frames
  '''
  import_all()
  import glob
  organize_calibration()
  expTimes = darksort()
  flatsort()
  master_bias()
  master_darks(expTimes)
  master_flats(expTimes)

def calibrate(targ_dict):
  sortraw(targ_dict['names'],targ_dict['fstrings'])
  for index in targ_dict['targ_folders']:
    folder = targ_dict['names'][index]
    os.chdir(folder)
    object_paths = sort_extras(targ_dict['imrange'][index])
    expTime = find_exptime(targ_dict['imrange'][index])
    if expTime == 'null': return
    db_subtract(expTime)
    flat_fielding()
    os.chdir('..')

def full_register(datadir,band,region=centcrop(4096,50),objname='null'):
  filtlist = getfiltlist(band,datadir+'/reduced/*.fit')
  offsets = register(filtlist,region)
  regtest = regstack(datadir,filtlist,offsets,25,objname)
