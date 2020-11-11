# main.py
# import logging
from image_moments import CreateImages
from create_moments import ClipCube
import warnings; warnings.filterwarnings("ignore")
import os, sys, glob
from matplotlib import pyplot as plt
from pathlib import Path
from datetime import datetime

class Logger(object):
    """A class to duplicate an output stream to stdout/err.
       This works in a manner very similar to the Unix ‘tee’ command."""
    # source: https://stackoverflow.com/questions/616645/how-to-duplicate-sys-stdout-to-a-log-file
    def __init__(self, filename, mode):
        self.stdout = sys.stdout
        self.file = open(filename, mode)
        sys.stdout = self

    def __del__(self):
        self.close()

    def __enter__(self):
        pass

    def __exit__(self, *args):
        self.close()

    def write(self, message):
        self.stdout.write(message)
        # self.stderr.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        # os.fsync(self.file.fileno())

    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None

        if self.file != None:
            self.file.close()
            self.file = None

class dummy_context_mgr():
    # Context manager that returns enter_result from __enter__, but otherwise does nothing
    def __enter__(self):
        return None
    def __exit__(self, exc_type, exc_value, traceback):
        return False


def main(logOutput=False):



    if logOutput:
        logfile = 'pipeline-'+'{:{dfmt}-{tfmt}}'.format(datetime.now(), dfmt='%Y%m%d', tfmt='%H:%M:%S') + '.log'
        print('')
    else:
        logfile = ""

    with Logger(logfile, 'w') if logOutput else dummy_context_mgr(): # execute either with or without logging, set by logOutput

        print(logfile)
        print("")
        starttime = datetime.now()
        print('START: {:{dfmt} {tfmt}}'.format(starttime, dfmt='%Y-%m-%d', tfmt='%H:%M:%S'))
        print("")
    
        ######################################################################################
        # User settings ######################################################################
        # Applied to all routines, otherwise use defaults
        
        # data release version to be used (and included file name, header information)
        version = '1.0'

        # VERTICO datacubes are produced with different round beam sizes:
        # 'native' = native cube resolution_str (typically 7" - 9")
        # '9arcsec' = 9" beam cubes
        # '15arcsec' = 15" beam cubes
        resolution_str = 'native'

        # Select the products you want the pipeline to produce
        # [ADD DESC HERE]
        # Product list:
        #     mom0_brightness - 
        #     mom0_density - 
        #     mom1 - 
        #     mom2 - 
        #     peakT - 
        #     mom0_noise - 
        #     mom1_noise - 
        #     pvd_major - 
        #     pvd_minor - 
        #     spectrum_vel - 
        #     spectrum_veloffset - 
        #     spectrum_freq - 
        #     radial_profile_brightness_arcsec - 
        #     radial_profile_brightness_kpc - 
        #     radial_profile_density_arcsec - 
        #     radial_profile_density_kpc - 

        productList = ['mom0', 
                        'mom1', 
                        'mom2', 
                        'peakT', 
                        'mom0_noise', 
                        'mom1_noise', 
                        'pvd', 
                        'spectrum'
                        'radial_profile']

        refresh = True # TB - I don't know what this does?
        overwrite = True # Overwrite when writing products and figures to file
        sunflag = True  # Compute the cube masks using method outlined in Sun+2018 (PHANGS)
        clip = True # Clip the input data cubes in velocity space
        tosave = True # Save products and figures to file
        pbcor = True # Use primary beam corrected cube
        TP = True # If True, products will use 7m+TP datacube where available. If False, 7m only cubes will be used
        alpha_co = 6.25

        # Summary
        print("USER SETTINGS:")
        print("data version", version)
        print("product list: " + ','.join(productList))
        print("refresh = ", refresh)
        print("overwrite = ", overwrite)
        print("sunflag = ", sunflag)
        print("clip = ", clip)
        print("tosave = ", tosave)
        print("pbcor = ", pbcor)
        print("resolution_str", resolution_str)
        print("alpha_co = {} Msol/pc^2 K km/s".format(alpha_co))
        print("")

        ######################################################################################
        # Local paths ########################################################################
        # The product pipeline requires two paths:
        # 1. inputDataPath, input directory containing processed data cubes in the following directory structure
        #         inputDataPath/[resolution_str]/[galaxyID]/
        # (note the resolution and galaxyID are set programmatically so not included in this path
        # 2. outputProductPath, output directory where products will be stored (if specified)

        inputDataPath = '/Users/thbrown/VERTICO/share/cubes/'

        outputProductPath = '/Users/thbrown/VERTICO/pipelineTesting/'

        # in/output path objects
        readpath = Path(inputDataPath + 'v'+version + '/' + resolution_str + '/')

        writepath = Path(outputProductPath + '/' + 'products.' + 'v'+version + '/' + resolution_str + '/')

        ######################################################################################
        # Input galaxy ID(s) #################################################################

        # all galaxies for which processed data exist
        galaxies = [dI for dI in os.listdir(readpath) if os.path.isdir(os.path.join(readpath,dI)) and ('_' not in dI)]  

        # galaxies = ['IC3392', 'NGC4064', 'NGC4189', 'NGC4192', 'NGC4216', 'NGC4222', 'NGC4294', 'NGC4299', 'NGC4302',
        #             'NGC4330', 'NGC4351', 'NGC4380', 'NGC4383', 'NGC4388', 'NGC4394', 'NGC4405', 'NGC4419', 'NGC4522',
        #             'NGC4532', 'NGC4533', 'NGC4568', 'NGC4606', 'NGC4607', 'NGC4651', 'NGC4713', 'NGC4808', 'NGC4396',
        #             'NGC4567', 'NGC4772', 'NGC4580', 'NGC4450', 'NGC4254', 'NGC4293', 'NGC4298', 'NGC4321', 'NGC4402',
        #             'NGC4424', 'NGC4457', 'NGC4535', 'NGC4536', 'NGC4548', 'NGC4569', 'NGC4579', 'NGC4654', 'NGC4689',
        #             'NGC4698', 'NGC4694']

        #galaxies = ['IC3392', 'NGC4064', 'NGC4189', 'NGC4192', 'NGC4216', 'NGC4222', 'NGC4294', 'NGC4299', 'NGC4302',
        #            'NGC4330', 'NGC4351', 'NGC4380', 'NGC4383', 'NGC4388', 'NGC4394', 'NGC4405', 'NGC4419', 'NGC4522',
        #            'NGC4532', 'NGC4533', 'NGC4568', 'NGC4606', 'NGC4607', 'NGC4651', 'NGC4713', 'NGC4808', 'NGC4396',
        #            'NGC4567', 'NGC4772', 'NGC4580', 'NGC4450', 'NGC4694', 'NGC4561']

        #galaxies = ['NGC4254', 'NGC4293', 'NGC4298', 'NGC4321', 'NGC4402',
        #            'NGC4424', 'NGC4457', 'NGC4535', 'NGC4536', 'NGC4548', 'NGC4569', 'NGC4579', 'NGC4654', 'NGC4689']

        #galaxies = ['NGC4064', 'NGC4222', 'NGC4294', 'NGC4330', 'NGC4388', 'NGC4394', 'NGC4402', 'NGC4405', 'NGC4419',
        #           ', 'NGC4533', 'NGC4567', 'NGC4606', 'NGC4607', 'NGC4772']  # These are the 7m only detections

        # galaxies = ['NGC4254']

        # galaxies = ['NGC4064']

        print('Input Galaxy ID(s): ' + ','.join(galaxies))
        print('')


        # if the Sun+18 masking method not used, ask user to input identifier to be used in filenames and headers
        if not sunflag:
            maskingmethod = input("Enter your masking method identifier (e.g., \"dame11\"):\n")
            writepath = writepath / maskingmethod

        # check to see if the paths exist, raising error if no readpath and making the write path if not.
        if readpath.exists():
            print("Cube directory found, reading from \n", readpath)
            print("")
        else:
            raise FileExistsError(readpath, " directory not found")

        if writepath.exists():
            print("Product directory exists, writing to \n", writepath)
            print("")
        else:
            writepath.mkdir(parents=True, exist_ok=False)
            print("Creating product directory. \n Writing to ", writepath)
            print("")
        
        # loop through the galaxies
        for i, galaxy in enumerate(galaxies):

            print("")
            print("------------------------------------------")
            print("##########################################")
            print("------------------------------------------")
            print("Beginning "+galaxy+"...")

            
            # create galaxy specific write path
            galaxywritepath = writepath / galaxy
            galaxywritepath.mkdir(parents=False, exist_ok=True)

            # galaxy file prefix for each product file
            galaxyfileprefix =  galaxy + '_7m+tp_co21_pbcorr_round_k_' + resolution_str + '_'

            # string save path + prefix
            savepath = str(galaxywritepath / galaxyfileprefix)

            if tosave:
                print("")
                print("Galaxy products will be saved in {}".format(galaxywritepath))
                print("with the file prefix: {}*".format(galaxyfileprefix))
                print("")

            # not sure why this is here?
            if galaxy == 'NGC4606' or galaxy == 'NGC4351':
                import matplotlib
                matplotlib.rcParams['text.usetex'] = False

            # assign the flat and pb corrected cube files and check they exist
            # Search for cube in data directory
            cubedir = readpath / str(galaxy +'/')
            if os.path.exists(cubedir):
                cubesearch = str(cubedir) + '/' + galaxy + '*_co21_*_round_k.fits' # use wildcard & glob to find cubes
                cubefiles = glob.glob(cubesearch)
            
            # select the 7m+TP cube if it's there AND the TP = True, else use the 7m data
            if TP and any("+tp" in s for s in cubefiles):
                cubefiles = [s for s in cubefiles if "+tp" in s]
            else:
                cubefiles = [s for s in cubefiles if "_7m_" in s]

            # pb corrected cube
            file_pbcorr = [s for s in cubefiles if 'pbcorr' in s][0]

            # flat, non-pb corrected cube (for rms calculation)
            file_uncorr = [s for s in cubefiles if 'flat' in s][0]
        
            print("")
            assert Path(file_pbcorr).exists(), "{} not found.".format(file_pbcorr)
            print("PB corrected cube: {}".format(file_pbcorr))
            assert Path(file_uncorr).exists(), "{} not found.".format(file_uncorr)
            print("Flat cube: {}".format(file_uncorr))
            print("")
            
            # clip cubes if True:
            if clip:

                print("")
                print("CLIPPING DATA CUBE:")
                C = ClipCube(galaxy
                            , file_pbcorr
                            , file_uncorr
                            , sun=sunflag
                            , savepath=savepath
                            , tosave=tosave)

                cube_corr, cube_uncorr = C.readfits()

                clipped_hdu, noisecube_hdu = C.do_clip(silent=False)
        


            #galaxy = 'NGC4222'
            #file_pbcorr = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/NGC4222_7m_co21_pbcorr_round_k.fits'
            #file_uncorr = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/NGC4222_7m_co21_flat_round_k.fits'
            #savepath = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/Products/Dame_method/'

            
            # try:
            #     data = fits.open(file_pbcorr)[0]
            #     pb = fits.open('/home/nikki/Documents/Data/VERTICO/ReducedData/' + galaxy + '/' +
            #                        galaxy + '_7m_co21_pb_rebin.fits')[0]

            #     if not pb.shape[1:3] == data.shape[1:3]:
            #         print('Wrong PB shape:')
            #         print(galaxy)
            #         continue
            # except:
            #     print('No PB info: ')
            #     print(galaxy)
            

            
            # if not pbcor:
            #     file_pbcorr = file_uncorr
            #     cube_corr = cube_uncorr.copy()
            #     if not os.path.exists(savepath_temp + 'PB_uncorrected/'):
            #         print("making product directory "+(savepath_temp + 'PB_uncorrected/')
            #         pathlib.Path(savepath_temp + 'PB_uncorrected/').mkdir(parents=True, exist_ok=True)
            #     if TP:
            #         savepath = savepath_temp + 'PB_uncorrected/' + galaxy + '_7m+tp_co21_flat_round_k.fits'
            #     else:
            #         savepath = savepath_temp + 'PB_uncorrected/' + galaxy + '_7m_co21_flat_round_k.fits'
            

            # # Have a first look at the cube to figure out some parameters
            # #plt.figure()
            # #plt.imshow(np.sum(cube_corr.data, axis=0), origin='lower')
            # #plt.figure()
            # #spec = np.sum(cube_corr.data, axis=(1, 2))
            # #plt.plot(spec)
            # #std = np.std(spec)
            # #x = np.arange(0, len(spec), 1)
            # #plt.plot(x, std * np.ones(len(x)))


            # Moment maps
            if 'mom0' in productList:
                print("")
                print("MOMENT 0:")

                for mom0units in ['K km/s', 'M_Sun/pc^2']:
                    CreateImages(galaxy, file_pbcorr, file_uncorr, alpha_co=alpha_co, savepath=savepath, refresh=refresh, overwrite=overwrite,
                            sun=sunflag, tosave=tosave).moment_zero(units=mom0units)

            if 'mom1' in productList:
                
                print("")
                print("MOMENT 1:")

                CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                        sun=sunflag, tosave=tosave).moment_1_2(moment=1)

            if 'mom2' in productList:
                
                print("")
                print("MOMENT 2:")

                CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sunflag, tosave=tosave).moment_1_2(moment=2)

            if any(c in productList for c in ('mom8', 'peakT')):

                print(productList)
                
                print("")
                print("MOMENT 8 (PEAK TEMPERATURE):")

                CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                        sun=sunflag, tosave=tosave).moment_zero(peak=True)
            
            # Uncertainty maps
            if 'mom0_noise' in productList:
                
                print("")
                print("MOMENT 0 UNCERTAINTY:")

                CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                        sun=sunflag, tosave=tosave).mom0_noise_maps()
            
            if 'mom1_noise' in productList:
                
                print("")
                print("MOMENT 1 & 2 UNCERTAINTY:")

                CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                        sun=sunflag, tosave=tosave).mom1_2_noise_maps()
            
            #if galaxy == 'NGC4561': continue

            # PVDs
            if 'pvd' in productList:

                print("")
                print("POSITION-VELOCITY DIAGRAM:")

                for axis in ['major', 'minor']:
                    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                            sun=sunflag, tosave=tosave).\
                                PVD(axis=axis, find_angle=False, check_slit=True)

            # Spectra
            if 'spectrum' in productList:

                print("")
                print("SPECTRUM:")

                for specUnits in ['velocity', 'vel_offset', 'frequency']:
                    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                            sun=sunflag, tosave=tosave).spectrum(x_axis=specUnits)

            # Radial profiles
            if 'radial_profile' in productList:

                print("")
                print("RADIAL PROFILE:")

                for xunits in ['kpc', 'arcsec']:
                    for yunits in ['K km/s', 'M_Sun pc^-2']:
                        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                                sun=sunflag, tosave=tosave).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                                    alpha_co=alpha_co, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                                    check_aperture=False)
         
    print("")
    print("FINISH: {:{tfmt}}".format(datetime.now(), dfmt='%Y%m%d', tfmt='%H:%M:%S'))
    elapsedTime = datetime.now() - starttime
    print("ELAPSED: ", str(elapsedTime))
    
    if logOutput == True:
        print("moved log file to ", writepath/logfile)
        os.rename("./"+logfile, writepath/logfile)

if __name__ == '__main__':
    
    main(logOutput=True)


