"""
Functions useful for data reduction
"""
import os
import sys
sys.path.append("/home/yamato/Application/analysis_scripts/")
import matplotlib.pyplot as plt
import numpy as np
import analysisUtils as au
import json
from casatasks import *
import casatools
import casatasks
import astropy.units as u
import astropy.constants as ac
import subprocess

tb = casatools.table()
msmd = casatools.msmetadata()
su = casatools.synthesisutils()

cms = ac.c.to(u.m/u.s).value

# from casatools import *

# def LSRKvel_to_chan(msfile, field, spw, restfreq, LSRKvelocity):
#     """
#     Identifies the channel(s) corresponding to input LSRK velocities.
#     Useful for choosing which channels to split out or flag if a line is expected to be present

#     Parameters
#     ==========
#     msfile: Name of measurement set (string)
#     spw: Spectral window number (int)
#     obsid: Observation ID corresponding to the selected spectral window
#     restfreq: Rest frequency in Hz (float)
#     LSRKvelocity: input velocity in LSRK frame in km/s (float or array of floats)

#     Returns
#     =======
#     Channel number most closely corresponding to input LSRK velocity
#     """
#     cc = 299792458. #speed of light in m/s

#     tb.open(msfile)
#     spw_col = tb.getcol('DATA_DESC_ID')
#     obs_col = tb.getcol('OBSERVATION_ID')

#     tb.close()
#     obsid = np.unique(obs_col[np.where(spw_col==spw)])
#     tb.open(msfile+'/SPECTRAL_WINDOW')
#     chanfreqs = tb.getcol('CHAN_FREQ', startrow = spw, nrow = 1)
#     tb.close()
#     tb.open(msfile+'/FIELD')
#     fieldnames = tb.getcol('NAME')
#     tb.close()
#     tb.open(msfile+'/OBSERVATION')
#     obstime = np.squeeze(tb.getcol('TIME_RANGE', startrow = obsid, nrow = 1))[0]
#     tb.close()
#     nchan = len(chanfreqs)
#     ms.open(msfile)
#     lsrkfreqs = ms.cvelfreqs(spwids = [spw], fieldids = np.where(fieldnames==field)[0][0], mode = 'channel', nchan = nchan, obstime = str(obstime)+'s', start = 0, outframe = 'LSRK')
#     chanvelocities = (restfreq-lsrkfreqs)/restfreq*cc/1.e3 #converted to LSRK velocities in km/s
#     ms.close()
#     if type(LSRKvelocity)==np.ndarray:
#         outchans = np.zeros_like(LSRKvelocity)
#         for i in range(len(LSRKvelocity)):
#             outchans[i] = np.argmin(np.abs(chanvelocities - LSRKvelocity[i]))
#         return outchans
#     else:
#         return np.argmin(np.abs(chanvelocities - LSRKvelocity))

# def get_flagchannels(ms_dict, output_prefix, exclude_spws = [], velocity_range = np.array([-20,20])):
#     """
#     Identify channels to flag based on provided velocity range of the line emission

#     Parameters
#     ==========
#     ms_dict: Dictionary of information about measurement set
#     output_prefix: Prefix for all output file names (string)
#     exclude_spws: spws to flag completely
#     velocity_range: Velocity range (in km/s) over which line emission has been identified, in the format np.array([min_velocity, max_velocity])

#     Returns
#     =======
#     String of channels to be flagged, in a format that can be passed to the spw parameter in CASA's flagdata task.
#     """

#     flagchannels_string = ''
#     for j,spw in enumerate(ms_dict['line_spws']):
#         if spw not in exclude_spws:
#             chans = LSRKvel_to_chan(ms_dict['vis'], ms_dict['field'], spw, ms_dict['line_freqs'][j] , velocity_range)
#             if j==0:
#                 flagchannels_string+='%d:%d~%d' % (spw, np.min([chans[0], chans[1]]), np.max([chans[0], chans[1]]))
#             else:
#                 flagchannels_string+=', %d:%d~%d' % (spw, np.min([chans[0], chans[1]]), np.max([chans[0], chans[1]]))
#         else:
#             print spw
#             if j==0: flagchannels_string+='%d' % (spw)
#             else: flagchannels_string+=', %d' % (spw)

#     print "# Flagchannels input string for %s: \'%s\'" % (ms_dict['name'], flagchannels_string)
#     return flagchannels_string

# def avg_cont(ms_dict, output_prefix, flagchannels = '', maxchanwidth = 125, datacolumn = 'data', contspws = None, width_array = None):
#     """
#     Produce spectrally averaged continuum measurement sets

#     Parameters
#     ==========
#     ms_dict: Dictionary of information about measurement set
#     output_prefix: Prefix for all output file names (string)
#     flagchannels: Argument to be passed for flagchannels parameter in flagdata task
#     maxchanwidth: Maximum width of channel (MHz). This is the value recommended by ALMA for Band 6 to avoid bandwidth smearing
#     datacolumn: Column to pull from for continuum averaging (usually will be 'data', but may sometimes be 'corrected' if there was flux rescaling applied)
#     contspws: Argument to be passed to CASA for the spw parameter in split. If not set, all SPWs will be selected by default. (string)
#     width_array: Argument to be passed to CASA for the width parameter in split. If not set, all SPWs will be selected by default. (array)
#     """
#     msfile = ms_dict['vis']
#     tb.open(msfile+'/SPECTRAL_WINDOW')
#     total_bw = tb.getcol('TOTAL_BANDWIDTH')
#     num_chan = tb.getcol('NUM_CHAN')
#     tb.close()

#     if width_array is None and contspws is None:
#         width_array = (num_chan/np.ceil(total_bw/(1.e6*maxchanwidth))).astype('int').tolist() #array of number of channels to average to form an output channel (to be passed to mstransform)
#         contspws = '%d~%d' % (0, len(total_bw)-1)#by default select all SPWs

#     elif (width_array is not None and contspws is None) or (width_array is None and contspws is not None):
#         print "If either contspws or width_array is set to a value, the other parameter has to be manually set as well"
#         return

#     if ms_dict['name']=='LB1':
#         timebin = '0s' #'6s'
#     else:
#         timebin = '0s' #default in CASA

#     #start of CASA commands

#     if len(flagchannels)==0:
#         outputvis = output_prefix+'_'+ms_dict['name']+'_initcont.ms*'
#         os.system('rm -rf '+outputvis)
#         split(vis=msfile,
#               field = ms_dict['field'],
#               spw = contspws,
#               outputvis = outputvis,
#               width = width_array,
#               timebin = timebin,
#               datacolumn=datacolumn,
#               intent = 'OBSERVE_TARGET#ON_SOURCE',
#               keepmms = True,
#               keepflags = False)
#     else:
#         if os.path.isdir(msfile+'.flagversions/flags.before_cont_flags'):
#             flagmanager(vis = msfile, mode = 'delete', versionname = 'before_cont_flags') # clear out old versions of the flags

#         print flagchannels
#         flagmanager(vis = msfile, mode = 'save', versionname = 'before_cont_flags', comment = 'Flag states before spectral lines are flagged') #save flag state before flagging spectral lines
#         flagdata(vis=msfile, mode='manual', spw=flagchannels, flagbackup=False, field = ms_dict['field']) #flag spectral lines

#         outputvis = output_prefix+'_'+ms_dict['name']+'_initcont.ms'
#         os.system('rm -rf '+outputvis)

#         split(vis=msfile,
#               field = ms_dict['field'],
#               spw = contspws,
#               outputvis = outputvis,
#               width = width_array,
#               timebin = timebin,
#               datacolumn=datacolumn,
#               intent = 'OBSERVE_TARGET#ON_SOURCE',
#               keepmms = True,
#               keepflags = False)

#         flagmanager(vis = msfile, mode = 'restore', versionname = 'before_cont_flags') #restore flagged spectral line channels

#     print "#Averaged continuum dataset saved to %s" % outputvis

def get_effective_bandwidth(vis, field, flagchannels):
    flagnchan = {}
    for spwchans in flagchannels.split(","):
        spw, chans = spwchans.split(":")
        flagnchan[spw] = 0
        for c in chans.split(";"):
            flagnchan[spw] += float(c.split("~")[1]) - float(c.split("~")[0]) + 1

    bandwidth = 0 
    msmd.open(vis)
    for spw in msmd.spwsforfield(field):
        bw = msmd.bandwidths(spw=spw)
        dchan = np.abs(msmd.chanwidths(spw=spw).mean())
        effbw = bw - flagnchan[str(spw)] * dchan
        bandwidth += effbw
    msmd.close()

    return bandwidth

def get_mean_frequency(vis, field):
    if type(vis) != list:
        vis = [vis]
    
    # get the mean frequencies for each MS and spw
    meanfreqs = []
    for ms in vis:
        msmd.open(ms)
        for i in msmd.spwsforfield(field):
            nu = msmd.meanfreq(spw=i, unit="Hz")
            meanfreqs.append(nu)
        msmd.close()

    # frequency averaged over all spws
    nu0 = np.mean(meanfreqs)

    return nu0


def get_imsize(vis, field, cellsize, D=12, pblimit=0.2):
    nu0 = get_mean_frequency(vis, field)

    # FoV based on antenna diameter (usually 12m) and frequency (in arcsec)
    FoV = 1.22 * cms / nu0 / D * 180 * 3600 / np.pi # note this is FWHM 
    FoV *= (-np.log(pblimit) / np.log(2)) ** 0.5
    nominal_imsize = FoV / float(cellsize.replace("arcsec", ""))

    # get optimum imsize for FFT
    imsize = su.getOptimumSize(int(nominal_imsize))

    return imsize

def tclean_continuum_wrapper(
    vis,
    imagename,
    # spw,
    imsize,
    cellsize,
    scales=None,
    phasecenter="",
    weighting="briggs",
    robust=0.5,
    niter=50000,
    mask="",
    noisethreshold=5.0,
    sidelobethreshold=3.0,
    minbeamfrac=0.3,
    lownoisethreshold=1.5,
    nterms=2,
    threshold="0.0Jy",
    nsigma=0.0,
    uvtaper=[],
    savemodel="none",
    parallel=True,
):

    ### mask setting
    if mask == "":
        usemask = "auto-multithresh"
    else:
        usemask = "user"

    for ext in [
        ".alpha*"
        ".image*",
        ".mask",
        ".model*",
        ".pb*",
        ".psf*",
        ".residual*",
        ".sumwt*",
        ".gridwt*",
    ]:
        os.system("rm -rf " + imagename + ext)
    tclean(
        vis=vis,
        imagename=imagename,
        # spw=spw,
        imsize=imsize,
        cell=cellsize,
        phasecenter=phasecenter,
        specmode="mfs",
        deconvolver="mtmfs",
        scales=scales,
        weighting=weighting,
        robust=robust,
        niter=niter,
        restoringbeam="common",
        mask=mask,
        usemask=usemask,
        noisethreshold=noisethreshold,
        sidelobethreshold=sidelobethreshold,
        lownoisethreshold=lownoisethreshold,
        minbeamfrac=minbeamfrac,
        nterms=nterms,
        threshold=threshold,
        nsigma=nsigma,
        parallel=parallel,
    )

    if savemodel == "modelcolumn":
        print("")
        print("Running tclean a second time to save the model...")
        tclean(
            vis=vis,
            imagename=imagename,
            # spw=spw,
            phasecenter=phasecenter,
            specmode="mfs",
            deconvolver="mtmfs",
            scales=scales,
            weighting=weighting,
            restoringbeam="common",
            robust=robust,
            imsize=imsize,
            cell=cellsize,
            niter=0,
            interactive=False,
            threshold=threshold,
            nsigma=nsigma,
            uvtaper=uvtaper,
            mask="",
            savemodel=savemodel,
            startmodel="",
            calcres=False,
            calcpsf=False,
            restart=True,
            parallel=False,
            nterms=nterms,
        )

def tclean_continuum_wrapper_old(
    vis,
    imagename,
    imsize,
    cellsize,
    scales='',
    smallscalebias=0.6,
    mask="",
    threshold="0.0mJy",
    interactive=False,
    robust=0.5,
    gain=0.3,
    niter=50000,
    cycleniter=300,
    uvtaper=[],
    savemodel="none",
    nsigma=0.0,
    parallel=True
):
    """
    Wrapper for tclean with keywords set to values desired for the Large Program imaging
    See the CASA 5.1.1 documentation for tclean to get the definitions of all the parameters
    """

    for ext in [
        ".image*",
        ".mask",
        ".model*",
        ".pb*",
        ".psf*",
        ".residual*",
        ".sumwt*",
        ".gridwt*",
    ]:
        os.system("rm -rf " + imagename + ext)
    tclean(
        vis=vis,
        imagename=imagename,
        specmode="mfs",
        deconvolver="mtmfs",
        scales=scales,
        weighting="briggs",
        robust=robust,
        gain=gain,
        imsize=imsize,
        cell=cellsize,
        smallscalebias=smallscalebias,  # set to CASA's default of 0.6 unless manually changed
        niter=niter,  # we want to end on the threshold
        interactive=interactive,
        threshold=threshold,
        nsigma=nsigma,
        cycleniter=cycleniter,
        cyclefactor=1,
        uvtaper=uvtaper,
        mask=mask,
        savemodel="none",
        parallel=parallel,
        nterms=1,
    )

    # this step is a workaround a bug in tclean that doesn't always save the model during multiscale clean. See the "Known Issues" section for CASA 5.1.1 on NRAO's website
    # still persists in casa 5.6 but now works in mpicasa (did not work in casa 5.4 ...)

    if savemodel == "modelcolumn":
        print("")
        print("Running tclean a second time to save the model...")
        tclean(
            vis=vis,
            imagename=imagename,
            specmode="mfs",
            deconvolver="mtmfs",
            scales=scales,
            weighting="briggs",
            robust=robust,
            gain=gain,
            imsize=imsize,
            cell=cellsize,
            smallscalebias=smallscalebias,  # set to CASA's default of 0.6 unless manually changed
            niter=0,
            interactive=False,
            threshold=threshold,
            nsigma=nsigma,
            cycleniter=cycleniter,
            cyclefactor=1,
            uvtaper=uvtaper,
            mask="",
            savemodel=savemodel,
            startmodel="",
            calcres=False,
            calcpsf=False,
            restart=True,
            parallel=False,
            nterms=1,
        )

def get_spectral_params(vrange, channelwidth):
    if vrange is not None:
        v1, v2 = vrange
        start = "{}km/s".format(v1)
    else:
        start = ""
    if channelwidth is not None:
        width = "{}km/s".format(channelwidth)
        if vrange is not None:
            nchan = int((v2 - v1) / channelwidth)
        else:
            nchan = -1
    else:
        width = ""
        nchan = -1

    return start, width, nchan

def tclean_spectral_line_wrapper(
    vis,
    imagename,
    spw,
    imsize,
    cellsize,
    scales='',
    phasecenter="",
    weighting="briggsbwtaper",
    perchanweightdensity=True,
    robust=0.5,
    niter=50000,
    mask="",
    noisethreshold=5.0,
    sidelobethreshold=3.0,
    minbeamfrac=0.3,
    lownoisethreshold=1.5,
    threshold="0.0Jy",
    nsigma=0.0,
    restfreq='',
    vrange=None,
    channelwidth=None,
    start='',
    width='',
    nchan=-1,
    parallel=True,
):
    ### mask setting
    if mask == "":
        usemask = "auto-multithresh"
    else:
        usemask = "user"

    # get params for spectral gridding
    if (vrange is not None) and (channelwidth is not None):
        start, width, nchan = get_spectral_params(vrange, channelwidth)
        print(start, width, nchan)

    for ext in [
        ".image",
        ".mask",
        ".model",
        ".pb",
        ".psf",
        ".residual",
        ".sumwt",
        ".gridwt",
    ]:
        os.system("rm -r " + imagename + ext)

    # first run tclean with niter=0, calcpsf=True, and calcres=False
    tclean(
        vis=vis,
        imagename=imagename,
        spw=spw,
        imsize=imsize,
        cell=cellsize,
        phasecenter=phasecenter,
        specmode="cube",
        start=start,
        width=width,
        nchan=nchan,
        restfreq=restfreq,
        outframe="LSRK",
        veltype="radio",
        deconvolver="multiscale",
        scales=scales,
        weighting=weighting,
        perchanweightdensity=perchanweightdensity,
        robust=robust,
        niter=0, # no iteration
        calcpsf=True, # only calculate psf
        calcres=False, # no residual calculation
        restoringbeam="common",
        mask=mask,
        usemask=usemask,
        noisethreshold=noisethreshold,
        sidelobethreshold=sidelobethreshold,
        lownoisethreshold=lownoisethreshold,
        minbeamfrac=minbeamfrac,
        threshold=threshold,
        nsigma=nsigma,
        parallel=parallel,
    )

    # check the beam shape through the channels and get rid of outlier
    ia.open(imagename+".psf")
    perchanbeam = ia.restoringbeam()["beams"]
    channel = sorted([int(i.replace("*", "")) for i in perchanbeam.keys()])
    perchanbeamarea = np.array([perchanbeam["*"+str(chan)]['*0']['major']['value']*perchanbeam["*"+str(chan)]['*0']['minor']['value'] for chan in channel]) # in arcsec2
    badchan = np.full(perchanbeamarea.shape, False)
    for _ in range(3): # heuristics
        median, stddev = np.nanmedian(perchanbeamarea[~badchan]), np.nanstd(perchanbeamarea[~badchan])
        badchan = np.abs(perchanbeamarea - median) >= 3*stddev # detecting outlier
    
    if np.any(badchan):
        print("Found bad beams that significantly deviated. Going to correct those at the spw edge...")
        # serach bad channels at the edge from channel=0
        for bad in badchan:
            if not bad:
                break
            else:
                if isinstance(start, int):
                    start += 1
                elif ("Hz" in start) or ("m/s" in start):
                    if width != '':
                        start_tmp = u.Quantity(start) + u.Quantity(width)
                        start = '{}{}'.format(start_tmp.value, start_tmp.unit.to_string())
                else:
                    start = 1
                nchan = nchan - 1 if nchan != -1 else len(badchan) - 1
        
        for bad in badchan[::-1]:
            if not bad:
                break
            else:
                nchan = nchan - 1 if nchan != -1 else len(badchan) - 1
                
        # badchan = np.where(bad)[0]
        # continuous = True if len(badchan) == 1 else np.all(np.diff(badchan) == 1)
        # if (0 in badchan) and continuous:
        #     print("Bad beams at the spw edge. Going to remove them...")
        #     start = len(badchan)
        # elif (len(bad)-1 in badchan) and continuous:
        #     print("Bad beams at the spw edge. Going to remove them...")
        #     nchan = len(bad) - len(badchan)
        # else:
        #     perchanbmaj = np.array([perchanbeam[chan]['*0']['major']['value'] for chan in perchanbeam.keys()])
        #     perchanbmin = np.array([perchanbeam[chan]['*0']['minor']['value'] for chan in perchanbeam.keys()])
        #     perchanbpa = np.array([perchanbeam[chan]['*0']['positionangle']['value'] for chan in perchanbeam.keys()])
        #     goodchan = np.where(~bad)
        #     bmaj, bmin, bpa = np.nanmean(perchanbmaj[goodchan]), np.nanmean(perchanbmin[goodchan]), np.nanmean(perchanbpa[goodchan])
        #     for chan in badchan:
        #         ia.setrestoringbeam(major=f"{bmaj}arcsec", minor=f"{bmin}arcsec", pa=f"{bpa}deg", channel=chan)
    ia.close()
    
    for ext in [".gridwt", ".pb", ".psf", ".sumwt"]:
        os.system("rm -r " + imagename + ext)

    # run the tclean a second time
    tclean(
        vis=vis,
        imagename=imagename,
        spw=spw,
        imsize=imsize,
        cell=cellsize,
        phasecenter=phasecenter,
        specmode="cube",
        start=start,
        width=width,
        nchan=nchan,
        restfreq=restfreq,
        outframe="LSRK",
        veltype="radio",
        deconvolver="multiscale",
        scales=scales,
        weighting=weighting,
        perchanweightdensity=perchanweightdensity,
        robust=robust,
        niter=niter,
        restoringbeam="common",
        # calcpsf=False, # assume that .psf exists -> NO! nchan is changed thus need to calculate again
        mask=mask,
        usemask=usemask,
        noisethreshold=noisethreshold,
        sidelobethreshold=sidelobethreshold,
        lownoisethreshold=lownoisethreshold,
        minbeamfrac=minbeamfrac,
        threshold=threshold,
        nsigma=nsigma,
        parallel=parallel,
    )


def calc_sensitivity(
    msfile,
    specmode="mfs",
    spw=[],
    chan=0,
    cellsize="0.03arcsec",
    imsize=700,
    robust=0.5,
    uvtaper="",
):
    """_summary_

    Parameters
    ----------
    msfile : str or list
        measurement set file name or list of them
    specmode : str, optional
        'specmode' parameter used in tclean, by default 'mfs'
    spw : list, optional
        spw selection, by default ''
    chan : int, optional
        channel selection, provide only one representative channel. Used for sensitivity estimates of line cube. by default 0
    cellsize : str, optional
        'cell' parameter used in tclean, by default '0.03arcsec'
    imsize : int, optional
        'imsize' parameter used in tclean, by default 700
    robust : float, optional
        'robust' parameter used in tclean, by default 0.5
    uvtaper : str, optional
        'uvtaper' parameter used in tclean, by default ''. Note that only in the unit of 'arcsec' is acceptable currently

    Returns
    -------
    float
        estimated sensitivity averaged over provided measurement set(s)
    """

    msfile = [msfile] if (type(msfile) != list) else msfile
    sensitivities = []
    print("Theoretical estimates of image sensitivity")
    for vis in msfile:
        print("# " + vis)
        if specmode == "mfs":
            spwstring = ""
            print("# specmode = 'mfs'")
            print("# spw = '{:s}'".format(spwstring))
        if specmode == "cube":
            print("# specmode = 'cube'")
            spwstring = (
                ",".join(str(i) for i in spw) + ":" + str(chan) + "~" + str(chan)
            )
            print("# spw = '{:s}'".format(spwstring))
        im.open(vis)
        im.selectvis(field="", spw=spwstring)
        im.defineimage(
            mode=specmode,
            stokes="I",
            spw=spw,
            cellx=cellsize,
            celly=cellsize,
            nx=imsize,
            ny=imsize,
        )
        im.weight(type="briggs", robust=robust)
        if uvtaper != "":
            if type(uvtaper) == list:
                bmaj, bmin, bpa = uvtaper
            elif type(uvtaper) == str:
                bmaj = uvtaper
                bmin = bmaj
                bpa = "0.0deg"
            print(f"# uvtaper = [{bmaj}, {bmin}, {bpa}]")
            im.filter(type="gaussian", bmaj=bmaj, bmin=bmin, bpa=bpa)
        sens = im.apparentsens()
        if sens[0]:
            print(sens)
            print(f"# Briggs sensitivity with robust = {robust}: {sens[1]} Jy")
            print(f"# Relative to natural weighting: {sens[2]}")
            sensitivities.append(sens[1])
    avesens = np.mean(sensitivities) / len(sensitivities) ** 0.5
    print(
        f"### Averaged sensitivity over {len(sensitivities)} measurement sets: {avesens} Jy"
    )
    return avesens


def image_cont_each_obs(
    msfile,
    scales,
    smallscalebias=0.6,
    mask="",
    threshold="0.0mJy",
    nsigma=0.0,
    imsize=None,
    cellsize=None,
    interactive=False,
    robust=0.5,
    gain=0.3,
    niter=50000,
    cycleniter=300,
):
    """
    Wrapper for tclean that will loop through all the observations in a measurement set and image them individual

    Parameters
    ==========
    ms_dict: Dictionary of information about measurement set
    prefix: Prefix for all output file names (string)

    See the CASA 5.1.1 documentation for tclean to get the definitions of all other parameters
    """

    imsize = 720 if imsize is None else imsize
    cellsize = "0.03arcsec" if cellsize is None else cellsize

    # start of CASA commands
    # for i in range(num_observations):
    imagename = msfile.replace(".ms", "")
    for ext in [
        ".image*",
        ".mask",
        ".model*",
        ".pb*",
        ".psf*",
        ".residual*",
        ".sumwt*",
    ]:
        os.system("rm -rf " + imagename + ext)
    tclean(
        vis=msfile,
        imagename=imagename,
        specmode="mfs",
        deconvolver="mtmfs",
        scales=scales,
        weighting="briggs",
        robust=robust,
        gain=gain,
        imsize=imsize,
        cell=cellsize,
        smallscalebias=smallscalebias,  # set to CASA's default of 0.6 unless manually changed
        niter=niter,  # we want to end on the threshold
        interactive=interactive,
        threshold=threshold,
        nsigma=nsigma,
        cycleniter=cycleniter,
        cyclefactor=1,
        mask=mask,
        nterms=1,
    )
    print("Image saved in %s.image.tt0" % (imagename))


def fit_gaussian(imagename, region, dooff=False):
    """
    Wrapper for imfit in CASA to fit a single Gaussian component to a selected region of the image
    Parameters
    ==========
    imagename: Name of CASA image (ending in .image) (string)
    region: CASA region format, e.g., 'circle[[200pix, 200pix], 3arcsec]' (string)
    dooff: boolean option to allow for fitting a zero-level offset
    """
    imfitdict = imfit(imagename=imagename, region=region, dooff=dooff)
    # Check if the source was resolved
    was_resolved = not imfitdict["deconvolved"]["component0"]["ispoint"]
    # Get the coordinate system
    coordsystem = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["refer"]
    # Get the parameters
    headerlist = imhead(imagename)
    phasecenter_ra, phasecenter_dec = headerlist["refval"][:2]
    peak_ra = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["m0"][
        "value"
    ]
    peak_dec = imfitdict["deconvolved"]["component0"]["shape"]["direction"]["m1"][
        "value"
    ]
    xcen, ycen = headerlist["refpix"][:2]
    deltax, deltay = headerlist["incr"][:2]
    peak_x = xcen + np.unwrap(np.array([0, peak_ra - phasecenter_ra]))[
        1
    ] / deltax * np.cos(phasecenter_dec)
    peak_y = ycen + (peak_dec - phasecenter_dec) / deltay
    # Print
    if coordsystem == "J2000":
        coords = au.rad2radec(imfitdict=imfitdict, hmsdms=True, delimiter=" ")
        print("#Peak of Gaussian component identified with imfit: J2000 %s" % coords)
    elif coordsystem == "ICRS":
        coords = au.rad2radec(imfitdict=imfitdict, hmsdms=True, delimiter=" ")
        print("#Peak of Gaussian component identified with imfit: ICRS %s" % coords)
        J2000coords = au.ICRSToJ2000(au.rad2radec(imfitdict=imfitdict, delimiter=" "))
        print("#Peak in J2000 coordinates: %s" % J2000coords)
    else:
        print(
            "#If the coordinates aren't in ICRS or J2000, then something weird is going on"
        )
    # If the object was resolved, print the inclination, PA, major and minor axis
    if was_resolved:
        PA = imfitdict["deconvolved"]["component0"]["shape"]["positionangle"]["value"]
        majoraxis = imfitdict["deconvolved"]["component0"]["shape"]["majoraxis"][
            "value"
        ]
        minoraxis = imfitdict["deconvolved"]["component0"]["shape"]["minoraxis"][
            "value"
        ]

        print("#PA of Gaussian component: %.2f deg" % PA)
        print(
            "#Inclination of Gaussian component: %.2f deg"
            % (np.arccos(minoraxis / majoraxis) * 180 / np.pi,)
        )
    print("#Pixel coordinates of peak: x = %.3f y = %.3f" % (peak_x, peak_y))
    if coordsystem == "J2000":
        return "J2000 " + coords
    elif coordsystem == "ICRS":
        return "J2000 " + au.convertColonDelimitersToHMSDMS(J2000coords).replace(
            ",", ""
        )


# def split_all_obs(msfile, nametemplate):
#     """

#     Split out individual observations in a measurement set

#     Parameters
#     ==========
#     msfile: Name of measurement set, ending in '.ms' (string)
#     nametemplate: Template name of output measurement sets for individual observations (string)
#     """
#     tb.open(msfile)
#     spw_col = tb.getcol('DATA_DESC_ID')
#     obs_col = tb.getcol('OBSERVATION_ID')
#     field_col = tb.getcol('FIELD_ID')
#     tb.close()

#     obs_ids = np.unique(obs_col)


#     #yes, it would be more logical to split out by observation id, but splitting out by observation id in practice leads to some issues with the metadata
#     for i in obs_ids:
#         spws = np.unique(spw_col[np.where(obs_col==i)])
#         fields = np.unique(field_col[np.where(obs_col==i)]) #sometimes the MS secretly has multiple field IDs lurking even if listobs only shows one field
#         if len(spws)==1:
#             spw = str(spws[0])
#         else:
#             spw = "%d~%d" % (spws[0], spws[-1])

#         if len(fields)==1:
#             field = str(fields[0])
#         else:
#             field = "%d~%d" % (fields[0], fields[-1])
#         #start of CASA commands
#         outputvis = nametemplate+'%d.ms' % i
#         os.system('rm -rf '+outputvis)
#         print "#Saving observation %d of %s to %s" % (i, msfile, outputvis)
#         split(vis=msfile,
#               spw = spw,
#               field = field,
#               outputvis = outputvis,
#               datacolumn='data')


def export_MS(msfile):
    """
    Spectrally averages visibilities to a single channel per SPW and exports to .npz file

    msfile: Name of CASA measurement set, ending in '.ms' (string)
    """
    filename = msfile
    # strip off the '.ms'
    MS_filename = filename.replace(".ms", "")

    # get information about spectral windows

    tb.open(MS_filename + ".ms/SPECTRAL_WINDOW")
    num_chan = tb.getcol("NUM_CHAN").tolist()
    tb.close()

    # spectral averaging (1 channel per SPW)

    os.system("rm -rf %s" % MS_filename + "_spavg.ms")
    split(
        vis=MS_filename + ".ms",
        width=np.max(num_chan),
        datacolumn="data",
        outputvis=MS_filename + "_spavg.ms",
    )

    # get the data tables
    tb.open(MS_filename + "_spavg.ms")
    data = np.squeeze(tb.getcol("DATA"))
    flag = np.squeeze(tb.getcol("FLAG"))
    uvw = tb.getcol("UVW")
    weight = tb.getcol("WEIGHT")
    spwid = tb.getcol("DATA_DESC_ID")
    tb.close()

    # get frequency information
    tb.open(MS_filename + "_spavg.ms/SPECTRAL_WINDOW")
    freqlist = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.close()

    # get rid of any flagged columns
    good = np.squeeze(np.any(flag, axis=0) == False)
    data = data[:, good]
    weight = weight[:, good]
    uvw = uvw[:, good]
    spwid = spwid[good]

    # Fedeles'data include only 1 spw (SB3&4)
    if freqlist.size == 1:
        freqlist = np.asarray([freqlist])

    # compute spatial frequencies in lambda units
    get_freq = lambda ispw: freqlist[ispw]
    freqs = get_freq(spwid)  # get spectral frequency corresponding to each datapoint
    u = uvw[0, :] * freqs / 2.9979e8
    v = uvw[1, :] * freqs / 2.9979e8

    # average the polarizations
    Re = np.sum(data.real * weight, axis=0) / np.sum(weight, axis=0)
    Im = np.sum(data.imag * weight, axis=0) / np.sum(weight, axis=0)
    Vis = Re + 1j * Im
    Wgt = np.sum(weight, axis=0)

    # output to npz file and delete intermediate measurement set
    os.system("rm -rf %s" % MS_filename + "_spavg.ms")
    os.system("rm -rf " + MS_filename + ".vis.npz")
    np.savez(MS_filename + ".vis", u=u, v=v, Vis=Vis, Wgt=Wgt)
    print("#Measurement set exported to %s" % (MS_filename + ".vis.npz",))
    return MS_filename + ".vis.npz"


def deproject_vis(
    data, bins=np.array([0.0]), incl=0.0, PA=0.0, offx=0.0, offy=0.0, errtype="mean"
):
    """
    Deprojects and azimuthally averages visibilities

    Parameters
    ==========
    data: Length-4 tuple of u,v, visibilities, and weight arrays
    bins: 1-D array of uv distance bins (kilolambda)
    incl: Inclination of disk (degrees)
    PA: Position angle of disk (degrees)
    offx: Horizontal offset of disk center from phase center (arcseconds)
    offy: Vertical offset of disk center from phase center (arcseconds)

    Returns
    =======
    uv distance bins (1D array), visibilities (1D array), errors on averaged visibilities (1D array)
    """

    # - read in, parse data
    u, v, vis, wgt = data
    # - convert keywords into relevant units
    inclr = np.radians(incl)
    PAr = 0.5 * np.pi - np.radians(PA)
    offx *= -np.pi / (180.0 * 3600.0)
    offy *= -np.pi / (180.0 * 3600.0)

    # - change to a deprojected, rotated coordinate system
    uprime = u * np.cos(PAr) + v * np.sin(PAr)
    vprime = (-u * np.sin(PAr) + v * np.cos(PAr)) * np.cos(inclr)
    rhop = np.sqrt(uprime**2 + vprime**2)

    # - phase shifts to account for offsets
    shifts = np.exp(-2.0 * np.pi * 1.0j * (u * -offx + v * -offy))
    visp = vis * shifts
    realp = visp.real
    imagp = visp.imag

    # - if requested, return a binned (averaged) representation
    if bins.size > 1.0:
        avbins = 1e3 * bins  # scale to lambda units (input in klambda)
        bwid = 0.5 * (avbins[1] - avbins[0])
        bvis = np.zeros_like(avbins, dtype="complex")
        berr = np.zeros_like(avbins, dtype="complex")
        for ib in np.arange(len(avbins)):
            inb = np.where((rhop >= avbins[ib] - bwid) & (rhop < avbins[ib] + bwid))
            if len(inb[0]) >= 5:
                bRe, eRemu = np.average(realp[inb], weights=wgt[inb], returned=True)
                eRese = np.std(realp[inb])
                bIm, eImmu = np.average(imagp[inb], weights=wgt[inb], returned=True)
                eImse = np.std(imagp[inb])
                bvis[ib] = bRe + 1j * bIm
                if errtype == "scat":
                    berr[ib] = eRese + 1j * eImse
                else:
                    berr[ib] = 1.0 / np.sqrt(eRemu) + 1j / np.sqrt(eImmu)
            else:
                bvis[ib] = 0 + 1j * 0
                berr[ib] = 0 + 1j * 0
            parser = np.where(berr.real != 0)
            output = avbins[parser], bvis[parser], berr[parser]
        return output

    # - if not, returned the unbinned representation
    output = rhop, realp + 1j * imagp, 1.0 / np.sqrt(wgt)

    return output


def plot_deprojected(
    filelist,
    incl=0,
    PA=0,
    offx=0,
    offy=0,
    fluxscale=None,
    uvbins=None,
    show_err=True,
    outfile="vis_profile.png",
):
    """
    Plots real and imaginary deprojected visibilities from a list of .npz files

    Parameters
    ==========
    filelist: List of names of .npz files storing visibility data
    incl: Inclination of disk (degrees)
    PA: Position angle of disk (degrees)
    offx: Horizontal offset of disk center from phase center (arcseconds)
    offy: Vertical offset of disk center from phase center (arcseconds)
    fluxscale: List of scaling factors to multiply the visibility values by before plotting. Default value is set to all ones.
    uvbins: Array of bins at which to plot the visibility values, in lambda. By default, the range plotted will be from 10 to 1000 kilolambda
    show_err: If True, plot error bars.
    outfile: filename in which the figure is saved
    """
    if fluxscale is None:
        fluxscale = np.ones(len(filelist))
    assert len(filelist) == len(fluxscale)

    if uvbins is None:
        uvbins = 10.0 + 10.0 * np.arange(100)

    minvis = np.zeros(len(filelist))
    maxvis = np.zeros(len(filelist))
    maxrho = np.zeros(len(filelist))
    fig, ax = plt.subplots(2, 1, sharex=True)
    for i, filename in enumerate(filelist):

        # read in the data
        inpf = np.load(filename)
        u = inpf["u"]
        v = inpf["v"]
        vis = fluxscale[i] * inpf["Vis"]
        wgt = inpf["Wgt"]

        # deproject the visibilities and do the annular averaging
        vp = deproject_vis(
            [u, v, vis, wgt], bins=uvbins, incl=incl, PA=PA, offx=offx, offy=offy
        )
        vp_rho, vp_vis, vp_sig = vp

        # calculate min, max of deprojected, averaged reals (for visualization)
        minvis[i] = np.min(vp_vis.real)
        maxvis[i] = np.max(vp_vis.real)
        maxrho[i] = np.max(vp_rho)

        # plot the profile
        if show_err:
            ax[0].errorbar(
                1e-3 * vp_rho,
                vp_vis.real,
                yerr=vp_sig.real,
                label=os.path.basename(filename),
                fmt=".",
            )
            ax[1].errorbar(
                1e-3 * vp_rho,
                vp_vis.imag,
                yerr=vp_sig.imag,
                label=os.path.basename(filename),
                fmt=".",
            )
        else:
            ax[0].plot(1e-3 * vp_rho, vp_vis.real, "o", markersize=2.8, label=filename)
            ax[1].plot(1e-3 * vp_rho, vp_vis.imag, "o", markersize=2.8, label=filename)

    allmaxvis = np.max(maxvis)
    allminvis = np.min(minvis)
    allmaxrho = np.max(maxrho)
    if (allminvis < 0) or (allminvis - 0.1 * allmaxvis < 0):
        ax[0].axis(
            [0, 1e-3 * allmaxrho * 1.1, allminvis - 0.1 * allmaxvis, 1.1 * allmaxvis]
        )
        ax[1].axis(
            [0, 1e-3 * allmaxrho * 1.1, allminvis - 0.1 * allmaxvis, 1.1 * allmaxvis]
        )
    else:
        ax[0].axis([0, 1e-3 * allmaxrho * 1.1, 0.0, 1.1 * allmaxvis])
        ax[1].axis([0, 1e-3 * allmaxrho * 1.1, 0.0, 1.1 * allmaxvis])

    ax[0].plot([0, 1e-3 * allmaxrho * 1.1], [0, 0], "--k")
    ax[1].plot([0, 1e-3 * allmaxrho * 1.1], [0, 0], "--k")
    plt.xlabel("deprojected baseline length [kilo$\lambda$]")
    ax[0].set_ylabel("average real [Jy]")
    ax[1].set_ylabel("average imag [Jy]")
    ax[0].legend()
    # plt.show(block = False)
    fig.savefig(outfile, bbox_inches="tight", pad_inche=0.01, dpi=500)


def estimate_flux_scale(
    reference,
    comparison,
    incl=0,
    PA=0,
    uvbins=None,
    offx=0,
    offy=0,
    outfile="ratio_profile.png",
):
    """
    Calculates the weighted average of the flux ratio between two observations of a source
    The minimum baseline compared is the longer of the minimum baselines in the individual datasets
    The longest baseline compared is either the shorter of the longest baselines in the individual datasets, or 800 kilolambda

    Parameters
    ==========
    reference: Name of .npz file holding the reference dataset (with the "correct" flux")
    comparison: Name of .npz file holding the comparison dataset (with the flux ratio being checked)
    filelist: List of names of .npz files storing visibility data
    incl: Inclination of disk (degrees)
    PA: Position angle of disk (degrees)
    offx: Horizontal offset of disk center from phase center (arcseconds)
    offy: Vertical offset of disk center from phase center (arcseconds)
    uvbins: Array of bins at which to compare the visibility values, in lambda.
            By default, the minimum baseline compared is the longer of the minimum baselines in the individual datasets.
            The longest baseline compared is either the shorter of the longest baselines in the individual datasets, or 800 kilolambda, whichever comes first.
    outfile: filename in which the ratio profile figure is saved
    """

    inpf = np.load(reference)
    u_ref = inpf["u"]
    v_ref = inpf["v"]
    vis_ref = inpf["Vis"]
    wgt_ref = inpf["Wgt"]

    inpf = np.load(comparison)
    u_comp = inpf["u"]
    v_comp = inpf["v"]
    vis_comp = inpf["Vis"]
    wgt_comp = inpf["Wgt"]

    uvdist_ref = np.sqrt(u_ref**2 + v_ref**2)
    uvdist_comp = np.sqrt(u_comp**2 + v_comp**2)

    mindist = np.max(np.array([np.min(uvdist_ref), np.min(uvdist_comp)]))
    maxdist = np.min(
        np.array([np.max(uvdist_ref), np.max(uvdist_ref), 8e5])
    )  # the maximum baseline we want to compare is the longest shared baseline or 800 kilolambda, whichever comes first (we don't want to go out to a baseline that's too long because phase decorrelation becomes a bigger issue at longer baselines.

    if uvbins is None:
        uvbins = mindist / 1.0e3 + 10.0 * np.arange(
            np.floor((maxdist - mindist) / 1.0e4)
        )

    # deproject the visibilities and do the annular averaging
    vp = deproject_vis(
        [u_ref, v_ref, vis_ref, wgt_ref],
        bins=uvbins,
        incl=incl,
        PA=PA,
        offx=offx,
        offy=offy,
    )
    ref_rho, ref_vis, ref_sig = vp

    # deproject the visibilities and do the annular averaging
    vp = deproject_vis(
        [u_comp, v_comp, vis_comp, wgt_comp],
        bins=uvbins,
        incl=incl,
        PA=PA,
        offx=offx,
        offy=offy,
    )
    comp_rho, comp_vis, comp_sig = vp

    maxlen = np.min(np.array([len(comp_rho), len(ref_rho)]))

    rho_intersection = np.intersect1d(
        ref_rho, comp_rho
    )  # we only want to compare overlapping baseline intervals

    comp_sig_intersection = comp_sig[
        np.where(np.in1d(comp_rho, rho_intersection))
    ].real  # they're the same for the real and imaginary components
    comp_vis_intersection = comp_vis[np.where(np.in1d(comp_rho, rho_intersection))]
    ref_sig_intersection = ref_sig[np.where(np.in1d(ref_rho, rho_intersection))].real
    ref_vis_intersection = ref_vis[np.where(np.in1d(ref_rho, rho_intersection))]

    ratio = np.abs(comp_vis_intersection) / np.abs(ref_vis_intersection)
    err = ratio * np.sqrt(
        (comp_sig_intersection / np.abs(comp_vis_intersection)) ** 2
        + (ref_sig_intersection / np.abs(ref_vis_intersection)) ** 2
    )

    w = 1 / err**2
    ratio_avg = np.sum(w * ratio) / np.sum(w)
    print(
        "#The ratio of the fluxes of %s to %s is %.5f"
        % (comparison, reference, ratio_avg)
    )
    print(
        "#The scaling factor for gencal is %.3f for your comparison measurement"
        % (np.sqrt(ratio_avg))
    )
    print(
        "#The error on the weighted mean ratio is %.3e, although it's likely that the weights in the measurement sets are too off by some constant factor"
        % (1 / np.sqrt(np.sum(w)),)
    )
    plt.figure()
    plt.errorbar(
        1e-3 * rho_intersection, ratio, yerr=err, fmt=".", label="Binned ratios"
    )
    plt.plot(
        1e-3 * rho_intersection,
        np.ones_like(ratio) * ratio_avg,
        label="weighted average",
    )
    plt.ylabel("Visibility amplitude ratios")
    plt.xlabel("UV distance (kilolambda)")
    plt.legend()
    # plt.show(block = False)
    plt.savefig(outfile, bbox_inches="tight", pad_inches=0.01, dpi=500)


def rescale_flux(vis, gencalparameter):
    """
    Rescale visibility fluxes using gencal, then split into a new measurement set

    Parameters:
    vis: Measurement set name, ending in ms (string)
    gencalparameter: List of flux rescaling parameters to be passed to 'parameter' for gencal task
    """
    caltable = "scale_" + vis.replace(".ms", ".gencal")
    os.system("rm -rf " + caltable)
    gencal(vis=vis, caltable=caltable, caltype="amp", parameter=gencalparameter)
    applycal(vis=vis, gaintable=caltable, calwt=True, flagbackup=True)
    vis_rescaled = vis.replace(".ms", "_rescaled.ms")
    print("#Splitting out rescaled values into new MS: %s" % (vis_rescaled,))
    os.system("rm -rf " + vis_rescaled + "*")
    split(vis=vis, outputvis=vis_rescaled, datacolumn="corrected")


def estimate_SNR(imagename, disk_mask, noise_mask):
    """
    Estimate peak SNR of source

    Parameters:
    imagename: Image name ending in '.image' (string)
    disk_mask: , in the CASa region format, e.g.
    noise_mask: Annulus to measure image rms, in the CASA region format, e.g. 'annulus[[500pix, 500pix],["1arcsec", "2arcsec"]]' (string)
    """
    headerlist = imhead(imagename, mode="list")
    beammajor = headerlist["beammajor"]["value"]
    beamminor = headerlist["beamminor"]["value"]
    beampa = headerlist["beampa"]["value"]
    disk_stats = imstat(imagename=imagename, region=disk_mask)
    disk_flux = disk_stats["flux"][0]
    peak_intensity = disk_stats["max"][0]
    rms = imstat(imagename=imagename, region=noise_mask)["rms"][0]
    SNR = peak_intensity / rms

    print("#%s" % imagename)
    print("#Beam %.3f arcsec x %.3f arcsec (%.2f deg)" % (beammajor, beamminor, beampa))
    print("#Flux inside disk mask: %.2f mJy" % (disk_flux * 1000,))
    print("#Peak intensity of source: %.2f mJy/beam" % (peak_intensity * 1000,))
    print("#rms: %.2e mJy/beam" % (rms * 1000,))
    print("#Peak SNR: %.2f" % (SNR,))

    return SNR, rms

### DO NOT WORK WELL IN CASA 6.2.1 ###
# def get_reference_antennas(vis, pipelinecasalog):
#     vis = os.path.basename(vis)
#     print(vis, pipelinecasalog)
#     _, totalscore, visname = au.readRefantScoresFromCasalog(pipelinecasalog, vis=vis)
#     assert vis == visname
#     name = list(totalscore.keys())
#     score = [totalscore[key] for key in totalscore.keys()]
#     return ",".join(np.array(name)[np.argsort(score)[::-1]])

def rank_reference_antennas(vis):
    # Get the antenna names and offsets.

    msmd.open(vis)
    names = msmd.antennanames()
    offset = [msmd.antennaoffset(name) for name in names]
    msmd.close()

    # Calculate the mean longitude and latitude.

    mean_longitude = np.mean(
        [offset[i]["longitude offset"]["value"] for i in range(len(names))]
    )
    mean_latitude = np.mean(
        [offset[i]["latitude offset"]["value"] for i in range(len(names))]
    )

    # Calculate the offsets from the center.

    offsets = [
        np.sqrt(
            (offset[i]["longitude offset"]["value"] - mean_longitude) ** 2
            + (offset[i]["latitude offset"]["value"] - mean_latitude) ** 2
        )
        for i in range(len(names))
    ]

    # Calculate the number of flags for each antenna.

    nflags = [
        tb.calc(
            "[select from "
            + vis
            + " where ANTENNA1=="
            + str(i)
            + " giving  [ntrue(FLAG)]]"
        )["0"].sum()
        for i in range(len(names))
    ]

    # Calculate a score based on those two.

    score = [
        offsets[i] / max(offsets) + nflags[i] / max(nflags) for i in range(len(names))
    ]

    # Print out the antenna scores.

    print("Antenna scores for " + vis)
    for i in np.argsort(score):
        print(names[i], score[i])

    # Return the antenna names sorted by score.

    return ",".join(np.array(names)[np.argsort(score)])

def fetch_scan_info(vis, field):
    obsdict = listobs(vis)
    scan_length = []
    integration_time = []
    for key in obsdict.keys():
        if "scan" in key and obsdict[key]["0"]["FieldName"] == field:
            scan_dt = (obsdict[key]["0"]["EndTime"] - obsdict[key]["0"]["BeginTime"]) * 86400.0 # in sec
            int_dt = obsdict[key]["0"]["IntegrationTime"] # in sec
            scan_length.append(scan_dt)
            integration_time.append(int_dt)
    return scan_length, integration_time

def get_solints(vis, field):
    scan_length, integration_time = fetch_scan_info(vis, field)
    if len(np.unique(np.round(integration_time, decimals=3))) != 1:
        print("Warning: Integration times are different among scans. Use the maximum value as a representative integration time.")
        integration_time = np.max(integration_time)
    else:
        integration_time = np.mean(integration_time)
    
    solint_list = ["inf"] # start from the scan length (inf)

    solint = np.max(scan_length) * 0.5
    while solint >= 2.0 * integration_time:
        ints_per_solint = solint / integration_time
        if not ints_per_solint.is_integer():
            delta = ints_per_solint - np.floor(ints_per_solint)
            solint -= delta * integration_time # to force the solint to be a integer multiple of the integration time
        solint_list.append("{:.2f}s".format(solint))
        solint *= 0.5
    
    # and finally with solint="int" (integration time)
    solint_list.append("int")

    print("Determined solints: ", solint_list)
    return solint_list


# def plot_gain(caltable, ):





# def get_station_numbers(msfile, antenna_name):
#     """
#     Get the station numbers for all observations in which the given antenna appears

#     Parameters
#     ==========
#     msfile: Name of measurement set (string)
#     antenna_name: Name of antenna (e.g. "DA48")
#     """
#     tb.open(msfile+'/ANTENNA')
#     ant_names = tb.getcol('NAME')
#     ant_stations = tb.getcol('STATION')
#     tb.close()

#     ant_numbers = np.where(ant_names == antenna_name)[0]

#     tb.open(msfile)
#     antenna1 = tb.getcol('ANTENNA1')
#     obsid = tb.getcol('OBSERVATION_ID')
#     tb.close()

#     for i in ant_numbers:
#         matching_obs = np.unique(obsid[np.where(antenna1==i)])
#         for j in matching_obs:
#             print "#Observation ID %d: %s@%s" % (j, antenna_name, ant_stations[i])


def fix_phasecenter(vis, outputvis, phase_center, common_dir, field):
    """
    Wrapper to shift the peaks to the phase center, and reassign
    the phase centers to a common direction.
    This is applied to all SUBMMS to force it to work on an Multi-MS.
    """

    os.system("rm -rf " + outputvis + "*")
    fixvis(vis=vis, outputvis=outputvis, field=field, phasecenter=phase_center)
    fixplanets(vis=outputvis, field=field, direction=common_dir)

    # for subvis in os.listdir(outputvis+'/SUBMSS'):
    #     subvisdir = outputvis+'/SUBMSS/'+subvis
    #     fixvis(vis=subvisdir, outputvis=subvisdir,field=field,phasecenter=phase_center)
    #     fixplanets(vis=subvisdir, field=field, direction=common_dir)


def plot_calibration_table(caltable, plotval="phase", sharey=False, show=False, outfile="gain.png"):

    tb.open(caltable)
    time = tb.getcol("TIME") # in sec
    start = np.min(time)
    time -= start
    antennas = tb.getcol("ANTENNA1")
    refant = np.unique(tb.getcol("ANTENNA2"))[0]
    solutions = tb.getcol("CPARAM").squeeze()
    flags = tb.getcol("FLAG").squeeze()
    tb.close()

    # fetch antenna name
    tb.open(caltable + "/ANTENNA")
    ant_name = tb.getcol("NAME")
    ant_sta = tb.getcol("STATION")
    tb.close()

    npanel = np.unique(antennas).size

    ncols = 6
    nrows = int(npanel/ncols*0.999) + 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3.8*ncols, 2.5*nrows), sharex=True, sharey=sharey)
    ymax = 0.0
    for i, (ax, ant) in enumerate(zip(axes.flatten(), np.unique(antennas))):
        t = time[antennas == ant]
        sol = solutions[antennas == ant]
        flag = flags[antennas == ant]
        color = np.array(["tab:blue"]*len(sol))
        color[flag] = "grey"

        if "ph" in plotval:
            val = np.angle(sol, deg=True)
        elif "amp" in plotval:
            val = np.abs(sol)

        ax.scatter(t, val, facecolors=color, alpha=0.5, edgecolors=None)
        if "ph" in plotval:
            ax.axhline(y=0.0, color="grey", ls="dashed")
        title = ant_name[ant] + "@" + ant_sta[ant]
        if ant == refant:
            title += " (reference)"
        ax.set_title(title, fontsize=10)


        # label stuff
        if i == int(ncols * (nrows - 1)):
            ylabel = "Gain phase (deg)" if "ph" in plotval else "Gain amp. (Jy)"
            ax.set(xlabel="$\Delta t$ (sec)", ylabel=ylabel)

        if "ph" in plotval:
            if sharey:
                ymax = np.max(np.abs(val)) if np.max(np.abs(val)) > ymax else ymax
            else:
                ymax = np.max(np.abs(val))
            ax.set(ylim=(-ymax*1.2, ymax*1.2))
    
    for j in range(i+1, len(axes.flatten())):
        ax = axes.flatten()[j]
        ax.set_axis_off()
        
    fig.suptitle(caltable)

    fig.savefig(outfile, bbox_inches="tight", pad_inches=0.01)

    if show:
        plt.show()

def plot_solution_SNR_hist(caltable, solint=None, combine=None, show=False, outfile="SNR.pdf"):

    tb.open(caltable)
    SNR = tb.getcol("SNR").squeeze()

    fig, ax = plt.subplots(figsize=(12, 9))
    label = "solint = '{:s}'".format(solint) if solint is not None else ""
    label += ", combine = {:s}".format(combine) if combine is not None else ""
    ax.hist(SNR, bins=100, alpha=0.3, label=label)
    ax.set(xlabel="S/N of solutions", ylabel="# of solutions")
    ax.legend()

    fig.savefig(outfile, bbox_inches="tight", pad_inches=0.01)

    if show:
        plt.show()

def run_selfcal_iteration(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=None,
    calmode="p",
    solint="inf",
    combine="",
    minsnr_p=2.0,
    minsnr_ap=4.0,
    minblperant=4,
    solnorm_ap=False,
    remove_following_iter=True,
    eval_image=True,
    disk_mask="",
    noise_mask="",
    tclean_kwargs=dict(),
    showplot=False,
):
    # prefix
    mode_str = "_" + mode if mode is not None else "" # used for namin convention
    data_params_selfcal_keys = ["mode", "table", "solint", "combine", "spwmap"]

    # remove products produced by the follwoing iterations
    if remove_following_iter:
        cmd = "rm -r " + processingdir + "*_selfcal" + mode_str + "_iter[{0:d}-99]*".format(iteration)
        subprocess.run(cmd, shell=True)
        cmd = "rm -r " + processingdir + "*_selfcal" + mode_str + "_final*"
        subprocess.run(cmd, shell=True)

    # get the list of target measurement sets
    vislist = []
    for i in data_params.keys():
        print("# " + i)
        if mode == "SB-only" and "LB" in i:
            print(f"Skip {i} since mode = 'SB-only'")
            continue
        vis = processingdir + prefix + f"_{i}_selfcal{mode_str}_iter{iteration}_{calmode}.ms"
        os.system("rm -r " + vis)
        if iteration == 1:
            # set data_params
            for key in data_params_selfcal_keys:
                data_params[i]["selfcal_" + key] = []

            # copy the original MS
            if mode == "combined" and "SB" in i:
                prevvis = data_params[i][visID + "_selfcal_SB-only"]
            else:
                prevvis = data_params[i][visID]
            print(f"MS to be calibrated: {prevvis}")
            os.system(f"cp -r {prevvis} {vis}")
        else:
            # update data_params
            for key in data_params_selfcal_keys:
                data_params[i]["selfcal_" + key] = data_params[i]["selfcal_" + key][0:iteration-1]

            # calmode of one step before, used to name the output products of the current iteration
            prevcalmode = data_params[i]["selfcal_mode"][-1]
            prevvis = processingdir + prefix + f"_{i}_selfcal{mode_str}_iter{iteration-1}_{prevcalmode}.ms"
            print(f"MS to be calibrated: {prevvis}")
            casatasks.split(vis=prevvis, outputvis=vis) # split out the corrected_data
            

        data_params[i][visID + "_selfcal"] = vis
        vislist.append(vis)

    # tclean to put the model visibilities to the MS
    if iteration == 1:
        imagename_prev = processingdir + prefix + f"_selfcal{mode_str}_initial"
    else:
        imagename_prev = processingdir + prefix + f"_selfcal{mode_str}_iter{iteration-1}_{prevcalmode}"
    tclean_continuum_wrapper(
        vis=vislist,
        imagename=imagename_prev,
        imsize=imsize,
        cellsize=cellsize,
        savemodel=tclean_kwargs.get("savemodel", "modelcolumn"),
        **tclean_kwargs
    )


    for i in data_params.keys():
        print("# " + i)
        if mode == "SB-only" and "LB" in i:
            print(f"Skip {i} since mode = 'SB-only'")
            continue
        
        # solve the gain
        print("Solving gain...")
        caltable = processingdir + prefix + f"_{i}_selfcal{mode_str}_iter{iteration}_{calmode}.g"

        # update data_params
        data_params[i]["selfcal_mode"].append(calmode)
        data_params[i]["selfcal_table"].append(caltable)
        data_params[i]["selfcal_solint"].append(solint)
        data_params[i]["selfcal_combine"].append(combine)

        # solve the gain 
        casatasks.gaincal(
            vis=data_params[i][visID + "_selfcal"],
            caltable=caltable,
            gaintype="T",
            refant=data_params[i]["reference_antennas"],
            calmode=calmode,
            solint=solint,
            spw=data_params[i]["spws"],
            minsnr=minsnr_p if calmode == "p" else minsnr_ap,
            minblperant=minblperant,
            combine=combine,
            solnorm=False if calmode == "p" else solnorm_ap,
        )

        # plotting stuffs
        print("Plotting gain solutions...")
        plot_solution_SNR_hist(
            caltable,
            solint=solint,
            combine=combine,
            outfile=caltable.replace(".g", "_SNR_hist.pdf"),
            show=showplot,
        )
        plot_calibration_table(
            caltable,
            plotval="phase" if calmode == "p" else "amplitude",
            outfile=caltable.replace(".g", "_gain.pdf"),
            show=showplot,
        )

        spwmap = (
            [0] * len(data_params[i]["spws"].split(",")) if "spw" in combine else [""]
        )
        # update data_params
        data_params[i]["selfcal_spwmap"].append(spwmap)

        # apply the solutions
        print("Applying solutions...")
        casatasks.applycal(
            vis=data_params[i][visID + "_selfcal"],
            gaintable=[caltable],
            interp="linearPD",
            calwt=True,
            spw=data_params[i]["spws"],
            spwmap=spwmap,
            applymode="calonly"
        )

        print(f"Iteration {iteration} of self-calibration (mode = '{calmode}') completed. Check the gain solutions.")
        print("Calibrated MS: " + vis)

    if eval_image:
        imagename = processingdir + prefix + f"_selfcal{mode_str}_iter{iteration}_{calmode}"
        tclean_continuum_wrapper(
            vis=vislist,
            imagename=imagename,
            imsize=imsize,
            cellsize=cellsize,
            **tclean_kwargs
        )

        print(f"###### pre selfcal iteration {iteration} ######")
        estimate_SNR(imagename_prev + ".image.tt0", disk_mask, noise_mask)
        print(f"###### post selfcal iteration {iteration} ######")
        snr, rms = estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)

        return snr, rms

def finalize_selfcal(
    data_params,
    processingdir,
    prefix,
    iteration,
    visID,
    imsize,
    cellsize,
    mode=None,
    calmode="ap",
    disk_mask="",
    noise_mask="",
    tclean_kwargs=dict(),
):
    # prefix
    mode_str = "_" + mode if mode is not None else "" # used for naming convention
    data_params_selfcal_keys = ["mode", "table", "solint", "combine", "spwmap"]

    # remove products produced by the follwoing iterations
    # os.system("rm -r *_selfcal" + mode_str + "_iter[{0:d}-99]*".format(iteration + 1))
    cmd = "rm -r " + processingdir + "*_selfcal" + mode_str + "_iter[{0:d}-99]*".format(iteration + 1)
    subprocess.run(cmd, shell=True)
    cmd = "rm -r " + processingdir + "*_selfcal" + mode_str + "_final*"
    subprocess.run(cmd, shell=True)

    # get the list of final measurement sets
    vislist = []
    for i in data_params.keys():
        print("# " + i)
        if mode == "SB-only" and "LB" in i:
            print(f"Skip {i} since mode = 'SB-only'")
            continue
        vis = processingdir + prefix + f"_{i}_selfcal{mode_str}_iter{iteration}_{calmode}.ms"
        outputvis = processingdir + prefix + f"_{i}_selfcal{mode_str}_final.ms"
        os.system("rm -r " + outputvis)
        print("Final self-calibrated MS: " + outputvis)
        casatasks.split(vis=vis, outputvis=outputvis, datacolumn="corrected")

        # backup data_params for each mode
        if mode is not None:
            data_params[i][visID + "_selfcal" + mode_str] = data_params[i][visID + "_selfcal"]
            for key in data_params_selfcal_keys:
                data_params[i]["selfcal_" + key + mode_str] = data_params[i]["selfcal_" + key].copy()
    
        # when mode = "combined", concatenate with the Sb-only metrics
        if mode == "combined" and "SB" in i:
            for key in data_params_selfcal_keys:
                data_params[i]["selfcal_" + key] = data_params[i]["selfcal_" + key + "_SB-only"] + data_params[i]["selfcal_" + key + "_combined"]

        data_params[i][visID + "_selfcal"] = outputvis
        vislist.append(outputvis)
    
    # final imaging
    imagename = processingdir + prefix + f"_selfcal{mode_str}_final"
    tclean_continuum_wrapper(
        vis=vislist,
        imagename=imagename,
        imsize=imsize,
        cellsize=cellsize,
        **tclean_kwargs
    )

    snr, rms = estimate_SNR(imagename + ".image.tt0", disk_mask, noise_mask)

    return snr, rms
    



def save_data_params(data_params, filename):
    with open(filename, "w") as f:
        json.dump(data_params, f, indent=4)

def load_data_params(filename):
    with open(filename, "r") as f:
        data_params = json.load(f)
    return data_params



