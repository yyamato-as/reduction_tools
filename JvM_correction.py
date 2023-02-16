# script to apply JvM correction
# author: Y. Yamato
# based on the original MAPS scropt (JvM_correction_casa6.py) 
# and Ian Czekala's MPol package for first null determination 
# (see https://github.com/MPoL-dev/MPoL/blob/main/src/mpol/gridding.py#L441) 
# also based on a script fron Gianni Cataldi

import shutil
import numpy as np
import casatools
import casatasks
import astropy.io.fits as fits
import os

pi = np.pi
rad_to_arcsec = 180 * 3600 / pi

def FWHM_to_sigma(FWHM):
    return FWHM / np.sqrt(8 * np.log(2))


def get_beam(header, verbose=False):
    try:
        beam_info = header["restoringbeam"]
    except KeyError:
        beam_info = header["perplanebeams"]["beams"]["*0"]["*0"]

    assert (
        beam_info["major"]["unit"] == beam_info["minor"]["unit"] == "arcsec"
    ), "Unit of beam is not arcsec."
    assert (
        beam_info["positionangle"]["unit"] == "deg"
    ), "Unit of beam position angle is not degree."

    bmaj = beam_info["major"]["value"]
    bmin = beam_info["minor"]["value"]
    bpa = beam_info["positionangle"]["value"]

    if verbose:
        print(
            "The CASA CLEAN beam is {:.3f} arcsec x {:.3f} arcsec ({:.2f} deg)".format(
                bmaj, bmin, bpa
            )
        )

    return bmaj, bmin, bpa

def get_axis(header, i=0, rpix=None):
    npix = header["shape"][i]
    dpix = header["incr"][i]
    if rpix is None:
        rpix = header["refpix"][i]  # casa image uses zero-based counting

    ax = dpix * (np.arange(npix) - rpix)

    return ax


def get_polar_coord(header, refpix=(None, None)):
    x = get_axis(header, i=0, rpix=refpix[0])
    y = get_axis(header, i=1, rpix=refpix[1])
    xx, yy = np.meshgrid(x, y, indexing="ij")
    r = np.sqrt(xx**2 + yy**2)
    theta = np.arctan2(yy, xx)
    return r, theta


def get_psf_peak_index(psf_array):
    return np.unravel_index(np.argmax(psf_array), psf_array.shape)


def get_dpix(header):
    assert np.abs(header["incr"][0]) == np.abs(
        header["incr"][1]
    ), "Image is not square."
    assert header["axisunits"][0] == header["axisunits"][1] == "rad", "Unit of the axis is not rad." 
    return np.abs(header["incr"][0]) * rad_to_arcsec


def get_nulled_psf(psf_array, r, theta, ntheta=24):

    # azimuthal coordinate
    da = 2 * pi / ntheta
    azim = np.arange(-pi, pi, da)

    nulled_psf = psf_array.copy()

    for a in azim:
        # some overlapping azimuth: see https://github.com/MPoL-dev/MPoL/blob/main/src/mpol/gridding.py#L441
        wedge = (theta >= a - 0.3 * da) & (theta <= a + 1.3 * da)
        wedge_is_neg = wedge & (psf_array < 0)
        first_null = np.min(r[wedge_is_neg])
        nulled_psf[wedge & (r >= first_null)] = 0.0

    return nulled_psf


def get_JvM_epsilon(psfname, ntheta=24):

    # get the psf array and header
    ia = casatools.image()
    ia.open(psfname)
    ref_chan_ind = int(ia.shape()[3]/2)
    header = ia.summary(list=False)
    ia.close()

    # workaround for memory error; loading the full cube sometimes throw a memory error for a large cube, 
    # thus instead do imsubimage to get a singe channel image which includes the representative psf image 
    # (assume psf is same along the spectral axis)
    psfname_temp = psfname.replace(".psf", "_temp.psf")
    casatasks.imsubimage(imagename=psfname, outfile=psfname_temp, chans=str(ref_chan_ind), overwrite=True)
    ia.open(psfname_temp)
    psf_array = ia.getregion().squeeze()
    ia.close()
    os.system("rm -r " + psfname_temp)

    if psf_array.ndim > 2:
        ref_chan_ind = int(psf_array.shape[2]/2) # assume psf is same along the spectral axis and use the nearly center channel; sometimes edge channel have not a good representation of psf over the channels? 
        psf_array = psf_array[:, :, ref_chan_ind]

    # get peak index of psf and polar coordinate
    peak_pix = get_psf_peak_index(psf_array)
    r, theta = get_polar_coord(header, refpix=peak_pix)

    # determine the null points in each azimuthal direction
    nulled_psf = get_nulled_psf(psf_array, r, theta, ntheta=ntheta)

    # get CLEAN beam info
    bmaj, bmin, _ = get_beam(header, verbose=True)

    Omega_cbeam = 2 * pi * FWHM_to_sigma(bmaj) * FWHM_to_sigma(bmin)
    print(f"clean beam area: {Omega_cbeam} arcsec2")
    Omega_psf = np.sum(nulled_psf) * get_dpix(header) ** 2
    print(f"psf area (inside first null): {Omega_psf} arcsec2")
    epsilon = Omega_cbeam / Omega_psf
    print(f"JvM epsilon: {epsilon}")

    return epsilon


def apply_JvM_correction(prefix, mfs=True, ntheta=24, pbcor=True):

    psfname = prefix + ".psf"
    modelname = prefix + ".model"
    residualname = prefix + ".residual"

    if mfs:
        psfname += ".tt0"
        modelname += ".tt0"
        residualname += ".tt0"

    # get JvM epsilon
    print("Calculating JvM epsilon...")
    epsilon = get_JvM_epsilon(psfname, ntheta=ntheta)
    print("Done.")

    # get beam info
    header = casatasks.imhead(psfname, mode="summary")
    bmaj, bmin, bpa = get_beam(header, verbose=False)

    print(f"Applying JvM correction to {prefix}...")
    # convolve the model image temporary
    model_convolved = prefix + "_model_convolved_temp.image"
    casatasks.imsmooth(
        imagename=modelname,
        major=f"{bmaj}arcsec",
        minor=f"{bmin}arcsec",
        pa=f"{bpa}deg",
        targetres=True,
        outfile=model_convolved,
        overwrite=True,
    )

    # apply JvM correction
    image_JvMcorr = prefix + ".JvMcorr.image"
    if mfs:
        image_JvMcorr += ".tt0"
    if os.path.exists(image_JvMcorr):
        shutil.rmtree(image_JvMcorr)
    casatasks.immath(
        imagename=[model_convolved, residualname],
        expr=f"IM0 + {epsilon}*IM1",
        imagemd=model_convolved,
        outfile=image_JvMcorr,
    )
    shutil.rmtree(model_convolved)
    corrected_image = [image_JvMcorr]

    print(f"JvM-corrected image is stored in {image_JvMcorr}.")

    print(f"Applying lowres JvM correction to {prefix}...")
    # convolve the model image temporaliry
    model_convolved = prefix + "_model_convolved_temp.image"
    casatasks.imsmooth(
        imagename=modelname,
        major=f"{bmaj/np.sqrt(epsilon)}arcsec",
        minor=f"{bmin/np.sqrt(epsilon)}arcsec",
        pa=f"{bpa}deg",
        targetres=True,
        outfile=model_convolved,
        overwrite=True,
    )

    # apply "lowres" JvM correction
    image_JvMcorr_lowres = prefix + ".JvMcorr.lowres.image"
    if mfs:
        image_JvMcorr_lowres += ".tt0"
    if os.path.exists(image_JvMcorr_lowres):
        shutil.rmtree(image_JvMcorr_lowres)
    casatasks.immath(
        imagename=[model_convolved, residualname],
        expr=f"IM0 + IM1",
        imagemd=model_convolved,
        outfile=image_JvMcorr_lowres,
    )
    shutil.rmtree(model_convolved)
    corrected_image.append(image_JvMcorr_lowres)

    print(f"lowres JvM-corrected image is stored in {image_JvMcorr}.")

    pbcorrected_image = []
    if pbcor:
        for image in corrected_image:
            image_JvMcorr_pbcor = image.replace(".image", ".pbcor")
            pbimage = prefix + ".pb.tt0" if mfs else prefix + ".pb"
            if os.path.exists(image_JvMcorr_pbcor):
                shutil.rmtree(image_JvMcorr_pbcor)
            casatasks.impbcor(
                imagename=image, pbimage=pbimage, outfile=image_JvMcorr_pbcor
            )
            pbcorrected_image.append(image_JvMcorr_pbcor)
            print(f"JvM-corrected, pb-corrected image is stored in {image_JvMcorr_pbcor}.")

    to_export = corrected_image + pbcorrected_image

    # export to fits
    for imagename in to_export:
        casatasks.exportfits(
            imagename=imagename,
            fitsimage=imagename + ".fits",
            dropstokes=True,
            overwrite=True,
        )
        # add epsilon to FITS header
        fits.setval(filename=imagename + ".fits", keyword="JVM_EPS", value=epsilon)
        print(f"{imagename} is exported to a FITS.")

    print(f"JvM correction done for {prefix}.")


# def plot_psf_profile(psfname, rmax=2, ax=None):

#     if not psfname.endswith(".fits"):
#         casatasks.exportfits(
#             imagename=psfname, fitsimage=psfname + "_temp.fits", dropdeg=True
#         )
#         psfname += "_temp.fits"

#     psfimage = FitsImage(psfname, xlim=(-rmax, rmax), ylim=(-rmax, rmax))
#     bincl = np.rad2deg(np.arccos(psfimage.bmin / psfimage.bmaj))
#     rbins = np.arange(0.0, rmax, psfimage.dpix)

#     # psf profile
#     r, I, dI = psfimage.radial_profile(PA=psfimage.bpa, incl=bincl, rbins=rbins)

#     # clean beam profile
#     cbeam = np.exp(-(r**2) / (2 * FWHM_to_sigma(psfimage.bmaj) ** 2))

#     if ax is None:
#         fig, ax = plt.subplots()

#     ax.plot(r, I, label="psf (dirty beam)")
#     ax.fill_between(r, I - 3 * dI, I + 3 * dI, alpha=0.25)
#     ax.plot(r, cbeam, color="black", ls="dashed", label="clean beam")
#     ax.set(xlabel="Radius [arcsec]", ylabel="Relative profile")

#     if psfname.endswith("_temp.fits"):
#         os.system("rm " + psfname)

#     return ax


# if __name__ == '__main__':

#     path = "/works/yamato/eDisk/L1489IRS/archive/v1_data/custom_images/"
#     prefix = path + "L1489IRS_SBLB_C18O_robust_1.0"

#     apply_JvM_correction(prefix, ntheta=24, mfs=False)