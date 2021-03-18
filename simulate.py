# Image filtering simulation
#
# Mikael Mieskolainen, 2021

import os
import copy
from tqdm import tqdm
from tools import *
import sys

def filter_and_compare(origfile, noisyfile, denoisedfile, M, N, threshold, step, scalenorm):

    def clip_norm(x):
        y = copy.deepcopy(x)
        y[y < 0] = 0
        y /= np.max(y[:])
        return y

    # Run image filter
    os.system(f"./dctfilter -q -i {'./out/' + noisyfile + '.bin'} -o {'./out/' + denoisedfile + '.bin'} -m {M} -n {N} -s {step} -t {threshold}")

    A_o = readimg(origfile + '.bin', M=M, N=N)
    A_n = readimg('./out/' + noisyfile + '.bin', M=M, N=N)
    A_d = readimg('./out/' + denoisedfile + '.bin', M=M, N=N)
    
    ## Max normalization of the scales
    if scalenorm:
        A_o = clip_norm(A_o)
        A_n = clip_norm(A_n)
        A_d = clip_norm(A_d)
    
    ## Compute metrics
    PSNR_n  = PSNR(ref=A_o, approx=A_n)
    PSNR_d  = PSNR(ref=A_o, approx=A_d)

    return PSNR_n, PSNR_d

def addnoise(orig, noisetype, sigma, clip):
    """ Helper function
    """
    if (noisetype == 'gaussian'):
        noisy     = add_gaussian_noise(img=orig, sigma=sigma, clip=clip)
    elif (noisetype == 'poisson'):
        noisy     = add_poisson_noise(img=orig, clip=clip)
    else:
        raise Exception(f'Unknown noise type {noisetype}')
    return noisy

def main():
    
    save_images = True

    # Seed the engine
    np.random.seed(123456)

    # Input
    M           = 512
    N           = 512
    origfile    = f'lena_{M}_{N}'

    # Number of random repetitions in the table values
    N_repeat    = 10

    # Filter parameters
    step        = 1
    tval        = np.linspace(0,150,150)  # Threshold values
    
    # Simulation
    noisetype   = 'gaussian'
    clip        = True
    sigma_val   = np.arange(5,55,5)       # Noise std values
    
    # Comparison metrics parameters
    scalenorm   = True
    
    
    # Text output
    txtfile     = f'./out/sim_{origfile}.out'
    OF          = open(txtfile, 'w')
    def dprint(text, end='\n'):
        """ Dual print to tex and stdout """
        print(text, end=end)
        OF.write(text + end)
    

    dprint('| sigma | PSNR +- std (noisy) [dB] | PSNR +- std (denoised) [dB] | optimal threshold |')
    dprint('|---|---|---|---|')


    # Threshold curves
    X = np.zeros((len(tval), len(sigma_val)))

    for k in range(len(sigma_val)):

        sigma = sigma_val[k]

        # Read image
        noisyfile = f'lena_{M}_{N}_noisy_{sigma:0.0f}'
        orig      = readimg(filename=origfile + '.bin', M=M, N=N)

        # Add noise
        noisy = addnoise(orig=orig, sigma=sigma, noisetype=noisetype, clip=clip)
        writeimg(img=noisy, filename='./out/' + noisyfile + '.bin')

        # -----------------------------------------------------
        # Noise threshold (oracle) loop
        psnrval      = np.zeros(len(tval))
        denoisedfile = f'lena_{M}_{N}_denoised_{sigma:0.0f}'
        for i in tqdm(range(len(tval)), file=sys.stdout):

            PSNR_n, PSNR_d = filter_and_compare(origfile=origfile, noisyfile=noisyfile,
                denoisedfile=denoisedfile, M=M, N=N, threshold=tval[i], step=step, scalenorm=scalenorm)
            psnrval[i]     = PSNR_d

        X[:,k] = psnrval

        # -----------------------------------------------------
        # Filter with PSNR optimal threshold

        best_t_ind        = np.argmax(psnrval)
        optimal_threshold = tval[best_t_ind]

        PSNR_n = np.zeros(N_repeat)
        PSNR_d = np.zeros(N_repeat)
        for i in range(N_repeat):

            # Add noise
            noisy = addnoise(orig=orig, sigma=sigma, noisetype=noisetype, clip=clip)
            writeimg(img=noisy, filename='./out/' + noisyfile + '.bin')

            PSNR_n[i], PSNR_d[i] = filter_and_compare(origfile=origfile, noisyfile=noisyfile,
                denoisedfile=denoisedfile, M=M, N=N, threshold=optimal_threshold, step=step, scalenorm=scalenorm)

        dprint(f'|{sigma:0.1f} | {np.mean(PSNR_n):0.2f} +- {np.std(PSNR_n):0.2f} | {np.mean(PSNR_d):0.2f} +- {np.std(PSNR_d):0.2f} | {optimal_threshold:0.1f} |')

        # -----------------------------------------------------
        A_o = readimg(origfile + '.bin', M=M, N=N)
        A_n = readimg('./out/' + noisyfile + '.bin', M=M, N=N)
        A_d = readimg('./out/' + denoisedfile + '.bin', M=M, N=N)

        ## Save images
        if save_images:

            im = Image.fromarray(A_o).convert("L")
            im.save(f'./img/{origfile}.png')

            im = Image.fromarray(A_n).convert("L")
            im.save(f'./img/{noisyfile}.png')

            im = Image.fromarray(A_d).convert("L")
            im.save(f'./img/{denoisedfile}.png')

        ## Plot
        fig,ax = plt.subplots(ncols=3, nrows=1, figsize=(25,35))
        cmap = 'viridis'

        ax[0].imshow(A_o, cmap=plt.get_cmap(cmap))
        ax[1].imshow(A_n, cmap=plt.get_cmap(cmap))
        ax[2].imshow(A_d, cmap=plt.get_cmap(cmap))

        ax[0].set_title('Original')
        ax[1].set_title(f'$\\sigma = {sigma:0.0f}$: Noisy PSNR = {np.mean(PSNR_n):0.2f} ')
        ax[2].set_title(f'Denoised PSNR = {np.mean(PSNR_d):0.2f}')

        #plt.show()
        plt.savefig(f'./img/compare_sigma_{sigma:0.0f}.png', bbox_inches='tight')
        plt.close()
    
    OF.close()

    # --------------------------------------------------------------------
    # Plot treshold behavior

    fig = plt.figure()

    for i in range(X.shape[1]):
        plt.plot(tval, X[:,i], label=f'$\\sigma = {sigma_val[i]:0.0f}$')

    plt.xlabel('threshold')
    plt.ylabel('PSNR (denoised) [dB]')
    plt.tight_layout()
    plt.legend(loc='upper right')
    plt.xlim([np.min(tval), np.max(tval)])
    plt.savefig(f'./img/threshold_curve.png', bbox_inches='tight')
    plt.close()


if __name__ == '__main__' :
   main()
