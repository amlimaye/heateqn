import os
import sys
import importlib
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def check_shlib_existence(scriptname):
    shlib_name = 'pyexample_%s.so' % (os.path.basename(scriptname).split('.py')[0])
    shlib_path = os.path.join('.', shlib_name)
    print(shlib_path)

    import_name = shlib_name[:-3]

    exists = True

    return exists, import_name

def load_module(scriptname):
    exists, import_name = check_shlib_existence(scriptname)
    if not exists:
        raise RuntimeError('Could not locate shared library')

    mod = importlib.import_module(import_name)
    return mod

def get_scaling(x, y):
    logx = np.log(x)
    logy = np.log(y)

    assert logx.shape == logy.shape
    first_half = int(logx.shape[0]/2.0)

    p = np.polyfit(logx[:first_half], logy[:first_half], 1)
    slope = p[0]

    return slope

def plot_results(results, outpath):
    final_times = [r[0] for r in results]
    dts = [[elem[0] for elem in r[1]] for r in results]
    gerrs = [[elem[1] for elem in r[1]] for r in results]

    exponents = []
    plt.figure()
    for final_time, dt, gerr in zip(final_times, dts, gerrs):
        plt.subplot(2, 1, 1)
        plt.plot(dt, gerr, 'o-', label='t=%0.2f' % final_time)
        plt.subplot(2, 1, 2)
        plt.loglog(dt, gerr, 'o-', label='t=%0.2f' % final_time)
        slope = get_scaling(dt, gerr)
        exponents.append(slope)

    expo = np.mean(exponents)
    expo_unc = np.std(exponents)

    plt.subplot(2, 1, 1)
    plt.legend()
    plt.xlabel(r'$\delta t$')
    plt.ylabel(r'$|y^e(t) - y(t)|$')
    plt.subplot(2, 1, 2)
    plt.text(0.05, 0.85, r'$\gamma = %0.2f \pm %0.2f$' % (expo, expo_unc),
             transform=plt.gca().transAxes)
    plt.xlabel(r'$\delta t$')
    plt.ylabel(r'$|y^e(t) - y(t)|$')

    plt.tight_layout()
    plt.savefig(outpath)

def main(scriptname):
    mod = load_module(scriptname)

    fwd_outpath = os.path.join(os.path.dirname(__file__), 'fwd_euler.png')
    bkwd_outpath = os.path.join(os.path.dirname(__file__), 'bkwd_euler.png')

    fwd_results = mod.forward_euler()
    bkwd_results = mod.backward_euler()

    plot_results(fwd_results, fwd_outpath)
    plot_results(bkwd_results, bkwd_outpath)

if __name__ == "__main__":
    main(sys.argv[0])
