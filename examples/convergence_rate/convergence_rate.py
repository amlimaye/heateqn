import os
import sys
import importlib
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

def plot_results(results, outpath):
    final_times = [r[0] for r in results]
    dts = [[elem[0] for elem in r[1]] for r in results]
    gerrs = [[elem[1] for elem in r[1]] for r in results]

    plt.figure()
    for final_time, dt, gerr in zip(final_times, dts, gerrs):
        plt.plot(dt, gerr, 'o-', label='t=%0.2f' % final_time)

    plt.legend()
    plt.xlabel(r'$\delta t$')
    plt.ylabel(r'$|y^e(t) - y(t)|$')
    plt.tight_layout()
    plt.savefig(outpath)

def main(scriptname):
    mod = load_module(scriptname)
    outpath = os.path.join(os.path.dirname(__file__), 'out.png')

    results = mod.run()
    plot_results(results, outpath)

if __name__ == "__main__":
    main(sys.argv[0])
