    def single_run(L,T,nequil,nmcs,seed=None):
    from f2py_jit import jit
    import matplotlib.pyplot as plt
    import numpy as np
    # import os
    # pwd = os.path.dirname(__file__)


    filename = 'ising.f90'
    module   = 'common.f90'

    # filename = os.path.join(pwd, filename)
    # module   = os.path.join(pwd, module)
    f90 = jit([module, filename])

    if seed is None: seed = np.random.random(8)

    eave,e2ave,mave,m2ave,abs_mave,cv,chi,acceptance_prob = f90.ising(L,T,nequil,nmcs,seed)
    return {
        'eave'              : eave,
        'e2ave'             : e2ave,
        'mave'              : mave,
        'm2ave'             : m2ave,
        'abs_mave'          : abs_mave,
        'cv'                : cv,
        'chi'               : chi,
        'acceptance_prob'   : acceptance_prob}



def prova():
    from f2py_jit import jit
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    pwd = os.path.dirname(__file__)
    from usemodule import plot_val_over_NMC, visual_lattice, find_Tc


    # filename = 'ising.f90'
    module   = 'common.f90'

    # filename = os.path.join(pwd, filename)
    module   = os.path.join(pwd, module)
    f90 = jit([module])

    L = 30
    N = L*L
    T = 2
    nequil  = int(1e4)
    # nequil  = 0
    nmcs    = int(1e5)
    # nequil  = int(nmcs/10)
    seed = [0,0,0,0,0,0,0,0]
    seed = np.random.random(8)


    cum = f90.common.initial(L,T,nequil,nmcs,seed)
    energies = []
    magnetis = []
    absmagns = []
    for i in range(nequil):
        f90.common.metropolis()
    for i in range(nmcs):
        f90.common.metropolis()
        f90.common.data(cum)
        if (i%100 == 0):
            if (i%int(nmcs/10) == 0): print(f'{int(i/nmcs*100):d}%')
            # energies.append(cum[0]/N/i)
            # magnetis.append(cum[2]/N/i)
            absmagns.append(cum[4]/N/i)
    lattice = f90.common.spin

    plot_val_over_NMC(absmagns,x_multipl=100)
    plt.show()

    visual_lattice(lattice)
    plt.show()
    f90.common.finish()



def ising_simulation():
    import numpy as np
    import matplotlib.pyplot as plt
    import pantarei as rei
    from usemodule import plot_val_over_NMC, visual_lattice, find_Tc, fit_data

    job = rei.Parallel(single_run)
    # job = rei.Parallel(single_run,clear_first=True)

    Tcs = []



    nequil = int(1e4)
    # nequil = 0
    nmcs = int(1e5)

    Ts      = np.linspace(2,2.7,20)
    Ts      = np.linspace(2,2.7,24)
    Ts      = np.linspace(1,4,24)

    Ts      = np.linspace(2.2,2.6,24)
    L = 4
    es      = []
    ms      = []
    cvs     = []
    chis    = []
    for T in Ts:
        data = job(L=L,T=T,nequil=nequil,nmcs=nmcs)
        es.append(data['eave'])
        ms.append(data['mave'])
        cvs.append(data['cv'])
        chis.append(data['chi'])
        print(f'T={T:.2f}')
    plt.plot(Ts,cvs, 'k.-', label=f'L = {L}')
    x,y = fit_data(Ts,cvs)
    plt.plot(x,y,'k-')
    Tcs.append(find_Tc(Ts,cvs))

    Ts      = np.linspace(2.2,2.4,24)
    L = 15
    es      = []
    ms      = []
    cvs     = []
    chis    = []
    for T in Ts:
        data = job(L=L,T=T,nequil=nequil,nmcs=nmcs)
        es.append(data['eave'])
        ms.append(data['mave'])
        cvs.append(data['cv'])
        chis.append(data['chi'])
        print(f'T={T:.2f}')
    plt.plot(Ts,cvs, 'r.-', label=f'L = {L}')
    x,y = fit_data(Ts,cvs)
    plt.plot(x,y,'r-')
    Tcs.append(find_Tc(Ts,cvs))


    L = 30
    es      = []
    ms      = []
    cvs     = []
    chis    = []
    for T in Ts:
        data = job(L=L,T=T,nequil=nequil,nmcs=nmcs)
        es.append(data['eave'])
        ms.append(data['mave'])
        cvs.append(data['cv'])
        chis.append(data['chi'])
        print(f'T={T:.2f}')
    plt.plot(Ts,cvs, 'b.-', label=f'L = {L}')
    x,y = fit_data(Ts,cvs)
    plt.plot(x,y,'b-')
    Tcs.append(find_Tc(Ts,cvs))


    L = 50
    es      = []
    ms      = []
    cvs     = []
    chis    = []
    for T in Ts:
        data = job(L=L,T=T,nequil=nequil,nmcs=nmcs)
        es.append(data['eave'])
        ms.append(data['mave'])
        cvs.append(data['cv'])
        chis.append(data['chi'])
        print(f'T={T:.2f}')
    plt.plot(Ts,cvs, 'm.-', label=f'L = {L}')
    x,y = fit_data(Ts,cvs)
    plt.plot(x,y,'m-')
    Tcs.append(find_Tc(Ts,cvs))

    # cv_derivative = []
    # for i in range(len(es)-2):
    #     cv_derivative.append((es[i+2]-es[i])/(Ts[i+2]-Ts[i]))
    # plt.plot(Ts,cvs,label='energy variation')
    # plt.plot(Ts[1:-1],cv_derivative,label='energy derivative')
    plt.title('cv over T')

    # deg = 4
    # coeff = np.polyfit(Ts,cvs,deg)
    # p = np.poly1d(coeff)
    # x  = np.linspace(Ts[0],Ts[-1],1000)
    # y = p(x)


    # from scipy.optimize import curve_fit
    # x = Ts
    # y = cvs
    # n = len(x)
    # # mean = sum(x*y)/n
    # # sigma = sum(y*(x-mean)**2)/n
    #
    # def gaus(x,a,x0,sigma):
    #     return a*np.exp(-(x-x0)**2/(2*sigma**2))
    #
    # def lorentzian( x, a, x0, gam ):
    #     return a * gam**2 / ( gam**2 + ( x - x0 )**2)
    #
    # # popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
    # x0 = x[np.argmax(y)]
    # gamma = 0.3
    # popt,pcov = curve_fit(lorentzian,x,y,p0=[1,x0,gamma])
    #
    # x  = np.linspace(Ts[0],Ts[-1],1000)
    # plt.plot(x,lorentzian(x,*popt),'r',label='fit')
    #
    # # plt.plot(x,y)
    #
    # Tc = x[np.argmax(y)]
    # Tc = x[np.argmax(lorentzian(x,*popt))]
    # plt.scatter(Tc,np.max(lorentzian(x,*popt)),label=f'Tc = {Tc:.3f} J/kb')
    #
    plt.grid()
    plt.legend()
    plt.show()
    # rei.report()

    Ls = [4,15,30,50]
    plt.loglog(Ls,Tcs,'.')
    # x,y = fit_data(Ls,[1/a for a in Tcs],deg=1)
    plt.show()




if __name__ == "__main__":
    ising_simulation()

