from brian2 import *
# https://www.jneurosci.org/content/jneuro/32/41/14064.full.pdf
# ^^ explanation of HH model wrt computational neuroscience


def simulate_HH_network(A,stimulus, duration, dt=.01, post_noise_sigma=0):
    '''
    #A is adjacency matrix
    stimulus is a list of tuples (time, neuron, current, length) [time/length have to be % by .5ms]
    duration in ms
    dt in ms (optional)
    post_noise_sigma is added after simulation variance (optional)

    returns: statemonitor with t and v
    '''
    
    # define parameters if not given
    parameters = {
        # Boltzmann function parameters
        'v_n_half': 12*mV,
        'v_m_half': 25*mV,
        'v_h_half': 3*mV,

        'k_n': 15*mV,
        'k_m': 9*mV,
        'k_h': -7*mV,

        # Gaussian function parameters
        'v_n_max': -14*mV,
        'v_m_max': 27*mV,
        'v_h_max': -2*mV,

        'sigma_n': 50*mV,
        'sigma_m': 30*mV,
        'sigma_h': 20*mV,

        'c_n_amp': 4.7*ms,
        'c_m_amp': 0.46*ms,
        'c_h_amp': 7.4*ms,

        'c_n_base': 1.1*ms,
        'c_m_base': 0.04*ms,
        'c_h_base': 1.2*ms,

        # conductances
        'g_K_bar': 36*mS / (cmeter**2),
        'g_Na_bar': 120*mS / (cmeter**2),
        'g_L': 0.3*mS / (cmeter**2),

        # reversal potentials
        'e_K': -12*mV,
        'e_Na': 120*mV,
        'e_L': 10.6*mV,

        # membrane capacitance
        'C_mem': 1*uF / cmeter**2,

        # initial membrane voltage
        'v_initial': 0*mV,

        # initial gating variable activations
        'm_initial': 0.05,
        'n_initial': 0.32,
        'h_initial': 0.60,

        #noise variance
        'sigma': 0, # i think we should add noise after simulation to avoid false fires and refactory periods affecting connections

    }
    # put arguments into parameters dictionary
    parameters["A"] = A
    parameters["defaultclock_dt"] = dt*ms
    parameters["duration"] = duration*ms
    parameters["neuron_count"] = A.shape[0]

    # define stimulus from custom input
    max_stim_time = np.max(stimulus[:,1]+stimulus[:,2])
    brian_stim = np.zeros([parameters["neuron_count"], int(max_stim_time/.5+1)])
    for stim in stimulus:
        start_ind = int(stim[1]/.5)
        end_ind = int(start_ind+stim[2]/.5)
        brian_stim[int(stim[0]), start_ind:end_ind] = stim[3]
    parameters["I_stim"] = TimedArray(np.transpose(brian_stim)*amp/meter**2, dt=.5*ms)

    eqs_HH = '''
        dv/dt = (I_ext - I_K - I_Na - I_L)/C_mem+sigma*xi*sqrt(1/second)*volt : volt
        I_K = g_K*(v - e_K) : ampere/meter**2
        I_Na = g_Na*(v - e_Na) : ampere/meter**2
        I_L = g_L*(v - e_L) : ampere/meter**2
        I_ext = I_stim(t,i)+I_syn : ampere/meter**2
        I_syn : ampere/meter**2
        g_K = g_K_bar*n**4 : siemens/meter**2
        g_Na = g_Na_bar*m**3*h : siemens/meter**2

        dm/dt = (m_inf - m)/tau_m : 1
        dn/dt = (n_inf - n)/tau_n : 1
        dh/dt = (h_inf - h)/tau_h : 1

        m_inf = 1/(1+exp((v_m_half-v)/k_m)) : 1
        n_inf = 1/(1+exp((v_n_half-v)/k_n)) : 1
        h_inf = 1/(1+exp((v_h_half-v)/k_h)) : 1

        tau_m = c_m_base + c_m_amp*exp(-(v_m_max - v)**2/sigma_m**2) : second
        tau_n = c_n_base + c_n_amp*exp(-(v_n_max - v)**2/sigma_n**2) : second
        tau_h = c_h_base + c_h_amp*exp(-(v_h_max - v)**2/sigma_h**2) : second
    '''
    # define neuron group based on paramters and HH model
    group = NeuronGroup(parameters['neuron_count'], eqs_HH, method='euler',namespace=parameters)

    # initial conditions
    group.v = randn(parameters["neuron_count"])*mV + parameters["v_initial"]
    group.m = parameters["m_initial"]
    group.n = parameters["n_initial"]
    group.h = parameters["h_initial"]

    # create synapse group based on synaptic current being proportional to input voltage
    eqs_syn = '''
                I_syn_post = weight*v_pre*siemens/meter**2 : ampere/meter**2 (summed)
                weight : 1
            '''
    syn = Synapses(group, group, eqs_syn)

    # create synaptic connections based on adjacency matrix
    sources, targets = parameters["A"].nonzero()
    syn.connect(i=sources, j=targets)
    weights = np.zeros(len(sources))
    for ii in range(len(sources)):
        weights[ii] = parameters["A"][sources[ii], targets[ii]]
    syn.weight[:] = weights

    # run simulation
    statemon = StateMonitor(group, ['v','I_ext'], record=True)
    defaultclock.dt = parameters["defaultclock_dt"]
    run(parameters["duration"])
    v = statemon.v
    t = statemon.t
    I_ext = statemon.I_ext
    v += np.random.normal(0,post_noise_sigma, v.shape)*volt

    v /= volt
    t /= second
    I_ext /= amp
    output = np.concatenate((np.asmatrix(t),v,I_ext))

    return output