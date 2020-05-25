import sys
sys.path.insert(0, './../../')
# from data_analysis import *

################################

def diffusion_map(filename):
    time_idx = 2
    particle_ids_idx = 1
    data_idx = [3,4]
    
    cell_ids,time_data,obs_data = read_observations(filename,
                                                    time_idx,
                                                    particle_ids_idx,
                                                    data_idx,
                                                    filling_values=1,
                                                    skip_header=1,
                                                    delimiter=",")
    
    print(cell_ids.shape)
    
    ###############################
    
    # build time lag vectors
    NLag = 30
    NLagStep =5
    NPmax = int(max(cell_ids))
    
    data,pdata,counter = build_timelag_vectors(obs_data, cell_ids, NLag = NLag, NLagStep = NLagStep, NPmax =NPmax)
    
    behaviors = np.zeros((data.shape[0],))
    
    ################################
    
    ind = np.sqrt(pdata[:,3]**2+pdata[:,4]**2) < 1
    pdata1 = pdata[ind,:]
    data1 = data[ind,:]
    
    visualize_data(pdata1)
    
    ################################
    
    n_evecs = 10
    epsilon = 5e-2
    state_dependent_eps = False
    
    evecs,evals,residuals,eps,density = reduce_dimension(data1,
                                                         epsilon = epsilon,
                                                         n_evecs=n_evecs,
                                                         state_dependent_eps=state_dependent_eps)
    
    ################################
    
    color = np.sqrt(pdata1[:,3]**2+pdata1[:,4]**2)
    
    visualize_reduced_dimension(evecs,residuals,color,PTS=[2,5],figsize=(8,3))
    
    ################################
    
    good_evecs = evecs[:,[1,2,3]]
    
    fig = plt.figure(figsize=(4,4))
    ax1 = fig.add_subplot(1,1,1, projection='3d')
    
    fig2,ax = plt.subplots(1,3,figsize=(8,2))
    
    minidx = 0
    x = good_evecs[:,minidx]
    y = good_evecs[:,minidx+1]
    z = good_evecs[:,minidx+2]
    
    w1 = color
    ax1.scatter(x,y,z, s=1, c=w1, cmap='jet')
    xlabel=ax1.set_xlabel('psi1')
    ylabel=ax1.set_ylabel('psi2')
    zlabel=ax1.set_zlabel('psi3')
    ax1.view_init(25,20)
    
    #fig.savefig(data_path + 'dmap_space_real1.pdf',bbox_extra_artists=[xlabel,ylabel,zlabel], bbox_inches='tight')
    
    cmap = 'jet'
    ax[0].scatter(pdata1[:,1],pdata1[:,2],s=10,c=good_evecs[:,0],cmap=cmap)
    ax[1].scatter(pdata1[:,1],pdata1[:,2],s=10,c=good_evecs[:,1],cmap=cmap)
    ax[2].scatter(pdata1[:,1],pdata1[:,2],s=10,c=good_evecs[:,2],cmap=cmap)
    
    plt.show() # Line added to show plot
    
    ################################
    
    r = radial_falloff(good_evecs[:,[0,1,2]])
    
    fig = plt.figure()
    
    bins = np.linspace(0,2,101)
    hist, bin_edges = np.histogram(r,bins=bins, density=True)
    bw = bin_edges[2]-bin_edges[1]
    plt.bar(bin_edges[:-1],hist,bw)
    plt.xlim([0,2])
    
    plt.show() # Line added to show plot



################################
# data_analysis.py #############
################################

import matplotlib

import matplotlib.pyplot as plt

import mpl_toolkits

from mpl_toolkits.mplot3d import Axes3D

import numpy as np

import matplotlib.animation as manimation

import matplotlib.tri as mtri

from scipy.integrate import odeint
from scipy.integrate import complex_ode
import scipy.spatial.distance
import scipy
import scipy.interpolate
import scipy.signal

# from dmap_sp import *
# from useful import *
# from riemannian_metric import *


verbosity_level = 10

def print_v(message,verbose=10,verbosity_level=10):
    if verbose <= verbosity_level:
        print(message)
        
def read_observations(data_path,time_idx,particle_ids_idx,data_idx,
                      filling_values=1, skip_header=1, delimiter=",",
                      verbose=10):
    """
    reads in the trajectory file and separates time steps, cell ids, and observations.
    data_path: path to the trajectory file.
    time_idx: index of the column in the file that contains the time information.
    particle_ids_idx: index of the column in the file that contains the id of the trajectories/particles.
    data_idx: indices of the columns that contain the observations (x,y positions).
    filling_values: if the data file contains NaN values, they will be replaced by this number.
    skip_header: 1 or 0, will skip the first row if set to 1.
    delimiter: this character has to be set to the delimiting character in the data file.
    verbose: if smaller than 11, prints output.
    """
    
    info_data = np.genfromtxt(data_path, missing_values="NaN",
                                filling_values=filling_values,
                                skip_header=skip_header,
                                delimiter=delimiter)
    info_data = np.array(info_data)
    print_v(" ".join(['data shape:', str(info_data.shape)]), verbose)
    
    cell_ids = info_data[:,particle_ids_idx]
    time_data = info_data[:,time_idx]
    obs_data = info_data[:,data_idx]

    return cell_ids,time_data,obs_data
    
def build_timelag_vectors(observations, traj_idx, NLag = 10, NLagStep = 5, NPmax = -1, onlyfirst=False, verbose=10):
    """
    builds time lag vectors from observation data.
    observations: the observation data, one row per observation, (x,y) values in the columns.
    traj_idx: the particle indices as returned by the read_observations method (cell_ids array).
    NLag: the length of the short trajectories that will be used in subsequent analysis.
    NLagStep: if onlyfirst=False, all trajectories will be split up into smaller segments, the starting position differs by NLagStep.
    NPmax: maximum number of rows in the observations to consider.
    onlyfirst: if set to True, will only the consider the first NLag (x,y) positions of each trajectory.
    verbose: if smaller than 10, will output text.
    """
    
    NLagStep = NLag//NLagStep
    
    data = np.zeros((observations.shape[0]//NLagStep//2,NLag*2))
    pdata = np.zeros((observations.shape[0],5))
    
    # print(data.shape)

    counter = 0
    if NPmax < 0:
        NP = int(np.max(traj_idx))
    else:
        NP = NPmax

    for i in range(1,NP): # loop over cells
        trajidx = traj_idx == i
        traj = observations[trajidx,:]
        traj = traj[0:NLag,:]
        
        if traj.shape[0] >= NLag:
            for k in range(1): #range(0,traj.shape[0]-NLag+1,NLagStep):
                cxy = np.mean(traj[:,k:(k+NLag)],axis=0)
                xy = np.array([[traj[n,0]-cxy[0], traj[n,1]-cxy[1]] for n in range(k,k+NLag)])
                sxy = np.std(xy,axis=0)
                data[counter,:] = xy.flatten()
                
                #if stat_function ~= -1:
                #    stat = stat_function(xy)
                #    pdata[counter,:] = stat
                #else:
                pdata[counter,0] = i
                pdata[counter,1] = cxy[0]
                pdata[counter,2] = cxy[1]
                pdata[counter,3] = sxy[0]
                pdata[counter,4] = sxy[1]
                counter = counter + 1

    data = data[0:(counter-1),:]
    pdata = pdata[0:(counter-1),:]
    behaviors = np.zeros((data.shape[0],))

    if counter == 0:
        print_v("ERROR: not enough data!!",verbose)

    print_v(data.shape)
    print_v(counter)
    # print_v(pdata)
    
    return data,pdata,counter

def visualize_data(pdata):
    """
    plots center positions, standard deviations, and number of points per trajectory in the data set.
    """
    fig,ax = plt.subplots(1,3,figsize=(6,3));
    ax[0].scatter(pdata[:,1], pdata[:,2],s=2,c='k')
    ax[0].set_xlabel('x [pixel]')
    ax[0].set_ylabel('y [pixel]')
    ax[0].set_aspect('equal')

    ax[1].scatter(pdata[:,3], pdata[:,4], s=5, c=(np.sqrt(pdata[:,3]**2+pdata[:,4]**2) > 50))
    ax[1].set_xlabel('sdev x [pixel]');
    ax[1].set_ylabel('sdev y [pixel]');
    ax[1].set_aspect('equal')

    ax[2].hist(pdata[:,0],30);
    ax[2].set_xlabel('cell#')
    ax[2].set_ylabel('#points')
    
    fig.tight_layout()

    plt.show() # Line added to show plot
    
def reduce_dimension(observations, n_evecs = 15, epsilon = 5e-2, verbosity=10, cutoff_dist=-1, state_dependent_eps=0):
    """
    computes an embedding of the given observations with Diffusion Maps.
    observations: the time lag vectors generated by the build_timelag_vectors method.
    n_evecs: dimension of the resulting embedding. Can be larger than the actual dimension of the intrinsic dataset.
    epsilon: bandwith parameter used to distinguish points in the observation space. Should be set to a value in [5e-4,1].
    verbosity: if smaller than 10, outputs text.
    cutoff_dist: if larger than 0, will only use the nearest neighbors in that distance for neighborhood computations.
    state_dependent_eps: if set to 1, changes epsilon intrinsically based on local density. Useful if large deviations in the observations density occur.
    """
    if observations.shape[0] < 10000:
        data_type = "compute dmatrix" # compute full distance matrix
    else:
        data_type = "points" # compute kd tree
    LB_flag = 1

    evecs,evals,density,M,distMatrix,eps = dmap_sp(observations,epsilon,n_evecs+1,LB_flag,data_type,
                                                  estimate_eps=True,
                                                   state_dependent_eps=state_dependent_eps,
                                                   exponent=2,
                                                   cutoff_num=0,
                                                   verbose=verbosity>5,
                                                   cutoff_dist=cutoff_dist)

    print_v('computing good embedding directions...', verbosity)

    eps_scale = 10
    res = compute_residuals_DMAPS(evecs,eps_scale)
    evecs = evecs[:,1::]
    evals = evals[1::]
    res = res[1::]
    
    for i in range(evecs.shape[1]):
        evecs[:,i] = evecs[:,i] - np.min(evecs[:,i])
        evecs[:,i] = (evecs[:,i] / np.max(evecs[:,i]))*2-1 # shift in [-1,1]
            
    print_v('done.', verbosity)
            
    return evecs,evals,res,eps,density

def visualize_reduced_dimension(evecs,residuals,color,PTS=[5,5],figsize=(8,8),evals=[-1]):
    """
    visualizes the results from the reduce_dimension method.
    """
    if len(evals) == 1:
        fig = plt.figure()
        plt.plot(residuals,'x')
    else:
        fig,ax = plt.subplots(1,2)
        ax[0].plot(residuals,'x')
        ax[0].set_title('residuals')
        ax[1].plot(evals,'x')
        ax[1].set_title('eigenvalues')
    good_evecs = evecs[:,residuals>.0]

    w1 = color

    cmin = np.min(w1)
    cmax = np.max(w1)

    fig,ax = plt.subplots(PTS[0],PTS[1],figsize=figsize,sharex=True,sharey=True)
    print(evecs.shape)
    for i in range(PTS[0]):
        for k in range(PTS[1]):
            ev1 = evecs[:,(i*PTS[1]+k)%evecs.shape[1]]
            ev2 = evecs[:,(i*PTS[1]+k+1)%evecs.shape[1]]
            sc1 = ax[i,k].scatter(ev1,ev2, 2, (w1), cmap='jet')
            ax[i,k].set_aspect('equal')
            ax[i,k].set_title(i*PTS[1]+k+1)
    fig.tight_layout()
    
    plt.show() # Line added to show plot

def radial_falloff(xyz):
    c = np.mean(xyz,axis=1)
    r = np.sqrt((xyz[:,0]-c[0])**2+(xyz[:,1]-c[1])**2+(xyz[:,2]-c[2])**2)
    return r

################################
# dmap_sp.py #############
################################


import numpy as np
import numpy.matlib
import scipy.spatial.distance
import scipy

import scipy
from time import time

import matplotlib
import matplotlib.pyplot as plt

def initialize_dmaps(data,epsilon_scale,n_evecs,LB_flag,data_type="points",estimate_eps=True,state_dependent_eps=0,exponent=2,cutoff_num=0,verbose=False,cutoff_dist=-1):
    if data_type == "points":
        tree = scipy.spatial.KDTree(data)
        if verbose:
            print('built kd tree.')
        t2 = time()

        d = data.shape[1]
        if cutoff_dist < 0:
            if verbose:
                print('max distance considered:', 10*epsilon_scale)
            distMatrix = tree.sparse_distance_matrix(tree, max_distance = 10*epsilon_scale, p = exponent)
        else:
            if verbose:
                print('max distance considered:', cutoff_dist)
            distMatrix = tree.sparse_distance_matrix(tree, max_distance = cutoff_dist, p = exponent)
        if verbose:
            print('built distance matrix in',int(time()-t2),'seconds.')
    elif data_type == "dmatrix":
        distMatrix = data
        D = data
    elif data_type == "compute dmatrix":
        distMatrix = scipy.spatial.distance.pdist(data, 'euclidean')
        distMatrix = scipy.spatial.distance.squareform(distMatrix)
        D = distMatrix
        distMatrix = scipy.sparse.csr_matrix(distMatrix)
    else:
        assert False,"data_type not recognized. exiting."
        
    if estimate_eps:
        epsilon = epsilon_scale * estimate_epsilon(distMatrix,exponent=exponent,cutoff_num=cutoff_num)
        print('estimated epsilon:', epsilon)
    else:
        epsilon = epsilon_scale
    
    csrdist = distMatrix.tocsr()
    A2 = csrdist.copy()
    
    if state_dependent_eps == 0:
        if verbose:
            print('epsilon',epsilon)
        epsilon_inv = 1/(4*epsilon)
        A = np.exp(-csrdist.data**exponent*(epsilon_inv))
        A2.data = A
    else:
        epsilon_inv = 1/(4.*epsilon)
        A = np.exp(-epsilon_inv*csrdist.data**exponent)
        A2.data = A
        D = np.array(A2.sum(axis=1)).flatten()
        
        print('density array', np.median(D))
        
        epsilon = (1/(D)**2) * epsilon
        epsilon_inv = scipy.sparse.csr_matrix(np.diag(1/(4*epsilon)))
        if verbose:
            print('epsilon',(np.median(epsilon)))
        A = -csrdist.data**exponent
        A2.data = A
        A2 = A2.dot(epsilon_inv)
        A2.data = np.exp(A2.data)
    A = A2
    
    if cutoff_num > 0:
        iDs = np.argsort(D,axis=1)
        for i in range(D.shape[1]):
            A[i,iDs[i,cutoff_num::]] = 0
    if cutoff_dist > 0:
        D = distMatrix.todense()
        A[D>cutoff_dist] = 0
        
    A = scipy.sparse.csr_matrix(A)
    #Adata = A.data
    #Adata[Adata < 1e-20] = 0
    #A.data = Adata
    
    if verbose:
        (x,y,z)=scipy.sparse.find(A)
        nsparse = np.array(z).shape[0]
        print('sparsity:', np.clip(int(10000*nsparse/(A.shape[0]**2))/10000,0,1))
    
    return A,distMatrix,epsilon

    
def SVD(data):
    # this is decomposing into "principal components"
    u,s,vh = np.linalg.svd(data)
    smat = np.zeros((u.shape[0],s.shape[0]), dtype=complex)
    smat[:s.shape[0], :s.shape[0]] = np.diag(s)
    return u,smat,vh
    
# dmaps sparse implementation, including k-d tree for distance matrix computation.
def dmap_sp(data,epsilon_scale,n_evecs,LB_flag,data_type="points",estimate_eps=True,state_dependent_eps=0,exponent=2,cutoff_num=0,verbose=False,cutoff_dist=-1):
    #
    # data          - data matrix, N rows (data points) with M columns (data dimension)
    # epsilon_scale - scale factor for the bandwidth of the Gaussian kernel.
    # n_evecs       - #evecs to compute by Arnoldi method
    # LB_flag       - 0 = FP normalization, 1 = LB (density free) normalization
    # data_type     - "points": treat as individual points and use a tree to compute the distances.
    #                 "dmatrix": treat as distance matrix.
    #                 "compute dmatrix": computes the matrix from data.
    # estimate_eps  - if True, epsilon_scale is multiplied with an estimation from data. False uses estimate_eps directly.
    # state_dependent_eps - if larger than 0, will use this number of nearest neighbors to define epsilon for each data point.
    # exponent      - exponent used in the similarity heat kernel. default is 2.
    # cutoff_num    - if larger than 0, and data_type="compute dmatrix", sets a distance cutoff after the given number of nearest neighbors.

    if verbose:
        print("initializing diffusion map computation, data is a " + str(data.shape[0]) + "X" + str(data.shape[1]) + " matrix.")
    A,distMatrix,epsilon = initialize_dmaps(data,epsilon_scale,n_evecs,LB_flag,data_type,estimate_eps,state_dependent_eps,exponent,cutoff_num,verbose,cutoff_dist)
    
    t = time()
    if verbose:
        print("starting diffusion map computation.")
    
    N = distMatrix.shape[0]

    D = np.array(A.sum(axis=1)).flatten()
    if((D == 0).any()):
        print('graph has disconnected subgraphs, using 1e15 distance')
        D = np.clip(D,1e-15,1e15)

    alpha = 1
    density = D

    D_inv = scipy.sparse.diags( D**(-alpha) )

    if (LB_flag==0): # FP norm
        M = D_inv.dot(A)
    else:            # LB norm
        A = D_inv.dot(A.dot(D_inv))

        D = np.array(A.sum(axis=1)).flatten()
        if((D == 0).any()):
            print('graph has disconnected subgraphs, using 1e15 distance')
            D = np.clip(D,1e-15,1e15)

        D_inv = scipy.sparse.diags( D**(-1/2) )

        M = D_inv.dot(A.dot(D_inv))

        M = M - 0.5 * (M-M.transpose())

    evals,evecs = scipy.sparse.linalg.eigsh(M,n_evecs)

    assert (evals.imag == 0).all()
    assert (evecs.imag == 0).all()

    if (LB_flag>0): # LB norm
        # transform back
        evecs = D_inv.dot(evecs)

    ix = np.argsort(np.abs(evals))
    ix = ix[::-1]
    evals = evals[ix]
    evecs = evecs[:,ix]
    
    for i in range(evecs.shape[1]):
        evecs[:,i] = evecs[:,i] - np.min(evecs[:,i])
        evecs[:,i] = (evecs[:,i] / np.max(evecs[:,i])) * 2 - 1
    
    if verbose:
        print('computed eigensystem in',int(time()-t),'seconds.')

    return evecs,evals,density,M,distMatrix,epsilon

####################################
def local_linear_regression(y, X, eps_med_scale):
    #
    # Code from Carmeline DSilva, adapted from MATLAB
    # 

    n = X.shape[0];

    K = scipy.spatial.distance.pdist(X, 'euclidean')
    K = scipy.spatial.distance.squareform(K);
    
    epssqr = np.median(K.flatten()**2)/eps_med_scale;
    W = np.exp(-K**2 / epssqr);

    L = np.zeros((n,n));
    for i in range(n):
        
        Xrep = np.matlib.repmat((X[i,:]), n, 1);
        Xones = np.ones((X.shape[0],1));
        
        Xx = np.hstack([Xones, X-Xrep]);
        
    #     Wx = diag(W(i,:));
    #     A = (Xx'*Wx*Xx)\(Xx'*Wx);
    
        # elementwise multiplication here:
        Xx2 = (Xx.transpose() * np.matlib.repmat(W[i,:], Xx.shape[1], 1));
        
        A,res,rank,s = np.linalg.lstsq(np.dot(Xx2,Xx), Xx2, rcond=None); # rcond=None added to support previous versions (numpy >=1.14.0)
        # print(A.shape)
        L[i,:] = A[0,:];

    fx = np.dot(L,y);
    
    stdy = np.std(y);
    omL = 1-np.diag(L);
    if(stdy == 0):
        stdy = 1;
        # print('error: stdy = 0');
    if((omL == 0).any()):
        omL[omL == 0] = 1;
        # print('error: omL = 0');
    
    res = np.sqrt(np.mean(((y-fx)/(omL))**2)) / np.std(y);
    return [fx, res];

def compute_residuals_DMAPS(V, eps_med_scale):
    #
    # Code from Carmeline DSilva, adapted from MATLAB
    # 
    # Computes the local linear regression error for each of the DMAPS
    # eigenvectors as a function of the previous eigenvectors
    # V are the DMAPS eigenvectors, stored in columns and ordered
    # V(:,1) is assumed to be the trivial constant eigenvector
    # eps_med_scale is the scale to use in the local linear regression kernel 
    # the kernel will be a Gaussian with width median(distances)/eps_med_scale
    # I typically take eps_med_scale = 3
    # res are the residuals of each of the fitted functions
    # res(0) is to be ignored, and res(1) will always be 1
    # res(i) is large/close to 1 if V(:,i) parameterizes a new direction in the
    # data, and res(i) is close to 0 if V(:,i) is a harmonic of a previous
    # eigenvector

    n = V.shape[1];
    res = np.zeros((n,));
    res[1] = 1;

    lt = time()
    for i in range(2,n):
        _,r = local_linear_regression(V[:,i], V[:, 1:i], eps_med_scale);
        res[i] = r
        if((time()-lt)>5):
            print('currently at i=',str(i),'of',str(n))
            lt = time()
    return res
    
    
def fwbw_dmaps(data,epsilon_scale,n_evecs,LB_flag,data_type="points",estimate_eps=True,state_dependent_eps=0,exponent=2,cutoff_num=0,verbose=False,cutoff_dist=-1):
    if verbose:
        print("initializing fw bw diffusion map computation, data is a " + str(data.shape[0]) + "X" + str(data.shape[1]) + " matrix.")
    A,distMatrix,epsilon = initialize_dmaps(data,epsilon_scale,n_evecs,LB_flag,data_type,estimate_eps,state_dependent_eps,exponent,cutoff_num,verbose,cutoff_dist)
    
    if verbose:
        print("starting fw bw diffusion map computation.")
    
    N = distMatrix.shape[0]

    q = 1./np.array(np.clip(A.sum(axis=1),1e-15,1e20)).flatten()
    assert ~((q == 0).any()), 'graph must not have disconnected points'
    q = scipy.sparse.diags(q)

    density = q

    Peps = q.dot(A)
    deps = scipy.sparse.diags(1./np.array(np.clip(Peps.sum(axis=0),1e-15,1e20)).flatten())
    
    B = deps.dot(Peps.transpose()).dot(Peps)
    if verbose:
        print("done.")
    
    return B

from joblib import Parallel, delayed
import multiprocessing

def dynamic_Laplacian(data_all,num_variables,dt,epsilon_scale,n_evecs,LB_flag,data_type="points",estimate_eps=True,state_dependent_eps=0,exponent=2,cutoff_num=0,verbose=False,cutoff_dist=-1):
    num_cores = multiprocessing.cpu_count()
    print('Performing ',str(data_all.shape[1]//2),' tasks in parallel.')
    B_all = Parallel(n_jobs=num_cores-1,verbose=9)(delayed(fwbw_dmaps)(data_all[:,(i*num_variables):((i+1)*num_variables)],epsilon_scale,n_evecs,LB_flag,data_type,estimate_eps,state_dependent_eps,exponent,cutoff_num,verbose,cutoff_dist)
                                    for i in range(data_all.shape[1]//2))
    
    Q = scipy.sparse.csr_matrix(np.zeros((data_all.shape[0],data_all.shape[0])))
    for i in range(data_all.shape[1]//2):
        Q = Q + B_all[i]
    Q = scipy.sparse.csr_matrix(1/dt * Q)
    return Q
    
def estimate_epsilon(distmatrix,emin=-5,emax=8,eN=40,exponent=2,cutoff_num=0):
    (x,y,z)=scipy.sparse.find(distmatrix)
    return np.median(z**exponent)
    if 1==0:
        epsilons = np.array([np.power(10,k) for k in np.linspace(emin,emax,eN)]);
        Asums = np.zeros((len(epsilons),));
        for i in range(len(epsilons)):
            epsilon = epsilons[i];
            A = np.exp(-distmatrix**exponent/(4*epsilon));
            Asum = np.median(A.flatten());
            Asums[i] = Asum;
        center = (np.max(Asums)-np.min(Asums)) / 2;
        ASMax = np.argsort(np.abs(Asums - center));
        return epsilons[ASMax[0]], epsilons, Asums;
