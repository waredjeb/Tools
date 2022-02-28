import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)
def plots3D(data, dataE, data_idx, seeds, save = ''):
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]
    sizes = [100*(i/max(dataE)) for i in dataE]
    ax = fig.add_subplot(1,1,1, projection = '3d')   
    ax.scatter(x, y, z, s = sizes, marker='o', alpha=1, c = data_idx) 
    # ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'red')
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    plt.savefig(save + ".png", bbox_inches='tight')


def plots3DwithProjectionSeeds(data, dataE, data_idx, seeds, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]

    sizes = [100*(i/max(dataE)) for i in dataE]
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    # cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(x, y, z, s = sizes, marker='o', alpha=.3, c = data_idx)
    ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'red')
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    ax.scatter(x, y, s = sizes, c = data_idx, alpha = .5)
    ax.scatter(sx, sy, s = 100, marker='*', color = 'red')
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    ax.scatter(x, z, s = sizes, c = data_idx, alpha = .5)
    ax.scatter(sx, sz, s = 100, marker='*', color = 'red')
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    im = ax.scatter(y, z, s = sizes, c = data_idx, alpha = .5)
    ax.scatter(sy, sz, s = 100, marker='*', color = 'red')

    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs Energy')
    # cbar.set_label("LCs energy", fontsize = 50)
    plt.savefig(save + ".png", bbox_inches='tight')
    # close(fig)
def plots3DwithProjectionFineTracksters3D(data, dataE, data_idx, data3D, dataE3D, data_idx3D, seeds, save = ''):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]

    x3D = [i.x() for i in data3D]
    y3D = [i.y() for i in data3D]
    z3D = [i.z() for i in data3D]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]

    # sizes = [100*(i/max(dataE)) for i in dataE]
    # sizes3D = [100*(i/max(dataE3D)) for i in dataE3D]
    sizes = [100 for i in dataE]
    sizes3D= [100+10 for i in dataE3D]
    # fig  = plt.figure()
    ax = fig.add_subplot(2,1,1, projection = '3d')    
    # cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(x, y, z, s = sizes, marker='o', alpha=0.8, c = data_idx , label = 'SimTracksters')
    # ax.scatter(x3D, y3D, z3D, s = sizes3D, marker='*', alpha=1, c = data_idx3D, label = 'FineSimTrackster')
    # ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'red', label = 'CLUE3DSeeds')
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    # cm = plt.cm.get_cmap('hot')
    noise_idx_l = -1
    ax.scatter(x3D, y3D, z3D, s = sizes3D, marker='o', alpha=0.8, c = data_idx3D , label = 'FineSimTracksters')
    # col = ['r','g','b']
    # colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    # ax = fig.add_subplot(2,2,2)
    # ax.scatter(x, y, s = sizes, alpha = 0.8, color = 'red')
    # ax.scatter(x3D, y3D, s = sizes3D, marker='*', alpha=1, c = data_idx3D)# c = data_idx3D)
    # # ax.scatter(sx, sy, s = 50, marker='*', color = 'red')
    # ax.set_xlabel('$X$')
    # ax.set_ylabel('$Y$')
    # ax.set_title("X-Y")


    # ax = fig.add_subplot(2,2,3  ) 
    # ax.scatter(x, z, s = sizes, alpha = 0.8 , color = 'red')
    # ax.scatter(x3D, z3D, s = sizes3D, marker='*', alpha= 1, c = data_idx3D)
    # # ax.scatter(sx, sz, s = 50, marker='*', color = 'red')
    
    # ax.set_xlabel('$X$')
    # ax.set_ylabel('$Z$')
    # ax.set_title("X-Z")

    
    # ax = fig.add_subplot(2,2,4)
    # im = ax.scatter(y, z, s = sizes, alpha = 0.8, color = 'red')
    # ax.scatter(y3D, z3D, s = sizes3D, marker='*', alpha=1     , c = data_idx3D)
    # # ax.scatter(sy, sz, s = 50, marker='*', color = 'red')

    # ax.set_xlabel('$Y$')
    # ax.set_ylabel('$Z$')
    # ax.set_title("Y-Z")

    # # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # # cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs Energy')
    # # cbar.set_label("LCs energy", fontsize = 50)
    plt.savefig(save + ".png", bbox_inches='tight')
    # close(fig)
def plots3DwithProjectionSeeds3D(data, dataE, data_idx, data3D, dataE3D, data_idx3D, seeds, xmin,xmax,ymin,ymax,zmin,zmax, title, save = '', map_ = False):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]

    x3D = [i.x() for i in data3D]
    y3D = [i.y() for i in data3D]
    z3D = [i.z() for i in data3D]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]

    # sizes = [100*(i/max(dataE)) for i in dataE]
    # sizes3D = [100*(i/max(dataE3D)) for i in dataE3D]
    sizes = [100 for i in dataE]
    sizes3D= [100+10 for i in dataE3D]
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('tab10')
    noise_idx_l = -1
    if(map_ == True):
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=1, c = data_idx, cmap= cm)
    else:
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=1, c = data_idx)
    # ax.scatter(x3D, y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = data_idx3D, label = 'CARecoTrackster')
    # ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'red', label = 'CLUE3DSeeds')
    plt.title(title)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin,zmax)
    # plt.legend()
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    if(map_ == True):
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 1)
    # ax.scatter(x3D, y3D, s = sizes3D, marker='x', alpha=0.5, c = data_idx3D)
    # ax.scatter(sx, sy, s = 50, marker='*', color = 'red')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # ax.set_zlim(zlim,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    if(map_ == True):
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 1)
    # ax.scatter(x3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = data_idx3D)
    # ax.scatter(sx, sz, s = 50, marker='*', color = 'red')
    ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    if(map_ == True):
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 1)
    # ax.scatter(y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = data_idx3D)
    # ax.scatter(sy, sz, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    ax.set_xlim(ymin, ymax)
    ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs Energy')
    # cbar.set_label("LCs energy", fontsize = 50)
    plt.savefig(save + ".png", bbox_inches='tight')
    # close(fig)


def plots3DwithProjectionSeeds3D_1(data, dataE, data_idx, data3D, dataE3D, data_idx3D, seeds, title, save = '', map_ = False):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]

    x3D = [i.x() for i in data3D]
    y3D = [i.y() for i in data3D]
    z3D = [i.z() for i in data3D]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]

    # sizes = [100*(i/max(dataE)) for i in dataE]
    # sizes3D = [100*(i/max(dataE3D)) for i in dataE3D]
    sizes = [100 for i in dataE]
    sizes3D= [100+10 for i in dataE3D]
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('tab10')
    noise_idx_l = -1
    if(map_ == True):
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=1, c = data_idx, cmap= cm)
    else:
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=0.8, c = data_idx, label = "SlimSimTracksters")
        ax.scatter(x3D, y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red', label = 'Missing LCs')
    # ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'red', label = 'CLUE3DSeeds')
    plt.title(title)
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_zlim(zmin,zmax)
    plt.legend()
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    if(map_ == True):
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(x3D, y3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
    # ax.scatter(sx, sy, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_zlim(zlim,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    if(map_ == True):
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(x3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
    # ax.scatter(sx, sz, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    if(map_ == True):
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
    # ax.scatter(sy, sz, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_xlim(ymin, ymax)
    # ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs Energy')
    # cbar.set_label("LCs energy", fontsize = 50)
    plt.savefig(save + ".png", bbox_inches='tight')
    # close(fig)

def plots3DwithProjectionSeeds3D_2(data, dataE, data_idx, data3D, dataE3D, data_idx3D, seeds, seedsE, title, save = '', map_ = False):
    '''
    Function to plot 3D representation with all the projections on the three planes
    - data: np.ndarray(4, num_rows). .iloc[:, :3] = coordinates, .iloc[:,3] = energy
    '''
    fig = plt.figure(figsize = (30,25))
    x = [i.x() for i in data]
    y = [i.y() for i in data]
    z = [i.z() for i in data]

    x3D = [i.x() for i in data3D]
    y3D = [i.y() for i in data3D]
    z3D = [i.z() for i in data3D]
    sx = [i.x() for i in seeds]
    sy = [i.y() for i in seeds]
    sz = [i.z() for i in seeds]
    sE = [e for e in seedsE]

    # sizes = [100*(i/max(dataE)) for i in dataE]
    # sizes3D = [100*(i/max(dataE3D)) for i in dataE3D]
    sizes = [100 for i in dataE]
    sizes3D= [100+10 for i in dataE3D]
    # fig  = plt.figure()
    ax = fig.add_subplot(2,2,1, projection = '3d')    
    cm = plt.cm.get_cmap('tab10')
    noise_idx_l = -1
    if(map_ == True):
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=1, c = data_idx, cmap= cm)
    else:
        ax.scatter(x, y, z, s = sizes, marker='o', alpha=0.8, c = data_idx, label = "SlimSimTracksters")
        ax.scatter(x3D, y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red', label = 'Missing LCs')
        ax.scatter(sx, sy, sz, s = 10, marker='*', color = 'black', label = 'Slim Seeds')
    plt.title(title)
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_zlim(zmin,zmax)
    plt.legend()
    col = ['r','g','b']
    colors = col
    ax.set_xlabel('$X$', rotation=150)
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$',rotation=60)
    ax = fig.add_subplot(2,2,2)
    if(map_ == True):
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, y, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(x3D, y3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
        ax.scatter(sx, sy, s = 10, marker='*', color = 'black', label = 'Slim Seeds')
    # ax.scatter(sx, sy, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_zlim(zlim,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_title("X-Y")


    ax = fig.add_subplot(2,2,3)
    if(map_ == True):
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(x, z, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(x3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
        ax.scatter(sx, sz, s = 10, marker='*', color = 'black', label = 'Slim Seeds')
    # ax.scatter(sx, sz, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_ylim(ymin, ymax)
    # ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Z$')
    ax.set_title("X-Z")

    
    ax = fig.add_subplot(2,2,4)
    if(map_ == True):
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 1,cmap= cm)
    else:
        ax.scatter(y, z, s = sizes, c = data_idx, alpha = 0.8)
        ax.scatter(y3D, z3D, s = sizes3D, marker='x', alpha=0.5, c = 'red')
        ax.scatter(sy, sz, s = 10, marker='*', color = 'black', label = 'Slim Seeds')
    # ax.scatter(sy, sz, s = 50, marker='*', color = 'red')
    # ax.set_xlim(xmin, xmax)
    # ax.set_xlim(ymin, ymax)
    # ax.set_ylim(zmin,zmax)
    ax.set_xlabel('$Y$')
    ax.set_ylabel('$Z$')
    ax.set_title("Y-Z")

    # cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(im, cax=cbar_ax, label = 'LCs Energy')
    # cbar.set_label("LCs energy", fontsize = 50)
    plt.savefig(save + ".png", bbox_inches='tight')
    # close(fig)
import os

#Create new directory
def mkdir_p(mypath):
    '''Function to create a new directory, if it not already exist
        - mypath : directory path
    '''
    from errno import EEXIST
    from os import makedirs,path
    try:
        makedirs(mypath)
    except OSError as exc:
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
