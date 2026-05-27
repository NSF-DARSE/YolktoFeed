def rc_update(rcParams, w= 3.375, h = 1.35, fs = 8, lw = 0.8, 
              ax_lw = 0.5, tick_M = 1.5, tick_m = 0.8, tex = 0,
              minor_x = 0, minor_y = 0, np_h=1, np_v=1, minor = 0,
              mks = 1, mkfs='full', mew=0.1):
    '''
    Parameters
    -------
    w,h         : width and height of the figure
    fs          : font size of labels and legend
    lw          : line width
    ax_lw       : axis line width
    tick_M/m    : minor tick size
    tex         : use latex boolean
    minor       : minor visible boolean
    Returns
    -------
    Updated rcParams
    '''
    #rcParams.update({'figure.figsize': [w*np_h,h*np_v]})
    #rcParams.update({'figure.dpi': 700}) 

    if minor: 
        minor_x = 1
        minor_y = 1
    if tex:
        rcParams.update({'text.usetex': True}) 

    #rcParams.update({'xtick.labelsize': fs})
    #rcParams.update({'ytick.labelsize': fs}) 
    #rcParams.update({'axes.labelsize' : fs})
    #rcParams.update({'legend.fontsize': fs})
    rcParams.update({'lines.markersize': mks})
    rcParams.update({'markers.fillstyle': mkfs})

    rcParams.update({'xtick.direction': 'in'}) 
    rcParams.update({'ytick.direction': 'in'}) 
    rcParams.update({'xtick.major.width': ax_lw}) 
    rcParams.update({'ytick.major.width': ax_lw}) 

    rcParams.update({'xtick.major.pad': 3})
    rcParams.update({'ytick.major.pad': 3})
    if minor_x:
        rcParams.update({'xtick.minor.visible':True})
        rcParams.update({'xtick.minor.size': tick_m}) 
        rcParams.update({'xtick.minor.width': ax_lw}) 
    if minor_y:
        rcParams.update({'ytick.minor.visible':True})  
        rcParams.update({'ytick.minor.size': tick_m}) 
        rcParams.update({'ytick.minor.width': ax_lw}) 
     
    rcParams.update({'axes.linewidth':  ax_lw})
    rcParams.update({'lines.linewidth': lw})
    rcParams.update({'lines.markeredgewidth': mew})
    
    rcParams.update({'xtick.major.size': tick_M})
    rcParams.update({'ytick.major.size': tick_M})  
    
    
    #rcParams.update({'xtick.top': False}) 
    #rcParams.update({'xtick.labeltop': False}) 
    #rcParams.update({'xtick.bottom': True})  
    #rcParams.update({'xtick.labelbottom': True})

    return 0

