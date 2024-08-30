import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from numpy import pi, angle, ndarray
import plotly.graph_objs as go
import matplotlib.colors as mcolors
import matplotlib as mpl


"""
Function to get extent of x,y (RA/DEC) coordinate system and to get extent of u,v coordinate system
"""
def get_extent_coord_system(header):
    convert_to_milliarc = 180/pi * 3600 * 1000
    xy_extent = [convert_to_milliarc*header.grid_RA[0,0],convert_to_milliarc*header.grid_RA[0,-1],convert_to_milliarc*header.grid_DEC[-1,0],convert_to_milliarc*header.grid_DEC[0,0]]
    
    uv_extent = [header.u[0,-1],header.u[0,0],header.v[0,0],header.v[-1,0]]

    return xy_extent, uv_extent

"""
Function to get dictionaries used in plotting for GUI
"""
def get_plot_dictionaries_gui(header, org, mdl2, anl2, anlDerivative,residuals):

    radius_earth = 6371000
    
    scale_factor = 1e12
    anlDerivative = anlDerivative*scale_factor      #Converting unit to ps

    xy_extent, uv_extent = get_extent_coord_system(header)

    plot_dict_mdl = {'data': [org,mdl2,org-mdl2],
                    'titles': ['Source Image', 'Model of Source Image','Residual'],
                    'circle': [False, False, False],
                    'x_axis': ['Relative RA (milliarcsec)', 'Relative RA (milliarcsec)','Relative RA (milliarcsec)'],
                    'y_axis':  ['Relative DEC (milliarcsec)','Relative DEC (milliarcsec)', 'Relative DEC (milliarcsec)'],
                    'cmap' : ['viridis','viridis','viridis'],
                    'extent': [xy_extent,xy_extent,xy_extent],
                    'circle_dict': {'x' : header.reference_pixel_RA, 'y' : header.reference_pixel_DEC, 'r': radius_earth}}
        
    plot_dict_anl = {'data': [org, mdl2, abs(anl2),angle(anl2), anlDerivative, residuals],
                'titles': ['Original', 'Model', 'Visibility Analytical Model', 'Phase Analytical Model (rad)','Group Delay Analytical Model','Residual Group Delay (ps)'],
                'circle': [False, False, True, True, True, True],
                'x_axis': ['Relative RA (milliarcsec)','Relative RA (milliarcsec)', 'U','U', 'U', 'U'],
                'y_axis':  ['Relative DEC (milliarcsec)','Relative DEC (milliarcsec)', 'V','V', 'V','V'],
                'cmap': ['viridis','viridis','viridis','viridis','viridis','viridis'],
                'extent': [xy_extent,xy_extent,uv_extent,uv_extent,uv_extent,uv_extent],
                'circle_dict': {'x' : header.reference_pixel_RA, 'y' : header.reference_pixel_DEC, 'r': radius_earth}}
    
    return plot_dict_mdl, plot_dict_anl


"""
Function to get dictionaries used in plotting for GUI
"""
def get_plot_dictionaries_gui_clean(header, org, mdl2, clean_visibility, dclean_visibility):

    radius_earth = 6371000

    xy_extent, uv_extent = get_extent_coord_system(header)

    cmap_custom = (mpl.colors.ListedColormap(['red', 'orange', 'yellow', 'lime','green','deepskyblue','blue','blueviolet','dimgrey','lightgrey'])
                    .with_extremes(over='lightgrey', under='red'))
        
    plot_dict_clean = {'data': [org, mdl2, abs(clean_visibility),angle(clean_visibility), dclean_visibility],
                'titles': ['Original', 'Model', 'Visibility Clean components', 'Phase Clean components (rad)','Group Delay Clean components'],
                'circle': [False, False, True, True, True, True],
                'x_axis': ['Relative RA (milliarcsec)','Relative RA (milliarcsec)', 'U','U', 'U'],
                'y_axis':  ['Relative DEC (milliarcsec)','Relative DEC (milliarcsec)', 'V','V', 'V'],
                'cmap': ['viridis','viridis',cmap_custom,cmap_custom,cmap_custom],
                'extent': [xy_extent,xy_extent,uv_extent,uv_extent,uv_extent],
                'circle_dict': {'x' : header.reference_pixel_RA, 'y' : header.reference_pixel_DEC, 'r': radius_earth}}
    
    return plot_dict_clean


"""
Function to get dictionaries used in plotting when running code using -r
"""
def get_plot_dictionaries(header, org, mdl2, anl2, anlDerivative,residuals):

    radius_earth = 6371000
    
    scale_factor = 1e12
    anlDerivative = anlDerivative*scale_factor      #Converting unit to ps

    xy_extent, uv_extent = get_extent_coord_system(header)
    
    plot_dict = {'data': [org, mdl2,org-mdl2, abs(anl2), angle(anl2), anlDerivative, residuals],
                'titles': ['Source Image fits file', 'Model','Residual', 'Visibility Analytical Model', 'Phase Analytical Model (rad)','Group Delay Analytical Model','Residual Group Delay (ps)'],
                'circle': [False, False, False, True, True, True, True],
                'x_axis': ['Relative RA (milliarcsec)','Relative RA (milliarcsec)','Relative RA (milliarcsec)', 'U','U', 'U', 'U'],
                'y_axis':  ['Relative DEC (milliarcsec)','Relative DEC (milliarcsec)','Relative DEC (milliarcsec)', 'V','V', 'V','V'],
                'cmap': ['viridis','viridis','viridis','viridis','viridis','viridis','viridis'],
                'extent': [xy_extent,xy_extent,xy_extent,uv_extent,uv_extent,uv_extent,uv_extent],
                'circle_dict': {'x' : header.reference_pixel_RA, 'y' : header.reference_pixel_DEC, 'r': radius_earth}}

    return plot_dict

"""
Function to plot dictionaries when running code using -r
"""
def plot_results(plot_dict, desired_min_value = None, desired_max_value = None):

    if not plot_dict or not plot_dict['data']:
        return

    num_cols = 4                           # Fixed number of columns (2)
    num_rows = 2                           # Number of rows based on total images

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(25, 10))
    axes = axes.flatten() if isinstance(axes, ndarray) else [axes]

    # Iterate over the data and titles in plot_dict to plot on the axes
    for i, (data, title,extent,cmap, circle,x_axis, y_axis) in enumerate(zip(plot_dict['data'], plot_dict['titles'],plot_dict['extent'],
                                                                plot_dict['cmap'], plot_dict['circle'],plot_dict['x_axis'],plot_dict['y_axis'])):
        ax = axes[i]
        norm = mcolors.Normalize(vmin=desired_min_value, vmax=desired_max_value)
        im = ax.imshow(data, cmap=cmap,extent = extent, origin='lower', norm=norm) # Plot the image
    
        if circle:
            circ = Circle((plot_dict['circle_dict']['x'], plot_dict['circle_dict']['y']), plot_dict['circle_dict']['r'], facecolor='None', edgecolor='w', lw=1)
            ax.add_patch(circ)
        
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        ax.set_title(title)
        fig.colorbar(im, ax=ax)

    for j in range(i + 1, len(axes)):
            axes[j].axis('off')

    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.3)
    plt.show()


def plot_two_surfaces_3d(u,v,data1,data2):

    scale_factor = 10e12
    data1 = data1*scale_factor
    data2 = data2*scale_factor

    custom_colorscale = [[0, 'rgba(0, 0, 0, 0.8)'],[1, 'rgba(0, 0, 0, 0.8)']]

    surface1 = go.Surface(x=u,y=v,z=data1, opacity=0.7,colorscale='Viridis',showscale=True)
    surface2 = go.Surface(x=u,y=v,z=data2,opacity=0.5,colorscale=custom_colorscale,showscale=False)

    fig = go.Figure(data=[surface1, surface2])

    fig.update_layout(scene=dict(xaxis_title='U',yaxis_title='V',zaxis_title='Picoseconds'),
                      title='3D Surface of group delay and fitted plane',width=1000,height=600)

    # Save the figure as an interactive HTML file
    fig.write_html('3d_plot_surface.html')
    fig.show()


def plot_3d(u,v,data,range_data, file_name = '3d_plot.html'):

    surface1 = go.Surface(x=u,y=v,z=data,opacity=0.7, colorscale='Viridis',showscale=True)

    fig = go.Figure(data=[surface1])

    fig.update_layout(scene=dict(xaxis_title='U',yaxis_title='V',zaxis_title='Picoseconds',zaxis=dict(range=range_data)),
        title='3D Surface of residual group delay',width=1000,height=600)

    # Save the figure as an interactive HTML file
    fig.write_html(file_name)
    fig.show()