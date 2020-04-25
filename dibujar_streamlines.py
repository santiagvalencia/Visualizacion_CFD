import numpy as np
import time
import vtk
import pyvista as pv

pv.set_plot_theme("document")
geometria = pv.PolyData('cylinder.stl')
t0 = time.time()
grid = pv.read('internal.vtu') # foamToVTK
print('> VTK file read in '+str(time.time()-t0)+' seconds')

def lineas_verticales(x, z, y_sup, puntos_por_linea, y_inf = 0.000001):
    print('> '+str(x.size*z.size*puntos_por_linea)+' lines to be drawn')
    mesh = pv.PolyData()
    c = 0
    t_inicial = time.time()
    for x0 in x:
        for z0 in z:
            t0 = time.time()
            c+=1
            pa = (x0, y_inf, z0)
            pb = (x0, y_sup, z0)
            streamlines = grid.streamlines(vectors = 'U', n_points = puntos_por_linea,
            integrator_type = 45, initial_step_length = 0.0001, step_unit = 'l',
            max_steps = 10000, max_error = 1e-10, terminal_speed=1e-18, pointa = pa, pointb = pb,
            min_step_length=0.000001, max_step_length=0.01, interpolator_type = 'point', integration_direction = 'both')
            mesh = mesh + streamlines
            print('> '+str(c*puntos_por_linea)+' lines done, last '+str(puntos_por_linea)+' done in '+str(np.around(time.time()-t0, decimals = 2))+' seconds')
    print('> '+str(x.size*z.size*puntos_por_linea)+' lines done in '+str(np.around(time.time()-t_inicial, decimals = 2))+' seconds')
    return mesh

def punto(centro):
    streamlines = grid.streamlines(vectors = 'U', source_radius = 0, source_center=centro, n_points = 1,
    integrator_type = 45, initial_step_length = 0.0001, step_unit = 'l',
    max_steps = 100000, max_error = 1e-10, terminal_speed=1e-18,
    min_step_length=0.000001, max_step_length=0.01, interpolator_type = 'point', integration_direction = 'both')
    return streamlines

# datos para mesh (líneas de horsheshoe vortex)
x = np.linspace(-0.076, -0.15, 4)
z = np.linspace(-0.5, 5, 11)#np.array([-0.0002, -0.0001, 0, 0.0001, 0.0002])
puntos_por_linea = 10
y_sup = 0.005
mesh = lineas_verticales(x, z, y_sup, puntos_por_linea)


plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars = 'U')
plotter.add_mesh(geometria)
plotter.add_axes()
plotter.show_grid(color = 'gray')
plotter.show()
