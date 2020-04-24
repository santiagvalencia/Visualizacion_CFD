import numpy as np
import time
import vtk
import pyvista as pv

#hi pepino

pv.set_plot_theme("document")
geometria = pv.PolyData('../Visualización/cylinder.stl')
t0 = time.time()
grid = pv.read('internal.vtu')
print('> VTK file read in '+str(time.time()-t0)+' seconds')

def lineas_verticales(x, z, y_sup, puntos_por_linea):
    print('> '+str(x.size*z.size*puntos_por_linea)+' lines to be drawn')
    mesh = pv.PolyData()
    c = 0
    t_inicial = time.time()
    for x0 in x:
        for z0 in z:
            t0 = time.time()
            c+=1
            pa = (x0, 0.000001, z0)
            pb = (x0, y_sup, z0)
            streamlines = grid.streamlines(vectors = 'U', n_points = puntos_por_linea,
            integrator_type = 45, initial_step_length = 0.0001, step_unit = 'l',
            max_steps = 10000, max_error = 1e-10, terminal_speed=1e-18, pointa = pa, pointb = pb,
            min_step_length=0.000001, max_step_length=0.01, interpolator_type = 'point')
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

x = np.linspace(-0.076, -0.15, 3)
z = np.array([0])#np.array([-0.0002, -0.0001, 0, 0.0001, 0.0002])#np.linspace(-0.02, 0.02, 5)
puntos_por_linea = 4
y_sup = 0.005

linea = punto((-0.11, 0.006, 0))
# list = dir(linea)
# for l in list:
#     print(l)
#mesh = lineas_verticales(x, z, y_sup, puntos_por_linea)
vort = linea['vorticity']
# vel = linea['U']
points = linea.points
#
# vorticidad = pv.PolyData(points)
# vorticidad.vectors = vort
#
# linea.set_active_scalars('vorticity')
# linea
#linea.set_active_vectors('vorticity')
#print(linea._active_vectors_info)
glyphs = linea.glyph(orient = True, scale = True, geom = pv.Arrow(), factor = 1)
#
plotter = pv.Plotter()
plotter.add_arrows(points[::5], vort[::5], mag = 0.01)
# #plotter.add_arrows(points[::5], vel[::5], mag = 10, color = 'black')
# # #plotter.add_mesh(mesh, color = 'gray')
plotter.add_mesh(geometria)
# # #plotter.add_mesh(linea.tube(radius = 0.001), scalars = 'vorticity')
# plotter.add_mesh(glyphs)
# # plotter.add_axes()
# # plotter.show_grid(color = 'gray')
plotter.show()
