import numpy as np
import time
import vtk
import pyvista as pv
import sys

######################### VTK & STL files #########################
pv.set_plot_theme("document")
geometria = pv.PolyData(sys.argv[2]) # Archivo .stl
t0 = time.time()
grid = pv.read(sys.argv[1]) # Introducir nombre de archivo vtk en linea de comando (foamToVTK)
print('> VTK file read in '+str(time.time()-t0)+' seconds')

######################### Inicializa Plot #########################
plotter = pv.Plotter()
plotter.add_mesh(geometria)
plotter.add_axes()
plotter.show_grid(color = 'gray')

######################### Funciones core ##########################
def punto(centro, field = 'U', direction = 'both', radio = 0.0, num = 1, max_steps = 100000, terminal_speed = 1e-18, initial_step_length = 0.0001, integrator_type = 45, min_step_length=0.000001):
    streamlines = grid.streamlines(vectors = field, source_radius = radio, source_center=centro, n_points = num,
    integrator_type = integrator_type, initial_step_length = initial_step_length, step_unit = 'l',
    max_steps = max_steps, max_error = 1e-10, terminal_speed=terminal_speed,
    min_step_length=min_step_length, max_step_length=0.01, interpolator_type = 'point', integration_direction = direction)
    return streamlines

def lineas_verticales(x, z, y_sup, puntos_por_linea, y_inf = 0.000001, field = 'U', direction = 'both'):
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
            streamlines = grid.streamlines(vectors = field, n_points = puntos_por_linea,
            integrator_type = 45, initial_step_length = 0.0001, step_unit = 'l',
            max_steps = 10000, max_error = 1e-10, terminal_speed=1e-18, pointa = pa, pointb = pb,
            min_step_length=0.000001, max_step_length=0.01, interpolator_type = 'point', integration_direction = direction)
            mesh = mesh + streamlines
            print('> '+str(c*puntos_por_linea)+' lines done, last '+str(puntos_por_linea)+' done in '+str(np.around(time.time()-t0, decimals = 2))+' seconds')
    print('> '+str(x.size*z.size*puntos_por_linea)+' lines done in '+str(np.around(time.time()-t_inicial, decimals = 2))+' seconds')
    return mesh

def helicity(mesh):
    helicity = np.zeros((mesh['U'].shape[0], 1))
    for i in range(mesh['U'].shape[0]):
        u = mesh['U'][i]
        w = mesh['vorticity'][i]
        helicity[i] = u.dot(w)/(np.linalg.norm(u)*np.linalg.norm(w))
    return helicity

def vectores(mesh, vort_scale = 0.005, vel_scale = 0.3, tol = 0.01, absolute = False):
    mesh.set_active_vectors('vorticity')
    mesh.point_arrays['mag_vort'] = np.linalg.norm(mesh.active_vectors, axis = 1)
    mesh.point_arrays['vec_vort'] = mesh.active_vectors
    arrows_vort = mesh.glyph(scale = 'mag_vort', orient = 'vec_vort', tolerance = tol, factor = vort_scale, absolute = absolute)
    mesh.set_active_vectors('U')
    mesh.point_arrays['mag_vel'] = np.linalg.norm(mesh.active_vectors, axis = 1)
    mesh.point_arrays['vec_vel'] = mesh.active_vectors
    arrows_vel = mesh.glyph(scale = 'mag_vel', orient = 'vec_vel', tolerance = tol, factor = vel_scale, absolute = absolute)
    return arrows_vort, arrows_vel

def coordinate_input():
    coord_x = float(input('Enter x-coordinate:  '))
    coord_y = float(input('Enter y-coordinate:  '))
    coord_z = float(input('Enter z-coordinate:  '))
    return (coord_x,coord_y,coord_z)

############################ Callables ############################
def plano_x(x, normal = 'x', field = 'U'):
    origin = (x, 0, 0)
    slc = grid.slice(normal = normal, origin = origin)
    w, u = vectores(slc, vel_scale = 1)
    if field == 'U':
        plotter.add_mesh(u, name = 'U', scalars = 'U')
    else:
        plotter.add_mesh(w, name = 'w', scalars = 'vorticity')

def plano_y(y, normal = 'y', field = 'U'):
    origin = (0, y, 0)
    slc = grid.slice(normal = normal, origin = origin)
    w, u = vectores(slc, vel_scale = 1)
    if field == 'U':
        plotter.add_mesh(u, name = 'U', scalars = 'U')
    else:
        plotter.add_mesh(w, name = 'w', scalars = 'vorticity')

def plano_z(z, normal = 'z', field = 'U'):
    origin = (0, 0, z)
    slc = grid.slice(normal = normal, origin = origin)
    w, u = vectores(slc, vel_scale = 1)
    if field == 'U':
        plotter.add_mesh(u, name = 'U', scalars = 'U')
    else:
        plotter.add_mesh(w, name = 'w', scalars = 'vorticity')

def linea_interactive(pointa, pointb):
    global plotter
    streamlines = grid.streamlines(vectors = 'U', n_points=30, pointa=pointa, pointb=pointb,
                                   integrator_type = 45, initial_step_length = 0.0001, step_unit = 'l',
                                   max_steps = 100000, max_error = 1e-10, terminal_speed=1e-18,min_step_length=0.000001, max_step_length=0.01, interpolator_type = 'point', integration_direction='forward')
    plotter.add_mesh(streamlines, name='streamlines', scalars = 'vorticity')

def punto_interactive(centro, centro2, length_forward=8, length_backward=8):
    global plotter, actor_text
    print('Centro vel: '+str(centro))
    print('Centro vort: '+str(centro2))
    streamlines = punto(centro, radio = 0.001, num = 10)
    streamlinesA = punto(centro2, field = 'vorticity', direction = 'forward', terminal_speed = 1e-20, max_steps = int(1e9), radio = 0.001, num = 10)
    streamlinesB = punto(centro2, field = 'vorticity', direction = 'backward', terminal_speed = 1e-20, max_steps = int(1e9), radio = 0.001, num = 10)
    mesh_S = streamlines + streamlinesA + streamlinesB
    centroA_old = (streamlinesA.points[-1][0],streamlinesA.points[-1][1],streamlinesA.points[-1][2])
    centroB_old = (streamlinesB.points[-1][0],streamlinesB.points[-1][1],streamlinesB.points[-1][2])
    for i in range(length_forward):
        streamlinesA_new = punto(centroA_old, field = 'vorticity', direction = 'forward', terminal_speed = 1e-20, max_steps = int(1e9), radio = 0.001, num = 10)
        centroA_old = (streamlinesA_new.points[-1][0],streamlinesA_new.points[-1][1],streamlinesA_new.points[-1][2])
        mesh_S = mesh_S + streamlinesA_new
    for j in range(length_backward):
        streamlinesB_new = punto(centroB_old, field = 'vorticity', direction = 'backward', terminal_speed = 1e-20, max_steps = int(1e9), radio = 0.001, num = 10)
        centroB_old = (streamlinesB_new.points[-1][0],streamlinesB_new.points[-1][1],streamlinesB_new.points[-1][2])
        mesh_S = mesh_S + streamlinesB_new
    plotter.add_mesh(mesh_S, name = 'streamlines', scalars = 'vorticity')
    plotter.remove_actor(actor_text)
    actor_text = plotter.add_text('Streamline: ' + str(centro) + '        Vorticity: ' + str(centro2))

def coordinate_interactive_x(value):
    global plotter, point, point_a, point_actor
    point_a=value
    point[0]=value
    plotter.remove_actor([point_actor])
    p = pv.Sphere(radius=0.01, center=point)
    point_actor = plotter.add_mesh(p, 'red')

def coordinate_interactive_y(value):
    global plotter, point, point_b, point_actor
    point_b=value
    point[1]=value
    plotter.remove_actor([point_actor])
    p = pv.Sphere(radius=0.01, center=point)
    point_actor = plotter.add_mesh(p, 'red')

def coordinate_interactive_z(value):
    global plotter, point, point_c, point_actor
    point_c=value
    point[2]=value
    plotter.remove_actor([point_actor])
    p = pv.Sphere(radius=0.01, center=point)
    point_actor = plotter.add_mesh(p, 'red')

def refresh(active):
    global plotter, point_a, point_b, point_c
    if active==True:
        streamline = punto((point_a, point_b, point_c), radio=0.001, num = 20, direction = 'both')
        plotter.add_mesh(streamline, name = 'streamlines', scalars = 'vorticity')

def reset_input_point(active):
    global plotter, point_actor
    if active==True:
        centro = coordinate_input()
        #mesh = punto(centro)
        mesh=punto(centro, radio=0.01, num = 50)
        plotter.remove_actor([point_actor])
        plotter.add_mesh(mesh, name = 'streamlines', scalars = 'vorticity')
        point_actor = plotter.add_mesh(pv.Sphere(radius=0.01, center=centro), color = 'red')

def lineas_verticales_interactive(z_c):
    global plotter, mesh_actor
    plotter.remove_actor(mesh_actor)
    x = np.linspace(-0.35, -0.30, 10)
    z = np.linspace(z_c-0.0001, z_c+0.0001, 3)
    y_sup = 0.1
    puntos_por_linea = 8
    mesh1 = lineas_verticales(x, z, y_sup, puntos_por_linea)
    mesh_actor = plotter.add_mesh(mesh1, scalars = 'vorticity')


###################### Funciones multimedia #######################
def animacion(centro, field = 'U', frames  = 1000, dt = 0.0001, max_steps = 20):
    mesh = punto(centro, field = field, initial_step_length = dt, integrator_type = 4, max_steps = max_steps, direction = 'forward')
    mesh.point_arrays['helicity'] = helicity(mesh)
    w, u = vectores(mesh, tol = 1.5*mesh.length, absolute = True)
    plotter.add_mesh(mesh, scalars = 'helicity')
    text_actor = plotter.add_text('h = '+str(mesh.point_arrays['helicity'][-1]), position = 'upper_right')
    u_actor = plotter.add_mesh(u, color = 'blue')
    w_actor = plotter.add_mesh(w, color = 'red')
    plotter.show(auto_close = False)
    plotter.open_gif('animacion.gif')
    for i in range(frames):
        plotter.remove_actor([u_actor, w_actor, text_actor])
        centro = mesh.points[-1]
        mesh_new = punto((centro[0], centro[1], centro[2]), initial_step_length = dt, integrator_type = 4, max_steps = max_steps, direction = 'forward')
        mesh_new.point_arrays['helicity'] = helicity(mesh_new)
        w, u = vectores(mesh_new, tol = 1.5*mesh_new.length, absolute = True)
        mesh+= mesh_new
        plotter.add_mesh(mesh, scalars = 'helicity')
        u_actor = plotter.add_mesh(u, color = 'blue')
        w_actor = plotter.add_mesh(w, color = 'red')
        text_actor = plotter.add_text('h = '+str(mesh_new.point_arrays['helicity'][-1]), position = 'upper_right')
        plotter.write_frame()
    plotter.close()

########################## Ejecutables ############################

########### Input punto ###########
# mesh=punto(coordinate_input())
# plotter.add_mesh(mesh, name = 'streamlines', scalars = 'vorticity')
# plotter.show()

########### Input punto interactive ###########
# point_actor = plotter.add_mesh(pv.Sphere(radius=0.01, center= (0, 0.5, 0)), color = 'red')
# plotter.add_checkbox_button_widget(reset_input_point, position=(200, 100), color_on='r')
# plotter.show()

########### Puntos (U y vorticity) interactive ###########
# actor_text = plotter.add_text('Streamline: ' + str((0, 0, 0)) + '        Vorticity: ' + str((0, 0, 0)))
# plotter.add_line_widget(callback=punto_interactive, use_vertices=True)
#plotter.show()

########### Linea interactive ###########
# plotter.add_line_widget(callback=linea_interactive, use_vertices=True)
# plotter.show()

########### Slider punto ###########
# point_a=0.5
# point_b=0.5
# point_c=0.5
# point=np.array((point_a,point_b,point_c))
# p = pv.Sphere(radius=0.01, center=point)
# point_actor = plotter.add_mesh(p, 'red')
# plotter.add_slider_widget(coordinate_interactive_x, [-0.8,0.8], value=point_a, title='Eje x', pointa=(0.75,0.9), pointb=(0.95,0.9))
# plotter.add_slider_widget(coordinate_interactive_y, [0,1], value=point_b, title='Eje y', pointa=(0.75,0.7), pointb=(0.95,0.7))
# plotter.add_slider_widget(coordinate_interactive_z, [-0.5,0.5], value=point_c, title='Eje z', pointa=(0.75,0.5), pointb=(0.95,0.5))
# plotter.add_checkbox_button_widget(refresh, position=(100, 100))
#plotter.show()

########### Slider plano ###########
# plotter.add_slider_widget(plano_x, [-1, 1], title = 'x')
# plotter.show()

########### Animación ###########
# centro = (-0.3899709, 0.15468751, 0.00157578)
# animacion(centro)


########################## Plots ###########################
plotter.show()
