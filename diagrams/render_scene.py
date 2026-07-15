"""Blender render of the FWA deployment scene.
Pipeline: run diagrams/make_deployment_3d.m in MATLAB first (replays the
drop construction and exports plots/deployment_scene.csv + the MATLAB 3D
figure), then render this scene headlessly:
  /Applications/Blender.app/Contents/MacOS/Blender -b -P diagrams/render_scene.py
Colors match the MATLAB figure legend one-to-one."""
import bpy, csv, math, random

CSV = '/Users/sdatta/FWA-with-CF/plots/deployment_scene.csv'
OUT = '/Users/sdatta/FWA-with-CF/plots/deployment_3d_blender.png'
ZX = 6.0  # vertical exaggeration, matches the MATLAB figure
random.seed(7)

# clean scene
bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete()
for block in (bpy.data.meshes, bpy.data.materials, bpy.data.lights, bpy.data.cameras):
    for item in list(block): block.remove(item)

def mat(name, rgb, rough=0.7, metal=0.0):
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']
    b.inputs['Base Color'].default_value = (*rgb, 1)
    b.inputs['Roughness'].default_value = rough
    b.inputs['Metallic'].default_value = metal
    return m

M = {'ground': mat('ground', (0.55, 0.53, 0.47), 0.95),
     'inner':  mat('inner',  (0.62, 0.60, 0.55), 0.95),
     'veg':    mat('veg',    (0.28, 0.44, 0.22), 0.95),
     'tree':   mat('tree',   (0.10, 0.28, 0.10), 0.9),
     'apt':    mat('apt',    (0.12, 0.28, 0.60), 0.6),
     'home':   mat('home',   (0.72, 0.65, 0.45), 0.8),
     'mall':   mat('mall',   (0.48, 0.28, 0.50), 0.6),
     'office': mat('office', (0.82, 0.35, 0.08), 0.5),
     'cpe':    mat('cpe',    (0.03, 0.08, 0.35), 0.3, 0.4),
     'mast':   mat('mast',   (0.15, 0.15, 0.17), 0.4, 0.8),
     'car':    mat('car',    (0.70, 0.03, 0.06), 0.3, 0.2),
     'ped':    mat('ped',    (0.25, 0.55, 0.12), 0.6)}

def box(x, y, w, h, m, z0=0.0):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, z0 + h/2))
    o = bpy.context.object; o.scale = (w, w, h); o.data.materials.append(m)
    return o

def disc(r, z, m):
    bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=1.0, location=(0, 0, z - 0.5), vertices=96)
    o = bpy.context.object; o.data.materials.append(m)
    return o

# large neutral ground plane so the cells sit in a landscape, not a void
bpy.ops.mesh.primitive_plane_add(size=12000, location=(0, 0, -0.5))
bpy.context.object.data.materials.append(mat('earth', (0.42, 0.40, 0.34), 0.95))

rows = list(csv.DictReader(open(CSV)))
gnbs = [r for r in rows if r['type'] == 'gnb']
veg  = [r for r in rows if r['type'] == 'veg'][0]
r_in, r_out = float(veg['height']), float(veg['width'])  # homeRadius, deployRange
sites_raw = [(float(g['x']), float(g['y'])) for g in gnbs]
# rotate the whole scene so the inter-site axis aligns with world x: the two
# cells then sit side by side and fill the frame instead of on a diagonal
theta = math.atan2(sites_raw[1][1] - sites_raw[0][1], sites_raw[1][0] - sites_raw[0][0])
ct, st = math.cos(-theta), math.sin(-theta)
def rot(x, y):
    return (x * ct - y * st, x * st + y * ct)
sites = [rot(x, y) for (x, y) in sites_raw]

# ground / vegetation / inner residential discs per site (stacked thin cylinders)
for (cx, cy) in sites:
    for (r, z, mm) in [(r_out, 0.0, M['ground']), (r_out, 0.3, M['veg']), (r_in, 0.6, M['inner'])]:
        bpy.ops.mesh.primitive_cylinder_add(radius=r, depth=0.3, location=(cx, cy, z), vertices=96)
        o = bpy.context.object; o.data.materials.append(mm)

# buildings + rooftop CPEs
FLOORS = {'apt': (2, 4, 28), 'home': (1, 2, 12), 'mall': (1, 1, 30), 'office': (2, 3, 24)}
for r in rows:
    t = r['type']
    if t in FLOORS:
        lo, hi, w = FLOORS[t]
        h = ZX * 3 * random.randint(lo, hi)
        x, y = rot(float(r['x']), float(r['y']))
        box(x, y, w, h, M[t])
        box(x, y, 7, 3, M['cpe'], z0=h)

# trees in the vegetation annulus of each site
for (cx, cy) in sites:
    for _ in range(220):
        rr = math.sqrt(r_in**2 + (r_out**2 - r_in**2) * random.random())
        aa = 2 * math.pi * random.random()
        tx, ty = cx + rr * math.cos(aa), cy + rr * math.sin(aa)
        hT = ZX * (7 + 5 * random.random())
        bpy.ops.mesh.primitive_cone_add(radius1=8, depth=hT, location=(tx, ty, hT/2), vertices=8)
        bpy.context.object.data.materials.append(M['tree'])

# gNB masts + panel head
for (cx, cy) in sites:
    hM = ZX * 35
    bpy.ops.mesh.primitive_cylinder_add(radius=3, depth=hM, location=(cx, cy, hM/2), vertices=12)
    bpy.context.object.data.materials.append(M['mast'])
    box(cx, cy, 16, 10, M['mast'], z0=hM)

# active devices
for r in rows:
    x, y = rot(float(r['x']), float(r['y']))
    if r['type'] == 'car':
        box(x, y, 11, ZX * 1.6, M['car'])
    elif r['type'] in ('ped', 'indoor'):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=5, location=(x, y, 6), segments=12, ring_count=8)
        bpy.context.object.data.materials.append(M['ped'])

# sun + sky
bpy.ops.object.light_add(type='SUN', location=(0, 0, 800))
sun = bpy.context.object; sun.data.energy = 5.5
sun.rotation_euler = (math.radians(50), 0, math.radians(-35))
bpy.data.worlds['World'].node_tree.nodes['Background'].inputs[0].default_value = (0.72, 0.80, 0.92, 1)
bpy.data.worlds['World'].node_tree.nodes['Background'].inputs[1].default_value = 0.55

# camera: tighter, more top-down oblique, centred on the two side-by-side cells
bpy.ops.object.empty_add(location=(0, 0, 18)); target = bpy.context.object
az, el, dist = math.radians(-15), math.radians(35), 3900
cx = dist * math.cos(el) * math.sin(az); cy = -dist * math.cos(el) * math.cos(az); cz = target.location[2] + dist * math.sin(el)
bpy.ops.object.camera_add(location=(cx, cy, cz))
cam = bpy.context.object; cam.data.clip_end = 20000; cam.data.lens = 48
tr = cam.constraints.new('TRACK_TO'); tr.target = target
bpy.context.scene.camera = cam

# render (wide aspect matched to the two side-by-side cells to minimise empty frame)
sc = bpy.context.scene
sc.render.engine = 'CYCLES'
sc.cycles.samples = 96
sc.cycles.use_denoising = True
sc.render.resolution_x, sc.render.resolution_y = 1800, 820
sc.render.filepath = OUT
bpy.ops.render.render(write_still=True)
print('RENDER DONE:', OUT)
