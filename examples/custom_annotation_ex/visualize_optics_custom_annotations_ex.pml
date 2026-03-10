# OPTICS Structure Annotation Script
load 1U19.pdb
hide everything
bg_color white
show cartoon, chain A
color white, chain A
set cartoon_transparency, 0.2, chain A

# --- Annotations ---
select site_124, chain A and resi 124
color red, site_124
show spheres, site_124
set sphere_scale, 1.2, site_124
label site_124 and n. ca, 'M124T (Missense)'
select site_296, chain A and resi 296
color yellow, site_296
show sticks, site_296
label site_296 and n. ca, 'Retinal Site'
select site_10, chain A and resi 10
color blue, site_10
set cartoon_transparency, 0.0, site_10
label site_10 and n. ca, 'N-Terminus'
select site_11, chain A and resi 11
color blue, site_11
set cartoon_transparency, 0.0, site_11
select site_12, chain A and resi 12
color blue, site_12
set cartoon_transparency, 0.0, site_12

# --- Final Settings ---
deselect
orient
ray 1000, 1000
