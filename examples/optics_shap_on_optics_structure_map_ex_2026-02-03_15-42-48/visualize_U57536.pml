load U57536_structure_heatmap.pdb
hide everything
show cartoon
spectrum b, blue_white_red, minimum=0, maximum=100
ramp_new scale, none, [0, 50, 100], [blue, white, red]
bg_color white
# NOTE: B-factors now represent ML Feature Importance (0=Low, 100=High)
